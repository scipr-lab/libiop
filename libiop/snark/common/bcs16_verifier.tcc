namespace libiop {

template<typename FieldT>
bcs16_verifier<FieldT>::bcs16_verifier(const bcs16_transformation_parameters<FieldT> &parameters,
                                               const bcs16_transformation_transcript<FieldT> &transcript) :
    bcs16_protocol<FieldT>(parameters),
    transcript_(transcript)
{
}

template<typename FieldT>
void bcs16_verifier<FieldT>::seal_interaction_registrations()
{
    enter_block("seal_interaction_registrations");
    bcs16_protocol<FieldT>::seal_interaction_registrations();

    /* Compute pseudorandom state (chaining together the MT roots) */
    hash_digest cur_state = "";

    std::size_t MTs_processed = 0;
    std::vector<size_t> MT_index_to_round;
    for (std::size_t round = 0; round < this->num_interaction_rounds_; ++round)
    {
        /* Update the pseudorandom state for the ORACLE messages. */
        const domain_to_oracles_map mapping = this->oracles_in_round(round);

        for (auto &kv : mapping)
        {
            MT_index_to_round.emplace_back(round);
            /* We aren't using the oracles at this point, but we know that
               one MT root corresponds to each oracle in this round. */
            UNUSED(kv);
            cur_state = this->parameters_.compression_hasher(cur_state,
                                                             this->transcript_.MT_roots_[MTs_processed++],
                                                             this->digest_len_bytes_);
        }

        /* Concatenate and then hash together EXPLICIT messages. */
        const std::size_t min_message_id = (round == 0 ? 0 : this->num_prover_messages_at_end_of_round_[round-1]);
        const std::size_t max_message_id = this->num_prover_messages_at_end_of_round_[round];

        std::vector<FieldT> message_concat = { FieldT::zero() };
        for (std::size_t message_id = min_message_id; message_id < max_message_id; ++message_id)
        {
            message_concat.insert(message_concat.end(),
                                  this->transcript_.prover_messages_[message_id].begin(),
                                  this->transcript_.prover_messages_[message_id].end());
        }

#ifdef DEBUG
        printf("Message concat (min_message_id=%zu, max_message_id=%zu, round=%zu:\n",
               min_message_id,
               max_message_id,
               round);
        for (auto &v : message_concat)
        {
            v.print();
        }
#endif // DEBUG

        /* Add the prover message hash as a "root" and update the pseudorandom state */
        const hash_digest message_hash = this->parameters_.field_hasher(message_concat, this->digest_len_bytes_);
        cur_state = this->parameters_.compression_hasher(cur_state, message_hash, this->digest_len_bytes_);

        this->pseudorandom_state_.emplace_back(cur_state);
    }

    /* Validate all MT queries relative to the transcript. */
    this->transcript_is_valid_ = true;

    for (std::size_t MT_idx = 0; MT_idx < this->transcript_.MT_roots_.size(); ++MT_idx)
    {
        auto &root = this->transcript_.MT_roots_[MT_idx];
        std::vector<std::size_t> &query_positions = this->transcript_.query_positions_[MT_idx];
        std::vector<std::size_t> &MT_leaf_positions = this->transcript_.MT_leaf_positions_[MT_idx];
        std::vector<std::vector<FieldT> > &query_responses = this->transcript_.query_responses_[MT_idx];
        auto &proof = this->transcript_.MT_multi_membership_proofs_[MT_idx];

        // Step 1) Serialize query responses into leaf responses
        std::vector<std::vector< FieldT> > MT_leaf_columns =
            this->query_responses_to_MT_leaf_responses(query_positions,
                                                       query_responses,
                                                       MT_index_to_round[MT_idx],
                                                       this->Merkle_trees_[MT_idx].num_leaves());
        // Step 2) hash each leaf column
        std::vector<hash_digest> column_hashes;
        for (auto &v : MT_leaf_columns)
        {
            column_hashes.emplace_back(this->parameters_.field_hasher(
                                       v,
                                       this->digest_len_bytes_));
        }
        // Step 3) validate proof

        const bool proof_is_valid = this->Merkle_trees_[MT_idx].validate_multi_membership_proof(
            root,
            MT_leaf_positions,
            column_hashes,
            proof);

        if (!proof_is_valid)
        {
            this->transcript_is_valid_ = false;
        }
    }

    /* Finally populate things for obtaining query responses */
    this->parse_query_responses_from_transcript();
    leave_block("seal_interaction_registrations");
}

template<typename FieldT>
void bcs16_verifier<FieldT>::parse_query_responses_from_transcript()
{
    std::size_t MTs_processed = 0;
    for (std::size_t round = 0; round < this->num_interaction_rounds_; ++round)
    {
        const domain_to_oracles_map mapping = this->oracles_in_round(round);
        const round_parameters<FieldT> round_params = this->get_round_parameters(round);

        /* For each oracle, pick out positions and values from the transcript
            and store them. */
        for (auto &kv : mapping)
        {
            std::size_t oracles_processed_for_MT = 0;
            for (auto &oh : kv.second)
            {
                for (std::size_t i = 0; i < this->transcript_.query_positions_[MTs_processed].size(); ++i)
                {
                    const std::size_t pos_idx = transcript_.query_positions_[MTs_processed][i];
                    const FieldT value = transcript_.query_responses_[MTs_processed][i][oracles_processed_for_MT];

                    this->oracle_id_and_pos_idx_to_value_[std::make_pair(oh.id(), pos_idx)] = value;
                }
                ++oracles_processed_for_MT;
            }

            ++MTs_processed;
        }
    }
}

template<typename FieldT>
std::vector<std::vector<FieldT> > bcs16_verifier<FieldT>::query_responses_to_MT_leaf_responses(
    std::vector<size_t> &query_positions,
    std::vector<std::vector<FieldT> > &query_responses,
    const size_t round,
    const size_t num_leaves)
{
    const round_parameters<FieldT> round_params = this->get_round_parameters(round);
    if (round_params.quotient_map_size_ == 1)
    {
        return query_responses;
    }
    // Initialize all the columns
    std::vector<std::vector<FieldT>> MT_leaf_columns(query_positions.size() / round_params.quotient_map_size_);
    const size_t leaf_size = query_responses[0].size() * round_params.quotient_map_size_;
    for (size_t i = 0; i < MT_leaf_columns.size(); i++)
    {
        MT_leaf_columns[i] = std::vector<FieldT>(leaf_size);
    }
    /** Elements within a given coset appear in order,
     * so we simply store the index for the next element of the coset,
     * and increment as we see new positions belonging to this coset. */
    std::vector<size_t> intra_coset_index(MT_leaf_columns.size(), 0);
    std::map<size_t, size_t> MT_leaf_pos_to_response_index;
    size_t next_response_index = 0;
    for (size_t i = 0; i < query_positions.size(); i++)
    {
        const size_t query_position = query_positions[i];
        const size_t MT_leaf_index = query_position_to_merkle_tree_position(query_position, num_leaves, round_params);
        std::map<size_t, size_t>::iterator it = MT_leaf_pos_to_response_index.find(MT_leaf_index);
        /* For supported domain types, new MT leaf positions appear in order of query positions.
         * If we don't yet know the index of this leaf within the queried for leaves,
         * we can find it by simply incrementing the prior leaf's index. */
        if (it == MT_leaf_pos_to_response_index.end()) {
            MT_leaf_pos_to_response_index[MT_leaf_index] = next_response_index;
            next_response_index++;
        }
        const size_t MT_response_index = MT_leaf_pos_to_response_index[MT_leaf_index];
        const size_t index_in_coset = intra_coset_index[MT_response_index];
        intra_coset_index[MT_response_index]++;
        for (size_t j = 0; j < query_responses[i].size(); j++)
        {
            const size_t oracle_index = j*round_params.quotient_map_size_;
            MT_leaf_columns[MT_response_index][oracle_index + index_in_coset] =
                query_responses[i][j];
        }
    }
    return MT_leaf_columns;
}

template<typename FieldT>
void bcs16_verifier<FieldT>::signal_prover_round_done()
{
    throw std::logic_error("Verifier IOP is not meant for proving.");
}

template<typename FieldT>
std::vector<FieldT> bcs16_verifier<FieldT>::obtain_verifier_random_message(
    const verifier_random_message_handle &random_message)
{
    const std::size_t message_length = this->verifier_random_message_registrations_[random_message.id()].size();

    /* Find the index of the round containing this message. */
    const std::size_t round = (std::lower_bound(this->num_verifier_random_messages_at_end_of_round_.begin(),
                                                this->num_verifier_random_messages_at_end_of_round_.end(),
                                                random_message.id() + 1) -
                               this->num_verifier_random_messages_at_end_of_round_.begin());

    /* Use the pseudorandom state from the PREVIOUS round (that's how the
       "random" message is generated, because the pseudorandom state is
       constructed at the END of each round.) */
    hash_digest prev_pseudorandom_state;
    if (round == 0)
    {
        prev_pseudorandom_state = "";
    }
    else
    {
        prev_pseudorandom_state = this->pseudorandom_state_[round - 1];
    }

    const std::vector<FieldT> result =
        this->parameters_.FieldT_randomness_extractor(prev_pseudorandom_state,
                                                      random_message.id(),
                                                      message_length);
    this->verifier_random_messages_[random_message.id()] = result;

#ifdef DEBUG
    printf("verifier: random message id=%zu, round=%zu\n", random_message.id(), round);
    for (auto &v : result)
    {
        v.print();
    }
#endif // DEBUG

    return result;
}

template<typename FieldT>
FieldT bcs16_verifier<FieldT>::get_oracle_evaluation_at_point(
    const oracle_handle_ptr &handle,
    const std::size_t evaluation_position,
    const bool record)
{
    UNUSED(record);

    if (std::dynamic_pointer_cast<virtual_oracle_handle>(handle))
    {
        /* If virtual oracle, we just defer to the original code (the
           get_oracle_evaluation_at_point function in iop_protocol),
           which call this function to get each of the constituent
           evaluations. */
        return bcs16_protocol<FieldT>::get_oracle_evaluation_at_point(handle, evaluation_position, false);
    }
    else if (std::dynamic_pointer_cast<oracle_handle>(handle))
    {
        /* If real oracle, use our saved values that we saved from
           the transcript. */
        auto it = this->oracle_id_and_pos_idx_to_value_.find(std::make_pair(handle->id(), evaluation_position));

#ifdef DEBUG
        printf("query: oracle %zu at position %zu\n", handle->id(), evaluation_position);
#endif // DEBUG

        if (it == this->oracle_id_and_pos_idx_to_value_.end())
        {
            throw std::logic_error("Got a request for a query position that's unavailable in the proof.");
        }
        else
        {
            return it->second;
        }

    }
    else
    {
        throw std::invalid_argument("oracle type not supported");
    }
}

template<typename FieldT>
std::vector<FieldT> bcs16_verifier<FieldT>::receive_prover_message(const prover_message_handle &message)
{
    return this->transcript_.prover_messages_[message.id()];
}

template<typename FieldT>
bool bcs16_verifier<FieldT>::transcript_is_valid() const
{
    return this->transcript_is_valid_;
}

} // namespace libiop
