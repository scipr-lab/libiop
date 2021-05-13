namespace libiop {

template<typename FieldT, typename MT_hash_type>
bcs_verifier<FieldT, MT_hash_type>::bcs_verifier(
    const bcs_transformation_parameters<FieldT, MT_hash_type> &parameters,
    const bcs_transformation_transcript<FieldT, MT_hash_type> &transcript) :
    bcs_protocol<FieldT, MT_hash_type>(parameters),
    transcript_(transcript),
    is_preprocessing_(false)
{
}

template<typename FieldT, typename MT_hash_type>
bcs_verifier<FieldT, MT_hash_type>::bcs_verifier(
    const bcs_transformation_parameters<FieldT, MT_hash_type> &parameters,
    const bcs_transformation_transcript<FieldT, MT_hash_type> &transcript,
    const bcs_verifier_index<FieldT, MT_hash_type> &index) :
    bcs_protocol<FieldT, MT_hash_type>(parameters),
    transcript_(transcript),
    is_preprocessing_(true),
    index_(index)
{
    /** We check that the indexer provided the correct number of roots and messages in
     *  seal_interaction_registrations()    */
    this->transcript_.MT_roots_.insert(
        this->transcript_.MT_roots_.begin(),
        index.index_MT_roots_.begin(),
        index.index_MT_roots_.end());
    this->transcript_.prover_messages_.insert(
        this->transcript_.prover_messages_.begin(),
        index.indexed_messages_.begin(),
        index.indexed_messages_.end());
}

template<typename FieldT, typename hash_digest_type>
void bcs_verifier<FieldT, hash_digest_type>::seal_interaction_registrations()
{
    libff::enter_block("verifier_seal_interaction_registrations");
    bcs_protocol<FieldT, hash_digest_type>::seal_interaction_registrations();

    this->transcript_is_valid_ = true;

    std::size_t processed_MTs = 0; // Updated at end of loop.
    for (std::size_t round = 0; round < this->num_interaction_rounds_; ++round)
    {
        /* Update the pseudorandom state for the oracle messages. */
        const std::size_t num_domains = this->num_domains_in_round(round);
        if (this->is_preprocessing_ && round == 0)
        {
            if (num_domains != this->index_.index_MT_roots_.size())
            {
                throw std::invalid_argument("Index had an incorrect number of MT roots");
            }
            if (this->num_prover_messages_at_end_of_round_[0]
                 != this->index_.indexed_messages_.size())
            {
                throw std::invalid_argument("Index had an incorrect number of prover messages");
            }
        }

        /* Absorb MT roots into hashchain. */
        // Each domain has one Merkle tree containing all the oracles.
        const std::vector<hash_digest_type> MT_roots_for_round(
            this->transcript_.MT_roots_.cbegin() + processed_MTs,
            this->transcript_.MT_roots_.cbegin() + processed_MTs + num_domains);
        this->run_hashchain_for_round(round, MT_roots_for_round, this->transcript_.prover_messages_);

        /* Validate all MT queries relative to the transcript. */
        for (std::size_t i = 0; i < num_domains; i++)
        {
            const auto &root = this->transcript_.MT_roots_[processed_MTs];
            std::vector<std::size_t> &query_positions = this->transcript_.query_positions_[processed_MTs];
            std::vector<std::size_t> &MT_leaf_positions = this->transcript_.MT_leaf_positions_[processed_MTs];
            std::vector<std::vector<FieldT> > &query_responses = this->transcript_.query_responses_[processed_MTs];
            const auto &proof = this->transcript_.MT_set_membership_proofs_[processed_MTs];

            // Step 1) serialize query responses into leafs
            std::vector<std::vector< FieldT> > MT_leaf_columns =
                this->query_responses_to_MT_leaf_responses(query_positions, query_responses, round);

            // Step 2) validate proof
            const bool proof_is_valid = this->Merkle_trees_[processed_MTs]
                .validate_set_membership_proof(root, MT_leaf_positions, MT_leaf_columns, proof);

            if (!proof_is_valid)
            {
                this->transcript_is_valid_ = false;
            }
            processed_MTs++;
        }
    }

    /* Check proof of work */

    hash_digest_type pow_challenge = this->hashchain_->squeeze_root_type();
    bool valid_pow = this->pow_.verify_pow(this->parameters_.compression_hasher, pow_challenge, this->transcript_.proof_of_work_);
    if (!valid_pow)
    {
        printf("Invalid pow\n");
        this->transcript_is_valid_ = false;
    }

    /* Finally populate things for obtaining query responses */
    this->parse_query_responses_from_transcript();
    libff::leave_block("verifier_seal_interaction_registrations");
}

template<typename FieldT, typename MT_hash_type>
void bcs_verifier<FieldT, MT_hash_type>::parse_query_responses_from_transcript()
{
    std::size_t processed_MTs = 0;
    for (std::size_t round = 0; round < this->num_interaction_rounds_; ++round)
    {
        const domain_to_oracles_map mapping = this->oracles_in_round_by_domain(round);
        const round_parameters<FieldT> round_params = this->get_round_parameters(round);

        /* For each oracle, pick out positions and values from the transcript
            and store them. */
        for (auto &kv : mapping)
        {
            std::size_t oracles_processed_for_MT = 0;
            for (auto &oh : kv.second)
            {
                for (std::size_t i = 0; i < this->transcript_.query_positions_[processed_MTs].size(); ++i)
                {
                    const std::size_t pos_idx = transcript_.query_positions_[processed_MTs][i];
                    const FieldT value = transcript_.query_responses_[processed_MTs][i][oracles_processed_for_MT];

                    this->oracle_id_and_pos_idx_to_value_[std::make_pair(oh.id(), pos_idx)] = value;
                }
                ++oracles_processed_for_MT;
            }

            ++processed_MTs;
        }
    }
}

template<typename FieldT, typename MT_hash_type>
std::vector<std::vector<FieldT> > bcs_verifier<FieldT, MT_hash_type>::query_responses_to_MT_leaf_responses(
    std::vector<size_t> &query_positions,
    std::vector<std::vector<FieldT> > &query_responses,
    const size_t round)
{
    const round_parameters<FieldT> round_params = this->get_round_parameters(round);
    if (round_params.quotient_map_size_ == 1)
    {
        return query_responses;
    }
    return this->Merkle_trees_[round].serialize_leaf_values_by_coset(
        query_positions,
        query_responses,
        round_params.quotient_map_size_
    );
}

template<typename FieldT, typename MT_hash_type>
void bcs_verifier<FieldT, MT_hash_type>::signal_prover_round_done()
{
    throw std::logic_error("Verifier IOP is not meant for proving.");
}

template<typename FieldT, typename MT_hash_type>
void bcs_verifier<FieldT, MT_hash_type>::signal_index_submissions_done()
{
    throw std::logic_error("Verifier IOP is not meant for indexing.");
}

template<typename FieldT, typename MT_hash_type>
std::vector<FieldT> bcs_verifier<FieldT, MT_hash_type>::obtain_verifier_random_message(
    const verifier_random_message_handle &random_message)
{
    return this->verifier_random_messages_[random_message.id()];
}

template<typename FieldT, typename MT_hash_type>
FieldT bcs_verifier<FieldT, MT_hash_type>::get_oracle_evaluation_at_point(
    const oracle_handle_ptr &handle,
    const std::size_t evaluation_position,
    const bool record)
{
    libff::UNUSED(record);

    if (std::dynamic_pointer_cast<virtual_oracle_handle>(handle))
    {
        /* If virtual oracle, we just defer to the original code (the
           get_oracle_evaluation_at_point function in iop_protocol),
           which call this function to get each of the constituent
           evaluations. */
        return bcs_protocol<FieldT, MT_hash_type>::get_oracle_evaluation_at_point(handle, evaluation_position, false);
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
        return it->second;
    }
    else
    {
        throw std::invalid_argument("oracle type not supported");
    }
}

template<typename FieldT, typename MT_hash_type>
std::vector<FieldT> bcs_verifier<FieldT, MT_hash_type>::receive_prover_message(const prover_message_handle &message)
{
    return this->transcript_.prover_messages_[message.id()];
}

template<typename FieldT, typename MT_hash_type>
bool bcs_verifier<FieldT, MT_hash_type>::transcript_is_valid() const
{
    return this->transcript_is_valid_;
}

} // namespace libiop
