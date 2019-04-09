namespace libiop {

template<typename FieldT>
bcs16_prover<FieldT>::bcs16_prover(const bcs16_transformation_parameters<FieldT> &parameters) :
    bcs16_protocol<FieldT>(parameters)
{
}

template<typename FieldT>
void bcs16_prover<FieldT>::signal_prover_round_done()
{
    enter_block("Finish prover round");
    iop_protocol<FieldT>::signal_prover_round_done();
    std::size_t ended_round = this->num_prover_rounds_done_-1;
    const domain_to_oracles_map mapping = this->oracles_in_round(ended_round);
    const round_parameters<FieldT> round_params = this->get_round_parameters(ended_round);

    /* First, go through all the oracle messages in this round and
       compress each one using a Merkle Tree. */
    hash_digest cur_state = (this->pseudorandom_state_.empty() ? "" : *this->pseudorandom_state_.rbegin());
    std::vector<hash_digest> roots;
    for (auto &kv : mapping)
    {
        std::vector<std::vector<FieldT> > all_evaluated_contents;
        for (auto &v : kv.second)
        {
            std::vector<FieldT> oracle_contents = this->oracles_[v.id()].evaluated_contents();
            this->serialize_leaf_data_by_round_params(oracle_contents,
                                                      all_evaluated_contents,
                                                      kv.first,
                                                      round_params);
        }
        enter_block("Construct Merkle tree");
        this->Merkle_trees_[this->MTs_processed_].construct(all_evaluated_contents);
        leave_block("Construct Merkle tree");
        cur_state = this->parameters_.compression_hasher(cur_state,
                                                         this->Merkle_trees_[this->MTs_processed_].get_root(),
                                                         this->digest_len_bytes_);

        ++this->MTs_processed_;
    }

    /* Hash explicitly sent prover messages */
    const std::size_t min_message_id =
        (this->num_prover_rounds_done_ == 1 ? 0 : this->num_prover_messages_at_end_of_round_[ended_round - 1]);
    const std::size_t max_message_id = this->num_prover_messages_at_end_of_round_[ended_round];

    std::vector<FieldT> message_concat = { FieldT::zero() };
    for (std::size_t message_id = min_message_id; message_id < max_message_id; ++message_id)
    {
        message_concat.insert(message_concat.end(),
                              this->prover_messages_[message_id].begin(),
                              this->prover_messages_[message_id].end());

    }
#ifdef DEBUG
    printf("Message concat (min_message_id=%zu, max_message_id=%zu, num_prover_rounds_done-1=%zu:\n",
           min_message_id,
           max_message_id,
           ended_round);
    for (auto &v : message_concat)
    {
        v.print();
    }
#endif // DEBUG

    const hash_digest message_hash = this->parameters_.field_hasher(message_concat, this->digest_len_bytes_);
    cur_state = this->parameters_.compression_hasher(cur_state, message_hash, this->digest_len_bytes_);

    /* Add the prover message hash as a "root" and update the pseudorandom state */
    this->pseudorandom_state_.emplace_back(cur_state);

    leave_block("Finish prover round");
}

/* Serializes the provided oracle's evaluated contents into multiple rows,
   and appends these rows to all_oracles_evaluated_contents.
   It is serialized into multiple rows according to the quotient map relation in round params.
   The default round parameters result in oracle_evaluated_contents being appended as a single row. */
template<typename FieldT>
void bcs16_prover<FieldT>::serialize_leaf_data_by_round_params(
    std::vector<FieldT> &oracle_evaluated_contents,
    std::vector<std::vector<FieldT>> &all_evaluated_contents,
    const domain_handle &evaluation_domain_handle,
    const round_parameters<FieldT> &round_params)
{
    if (round_params.quotient_map_size_ == 1) {
        all_evaluated_contents.emplace_back(oracle_evaluated_contents);
        return;
    }
    std::size_t start_row = all_evaluated_contents.size();
    /** Add a vector initialized to 0 to all_evaluated_contents,
     *  quotient_map_size times. */
    std::size_t num_cosets = oracle_evaluated_contents.size() / round_params.quotient_map_size_;
    for (std::size_t i = 0; i < round_params.quotient_map_size_; i++) {
        std::vector<FieldT> vector_for_ith_element_of_each_coset(num_cosets);
        all_evaluated_contents.emplace_back(vector_for_ith_element_of_each_coset);
    }
    /** Each leaf will now contain the evaluations of every element of the coset being queried.
     *  This is done such that within a leaf's column, the elements of the coset appear in order. */
    field_subset<FieldT> evaluation_domain = this->get_domain(evaluation_domain_handle);
    if (evaluation_domain.type() != round_params.quotient_map_type_)
    {
        throw std::invalid_argument("round params' domain type does not match the oracle evaluation's domain.");
    }
    if (evaluation_domain.type() == affine_subspace_type) {
        /** This assumes that the quotient map domain's basis vectors are the first
         *  log_2(|Q|) basis vectors of the evaluation domain.
         *
         *  When this is the case, elements of the same coset are consecutive
         *  elements of the evaluation domain.
         *  This follows from how we index subspaces.
         */
        std::size_t coset_index = 0; /* Which element within a single coset are we considering */
        for (std::size_t i = 0; i < oracle_evaluated_contents.size(); i++) {
            all_evaluated_contents[start_row + coset_index][i / round_params.quotient_map_size_]
                = oracle_evaluated_contents[i];
            coset_index = (coset_index + 1) % round_params.quotient_map_size_;
        }
    } else if (evaluation_domain.type() == multiplicative_coset_type) {
        /** In the multiplicative setting,
         *  Let i be the index of element x \in evaluation_domain.
         *  Then [x] has elements with indices i, i + |H|/|Q|, i + 2|H|/|Q|, ... i + (|Q| - 1) |H|/|Q|.
         *  These are then grouped into the same leaf.
         *
         *  For cache efficiency, this is done by iterating through every 1st element of a coset, then 2nd element, etc.
         */
        for (std::size_t i = 0; i < round_params.quotient_map_size_; i++)
        {
            // TODO: Replace this loop with std::move or just direct copying for better memory efficiency
            for (std::size_t j = 0; j < num_cosets; j++)
            {
                all_evaluated_contents[start_row + i][j] =
                    oracle_evaluated_contents[i*num_cosets + j];
            }
        }
    } else {
        throw std::invalid_argument("BCS16 IOP - Unknown field_subset type for quotient map domain");
    }
}

/* Each "random" verifier message is deterministically constructed from the previous
   one, so that both prover and verifier can reconstruct. */
template<typename FieldT>
std::vector<FieldT> bcs16_prover<FieldT>::obtain_verifier_random_message(
    const verifier_random_message_handle &random_message)
{
    /* TODO: technical debt */
    iop_protocol<FieldT>::obtain_verifier_random_message(random_message);

    const std::size_t message_length = this->verifier_random_message_registrations_[random_message.id()].size();
    /* Find the index of the round containing this message. */
    const std::size_t round = (std::lower_bound(this->num_verifier_random_messages_at_end_of_round_.begin(),
                                                this->num_verifier_random_messages_at_end_of_round_.end(),
                                                random_message.id() + 1) -
                               this->num_verifier_random_messages_at_end_of_round_.begin());
#ifdef DEBUG
    printf("prover: random message id=%zu, round=%zu\n", random_message.id(), round);
#endif

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
    for (auto &v : result)
    {
        v.print();
    }
#endif // DEBUG

    return result;
}

template<typename FieldT>
FieldT bcs16_prover<FieldT>::obtain_query_response(const query_handle &query)
{
    /* Just defer to oracles. */
    return iop_protocol<FieldT>::obtain_query_response(query);
}

template<typename FieldT>
bcs16_transformation_transcript<FieldT> bcs16_prover<FieldT>::get_transcript()
{
    bcs16_transformation_transcript<FieldT> result;

    /*
      Easy part: fill in (explicit) prover messages and MT roots.
    */
    result.prover_messages_ = this->prover_messages_;

    for (auto &MT : this->Merkle_trees_)
    {
        result.MT_roots_.emplace_back(MT.get_root());
    }

    /*
      Make sure all queries are hit.
    */
    for (std::size_t query_id = 0; query_id < this->query_registrations_.size(); ++query_id)
    {
        this->obtain_query_response(query_handle(query_id));
        /* this will populate oracle_id_to_query_positions_ (a hack used
           below in constructing the transcript to collect query positions) */
    }

    /*
      Go over all MTs and prepare multi membership proofs. (Now that we
      know both what the oracle messages are and where they'll be queried,
      we can find authentication paths.)
    */
    result.total_depth_without_pruning = 0;
    std::size_t MT_idx = 0;
    for (std::size_t round = 0; round < this->num_interaction_rounds_; ++round)
    {
        const domain_to_oracles_map mapping = this->oracles_in_round(round);
        const round_parameters<FieldT> round_params = this->get_round_parameters(round);

        /* Go over each oracle message in this round. */
        for (auto &kv : mapping)
        {
            /* Where is this oracle going to be queried? */
            std::set<std::size_t> query_positions_set;
            std::set<std::size_t> MT_leaf_positions_set;
            const std::size_t num_leaves = this->domains_[kv.first.id()].num_elements() / round_params.quotient_map_size_;
            for (auto oracle_h : kv.second)
            {
                for (auto pos : this->oracle_id_to_query_positions_[oracle_h.id()])
                {
                    query_positions_set.insert(pos);
                    MT_leaf_positions_set.insert(
                        query_position_to_merkle_tree_position(pos, num_leaves, round_params));
#ifdef DEBUG
                    printf("record: oracle %zu at position %zu\n", oracle_h.id(), pos);
#endif // DEBUG
                }
            }

            std::vector<std::size_t> query_positions(query_positions_set.begin(), query_positions_set.end());
            std::vector<std::size_t> MT_leaf_positions(MT_leaf_positions_set.begin(), MT_leaf_positions_set.end());

            /* What are the values at those query positions? */
            std::vector<std::vector<FieldT> > values;
            for (auto pos : query_positions)
            {
                std::vector<FieldT> column;
                for (auto oracle_h : kv.second)
                {
                    column.emplace_back(this->get_oracle_evaluation_at_point(std::make_shared<oracle_handle>(oracle_h), pos));
                }
                values.emplace_back(column);
            }

            result.total_depth_without_pruning += MT_leaf_positions.size() * this->Merkle_trees_[MT_idx].depth();

            result.query_positions_.emplace_back(query_positions);
            result.MT_leaf_positions_.emplace_back(MT_leaf_positions);
            result.query_responses_.emplace_back(values);

            /* Generate a combined authentication path for all those queries. */
            const merkle_tree_multi_membership_proof proof =
                this->Merkle_trees_[MT_idx].get_multi_membership_proof(MT_leaf_positions);

            result.MT_multi_membership_proofs_.emplace_back(proof);

            ++MT_idx;
        }
    }

    return result;
}

template<typename FieldT>
std::size_t bcs16_prover<FieldT>::MT_size() const
{
    std::size_t MT_size = 0;
    for (auto &MT : this->Merkle_trees_)
    {
        MT_size += MT.num_total_bytes();
    }

    return MT_size;
}

template<typename FieldT>
std::size_t bcs16_prover<FieldT>::state_size() const
{
    return (this->num_bytes_across_all_oracles() + this->MT_size());
}

template<typename FieldT>
void bcs16_prover<FieldT>::describe_sizes() const
{
    print_indent(); printf("* Total size of proof oracles (bytes): %zu\n", this->num_bytes_across_all_oracles());
    print_indent(); printf("* Total size of Merkle tree (bytes): %zu\n", this->MT_size());
    print_indent(); printf("* Total size of prover state (bytes): %zu\n", this->state_size());
}

} // namespace libiop
