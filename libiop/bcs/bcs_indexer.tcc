namespace libiop {

template<typename FieldT, typename MT_hash_type>
bcs_indexer<FieldT, MT_hash_type>::bcs_indexer(
    const bcs_transformation_parameters<FieldT, MT_hash_type> &parameters) :
    bcs_protocol<FieldT, MT_hash_type>(parameters)
{
}

template<typename FieldT, typename MT_hash_type>
void bcs_indexer<FieldT, MT_hash_type>::signal_prover_round_done()
{
    throw std::invalid_argument("Indexer should not be used for proving"
        " (signal prover round done called)");
}

template<typename FieldT, typename MT_hash_type>
void bcs_indexer<FieldT, MT_hash_type>::signal_index_submissions_done()
{
    libff::enter_block("Merkelize indexed oracles");
    iop_protocol<FieldT>::signal_prover_round_done();
    std::size_t ended_round = this->num_prover_rounds_done_-1;
    if (ended_round != 0)
    {
        throw std::invalid_argument("Index submissions should be round 0");
    }
    const domain_to_oracles_map mapping = this->oracles_in_round_by_domain(ended_round);
    const round_parameters<FieldT> round_params = this->get_round_parameters(ended_round);

    /* First, go through all the oracle messages in this round and
       compress each one using a Merkle Tree. */
    std::vector<MT_hash_type> roots;
    for (auto &kv : mapping)
    {
        std::vector<std::shared_ptr<std::vector<FieldT>>> all_evaluated_contents;
        for (auto &v : kv.second)
        {
            std::shared_ptr<std::vector<FieldT>> oracle_contents = this->oracles_[v.id()].evaluated_contents();
            all_evaluated_contents.emplace_back(oracle_contents);
        }
        libff::enter_block("Construct Merkle tree");
        this->Merkle_trees_[this->MTs_processed_].construct_with_leaves_serialized_by_cosets(
            all_evaluated_contents, round_params.quotient_map_size_);
        libff::leave_block("Construct Merkle tree");

        ++this->MTs_processed_;
        /* Now make the oracles in a form suitable for creating an index */
        for (auto &v : kv.second)
        {
            std::shared_ptr<std::vector<FieldT>> oracle_contents = this->oracles_[v.id()].evaluated_contents();
            this->indexed_oracles_.emplace_back(*oracle_contents.get());
            this->oracles_[v.id()].erase_contents();
        }
    }

    libff::leave_block("Merkelize indexed oracles");
}

template<typename FieldT, typename MT_hash_type>
std::vector<FieldT> bcs_indexer<FieldT, MT_hash_type>::obtain_verifier_random_message(
    const verifier_random_message_handle &random_message)
{
    throw std::invalid_argument("Should not be calling this on the indexing IOP");
}

template<typename FieldT, typename MT_hash_type>
FieldT bcs_indexer<FieldT, MT_hash_type>::obtain_query_response(const query_handle &query)
{
    throw std::invalid_argument("Should not be calling this on the indexing IOP");
}

template<typename FieldT, typename MT_hash_type>
bcs_verifier_index<FieldT, MT_hash_type> bcs_indexer<FieldT, MT_hash_type>::get_verifier_index()
{
    bcs_verifier_index<FieldT, MT_hash_type> index;
    for (size_t i = 0; i < this->MTs_processed_; i++)
    {
        index.index_MT_roots_.emplace_back(this->Merkle_trees_[i].get_root());
    }
    index.indexed_messages_ = this->prover_messages_;
    return index;
}

template<typename FieldT, typename MT_hash_type>
bcs_prover_index<FieldT, MT_hash_type> bcs_indexer<FieldT, MT_hash_type>::get_bcs_prover_index()
{
    if (this->get_prover_index_has_been_called_)
    {
        printf("Prover index has already been extracted from object "
            "due to memory optimizations, this operation can only be done once.\n");
        throw std::invalid_argument("Prover index has already been extracted");
    }
    bcs_prover_index<FieldT, MT_hash_type> index;
    index.index_MTs_ = this->Merkle_trees_;
    index.index_MTs_.erase(index.index_MTs_.begin() + this->MTs_processed_,
        index.index_MTs_.end());
    index.indexed_messages_ = this->prover_messages_;
    index.indexed_messages_.erase(
        index.indexed_messages_.begin() + this->num_prover_messages_at_end_of_round_[0],
        index.indexed_messages_.end());
    std::swap(index.iop_index_.all_oracle_evals_, this->indexed_oracles_);
    index.iop_index_.prover_messages_ = index.indexed_messages_;
    this->get_prover_index_has_been_called_ = true;
    return index;
}

} // namespace libiop
