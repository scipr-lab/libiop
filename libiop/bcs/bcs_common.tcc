#include <algorithm>
#include <numeric>
#include <set>
#include <libff/algebra/field_utils/bigint.hpp>
#include <libff/common/profiling.hpp>

#include "libiop/algebra/fft.hpp"

namespace libiop {

template<typename FieldT, typename MT_root_hash>
std::size_t bcs_transformation_transcript<FieldT, MT_root_hash>::IOP_size_in_bytes() const
{
    const size_t field_size =
        (libff::log_of_field_size_helper<FieldT>(FieldT::zero()) + 7) / 8;
    const std::size_t prover_messages_length =
        std::accumulate(this->prover_messages_.begin(),
                        this->prover_messages_.end(),
                        0,
                        [] (const std::size_t av,
                            const std::vector<FieldT> &msg) {
                            return av + msg.size();
                        });
    const std::size_t prover_messages_size =
        field_size * prover_messages_length;

    const std::size_t query_responses_length =
        std::accumulate(this->query_responses_.begin(),
                        this->query_responses_.end(),
                        0,
                        [] (const std::size_t av,
                            const std::vector<std::vector<FieldT> > &resp) {
                            return av + (resp.empty() ? 0 : resp.size() * resp[0].size());
                        });
    const std::size_t query_responses_size =
        field_size * query_responses_length;

    return (prover_messages_size +
            query_responses_size);
}

template<typename FieldT, typename MT_root_hash>
std::size_t bcs_transformation_transcript<FieldT, MT_root_hash>::BCS_size_in_bytes() const
{
    const std::size_t MT_roots_size =
        std::accumulate(this->MT_roots_.begin(),
                        this->MT_roots_.end(),
                        0,
                        [] (const std::size_t av, const MT_root_hash &h) {
                            return av + get_hash_size<MT_root_hash>(h);
                        });

    const std::size_t MT_set_membership_proofs_size =
        std::accumulate(this->MT_set_membership_proofs_.begin(),
                        this->MT_set_membership_proofs_.end(),
                        0,
                        [] (const std::size_t av,
                            const merkle_tree_set_membership_proof<MT_root_hash> &pi) {
                            return av + pi.size_in_bytes();
                        });
    
    const std::size_t pow_size = get_hash_size<MT_root_hash>(this->proof_of_work_);

    return (MT_roots_size +
            MT_set_membership_proofs_size + 
            pow_size);
}

template<typename FieldT, typename MT_root_hash>
std::size_t bcs_transformation_transcript<FieldT, MT_root_hash>::size_in_bytes() const
{
    return (this->IOP_size_in_bytes() +
            this->BCS_size_in_bytes());
}

template<typename FieldT, typename MT_root_hash>
std::size_t bcs_transformation_transcript<FieldT, MT_root_hash>::BCS_size_in_bytes_without_pruning() const
{
    const std::size_t MT_roots_size =
        std::accumulate(this->MT_roots_.begin(),
                        this->MT_roots_.end(),
                        0,
                        [] (const std::size_t av, const MT_root_hash &h) {
                            return av + get_hash_size<MT_root_hash>(h);
                        });

    const std::size_t digest_size_bytes = get_hash_size<MT_root_hash>(this->MT_roots_[0]);
    const std::size_t pow_size = get_hash_size<MT_root_hash>(this->proof_of_work_);

    return (MT_roots_size + pow_size + digest_size_bytes * total_depth_without_pruning);
}

/** TODO: Come back and cleanup this serialization code. 
 * Its really messy, and quite unnecessarily so.
 */
template<typename FieldT>
std::ostream& serialize_FieldT(
    std::ostream &out, const FieldT &v)
{
    bigint<FieldT::num_limbs> b = v.as_bigint();
    for (size_t j = 0; j < FieldT::num_limbs; j++)
    {
        out << b.data[j];
        out << ",";
    }
    return out;
}

template<typename FieldT>
FieldT deserialize_FieldT(
    std::istream &in)
{
    bigint<FieldT::num_limbs> b;
    char delimiter;
    for (size_t j = 0; j < FieldT::num_limbs; j++)
    {
        in >> b.data[j];
        in >> delimiter;
        assert(delimiter == char(','));
    }
    return FieldT(b);
}

template<typename FieldT>
std::ostream& serialize_Field_Elem_vec(
    std::ostream &out, const std::vector<FieldT> &v)
{
    out << v.size();
    out << ",";
    for (size_t i = 0; i < v.size(); i++)
    {
        serialize_FieldT(out, v[i]);
    }
    return out;
}

template<typename FieldT>
std::istream& deserialize_Field_Elem_vec(
    std::istream &in, std::vector<FieldT> &v)
{
    size_t size;
    in >> size;
    char delimiter;
    in >> delimiter;
    assert(delimiter == char(','));
    for (size_t i = 0; i < size; i++)
    {
        FieldT e = deserialize_FieldT<FieldT>(in);
        v.emplace_back(e);
    }
    return in;
}

template<typename FieldT>
std::ostream& serialize_Field_Elem_vec_of_vec(
    std::ostream &out, const std::vector<std::vector<FieldT>> &v)
{
    out << v.size();
    out << ",";
    for (size_t i = 0; i < v.size(); i++)
    {
        serialize_Field_Elem_vec<FieldT>(out, v[i]);
    }
    return out;
}

template<typename FieldT>
std::istream& deserialize_Field_Elem_vec_of_vec(
    std::istream &in, std::vector<std::vector<FieldT>> &v)
{
    size_t size;
    in >> size;
    char delimiter;
    in >> delimiter;
    assert(delimiter == char(','));
    for (size_t i = 0; i < size; i++)
    {
        std::vector<FieldT> vec;
        deserialize_Field_Elem_vec(in, vec);
        v.emplace_back(vec);
    }
    return in;
}

/* Currently only works for non-zk case. */
template<typename FieldT>
std::ostream& serialize_vec_of_MT_proofs(
    std::ostream &out, const std::vector<merkle_tree_set_membership_proof<FieldT>> &v)
{
    out << v.size();
    out << ",";
    for (size_t i = 0; i < v.size(); i++)
    {
        out << v[i].auxiliary_hashes.size();
        out << ",";
        for (size_t j = 0; j < v[i].auxiliary_hashes.size(); j++)
        {
            serialize_FieldT<FieldT>(out, v[i].auxiliary_hashes[j]);
        }
    }
    return out;
}

template<typename FieldT>
std::istream& deserialize_vec_of_MT_proofs(
    std::istream &in, std::vector<merkle_tree_set_membership_proof<FieldT>> &v)
{
    size_t size;
    in >> size;
    char delimiter;
    in >> delimiter;
    assert(delimiter == char(','));
    for (size_t i = 0; i < size; i++)
    {
        size_t aux_hash_size;
        in >> aux_hash_size;
        in >> delimiter;
        assert(delimiter == char(','));
        std::vector<FieldT> auxiliary_hashes;
        for (size_t j = 0; j < aux_hash_size; j++)
        {
            auxiliary_hashes.emplace_back(deserialize_FieldT<FieldT>(in));
        }
        merkle_tree_set_membership_proof<FieldT> prf;
        prf.auxiliary_hashes = auxiliary_hashes;
        v.emplace_back(prf);
    }
    return in;
}

template<typename FieldT>
std::istream& deserialize_Field_Elem_vec_of_vec_of_vec(
    std::istream &in, std::vector<std::vector<std::vector<FieldT>>> &v)
{
    size_t size;
    in >> size;
    char delimiter;
    in >> delimiter;
    assert(delimiter == char(','));
    for (size_t i = 0; i < size; i++)
    {
        std::vector<std::vector<FieldT>> vec;
        deserialize_Field_Elem_vec_of_vec(in, vec);
        v.emplace_back(vec);
    }
    return in;
}

template<typename FieldT>
std::ostream& serialize_Field_Elem_vec_of_vec_of_vec(
    std::ostream &out, const std::vector<std::vector<std::vector<FieldT>>> &v)
{
    out << v.size();
    out << ",";
    for (size_t i = 0; i < v.size(); i++)
    {
        serialize_Field_Elem_vec_of_vec<FieldT>(out, v[i]);
    }
    return out;
}

std::ostream& serialize_size_t_vec_of_vec(
    std::ostream &out, const std::vector<std::vector<size_t>> &v)
{
    out << v.size();
    out << ",";
    for (size_t i = 0; i < v.size(); i++)
    {
        out << v[i].size();
        out << ",";
        for (size_t j = 0; j < v[i].size(); j++)
        {
            out << v[i][j];
            out << ",";
        }
    }
    return out;
}

// TODO: Left off here
std::istream& deserialize_size_t_vec_of_vec(
    std::istream &in, std::vector<std::vector<size_t>> &v)
{
    size_t size;
    in >> size;
    char delimiter;
    in >> delimiter;
    assert(delimiter == char(','));
    for (size_t i = 0; i < size; i++)
    {
        std::vector<size_t> vec;
        size_t vec_len;
        in >> vec_len;
        in >> delimiter;
        assert(delimiter == char(','));
        for (size_t j = 0; j < vec_len; j++)
        {
            size_t cur;
            in >> cur;
            vec.emplace_back(cur);
            in >> delimiter;
            assert(delimiter == char(','));
        }
        v.emplace_back(vec);
    }
    return in;
}

template<typename FieldT, typename MT_hash_type>
std::ostream& serialize_transcript_internal(
    typename libff::enable_if<libff::is_multiplicative<FieldT>::value, FieldT>::type,
    typename libff::enable_if<std::is_same<MT_hash_type, FieldT>::value, FieldT>::type,
    std::ostream &out, const bcs_transformation_transcript<FieldT, MT_hash_type> &t)
{
    // algebraic hash on multiplicative field
    serialize_Field_Elem_vec_of_vec<FieldT>(out, t.prover_messages_);
    out << "\n";
    serialize_Field_Elem_vec<FieldT>(out, t.MT_roots_);
    out << "\n";
    serialize_size_t_vec_of_vec(out, t.query_positions_);
    out << "\n";
    serialize_Field_Elem_vec_of_vec_of_vec(out, t.query_responses_);
    out << "\n";
    serialize_size_t_vec_of_vec(out, t.MT_leaf_positions_);
    out << "\n";
    serialize_vec_of_MT_proofs(out, t.MT_set_membership_proofs_);
    out << "\n";
    return out;
}

template<typename FieldT, typename MT_hash_type>
std::istream& deserialize_transcript_internal(
    typename libff::enable_if<libff::is_multiplicative<FieldT>::value, FieldT>::type,
    typename libff::enable_if<std::is_same<MT_hash_type, FieldT>::value, FieldT>::type,
    std::istream &in, bcs_transformation_transcript<FieldT, MT_hash_type> &t)
{
    // algebraic hash on multiplicative field
    char delimiter;
    deserialize_Field_Elem_vec_of_vec<FieldT>(in, t.prover_messages_);
    in.get(delimiter);
    assert(delimiter == '\n');
    deserialize_Field_Elem_vec<FieldT>(in, t.MT_roots_);
    in.get(delimiter);
    assert(delimiter == '\n');
    deserialize_size_t_vec_of_vec(in, t.query_positions_);
    in.get(delimiter);
    assert(delimiter == '\n');
    deserialize_Field_Elem_vec_of_vec_of_vec(in, t.query_responses_);
    in.get(delimiter);
    assert(delimiter == '\n');
    deserialize_size_t_vec_of_vec(in, t.MT_leaf_positions_);
    in.get(delimiter);
    assert(delimiter == '\n');
    deserialize_vec_of_MT_proofs(in, t.MT_set_membership_proofs_);
    in.get(delimiter);
    assert(delimiter == '\n');
    return in;
}

template<typename FieldT, typename MT_hash_type>
std::ostream& serialize_transcript_internal(
    typename libff::enable_if<!libff::is_multiplicative<FieldT>::value, FieldT>::type,
    typename libff::enable_if<!std::is_same<MT_hash_type, FieldT>::value, FieldT>::type,
    std::ostream &out, const bcs_transformation_transcript<FieldT, MT_hash_type> &t)
{
    printf("Serializing on binary fields or non-algebraic hashes is not implemented\n");
    return out;
}

template<typename FieldT, typename MT_hash_type>
std::istream& deserialize_transcript_internal(
    typename libff::enable_if<!libff::is_multiplicative<FieldT>::value, FieldT>::type,
    typename libff::enable_if<!std::is_same<MT_hash_type, FieldT>::value, FieldT>::type,
    std::istream &in, bcs_transformation_transcript<FieldT, MT_hash_type> &t)
{
    printf("Deserializing on binary fields or non-algebraic hashes is not implemented\n");
    return in;
}

template<typename FieldT, typename MT_hash_type>
std::ostream& bcs_transformation_transcript<FieldT, MT_hash_type>::serialize(std::ostream &out) const
{
    return serialize_transcript_internal<FieldT, MT_hash_type>(FieldT::zero(), FieldT::zero(), out, *this);
}

template<typename FieldT, typename MT_hash_type>
std::istream& bcs_transformation_transcript<FieldT, MT_hash_type>::deserialize(std::istream &in)
{
    return deserialize_transcript_internal<FieldT, MT_hash_type>(FieldT::zero(), FieldT::zero(), in, *this);
}

template<typename FieldT, typename MT_root_hash>
std::size_t bcs_transformation_transcript<FieldT, MT_root_hash>::size_in_bytes_without_pruning() const
{
    return (this->IOP_size_in_bytes() +
            this->BCS_size_in_bytes_without_pruning());
}

template<typename FieldT, typename MT_root_hash>
bcs_protocol<FieldT, MT_root_hash>::bcs_protocol(
    const bcs_transformation_parameters<FieldT, MT_root_hash> &parameters) :
    iop_protocol<FieldT>(),
    parameters_(parameters)
{
    this->digest_len_bytes_ = 2 * (this->parameters_.security_parameter / 8);
    this->hashchain_ = this->parameters_.hashchain_->new_hashchain();
    printf("\nBCS parameters\n");
    libff::print_indent(); printf("* digest_len (bytes) = %zu\n", this->digest_len_bytes_);
    libff::print_indent(); printf("* digest_len (bits) = %zu\n", 8 * this->digest_len_bytes_);
    libff::print_indent(); printf("* hash_type = %s\n", bcs_hash_type_names[parameters.hash_enum]);
}

template<typename FieldT, typename MT_root_hash>
void bcs_protocol<FieldT, MT_root_hash>::register_proof_of_work()
{
    // The verifier squeezes out a challenge, the prover creates a pow,
    // and then that gets absorbed into the hash chain.
    // this->pow_challenge_handle_ = this->register_verifier_random_message(1);
    // this->pow_proof_handle_ = this->register_prover_message(1);
    this->pow_ = pow<FieldT, MT_root_hash>(this->parameters_.pow_params_, this->digest_len_bytes_);
}

template<typename FieldT, typename MT_root_hash>
void bcs_protocol<FieldT, MT_root_hash>::seal_interaction_registrations()
{
    // Do a proof of work as long as we are not in the indexer
    if (!(this->is_holographic_ && this->num_interaction_rounds_ == 1))
    {
        // To do the proof of work, we add a message-only round at the end of the interaction.
        // This means if the prover wants to try new queries, they must redo the proof of work.
        this->register_proof_of_work();
    }

    iop_protocol<FieldT>::seal_interaction_registrations();

    /* Now that all the interactions have been registered, we know how
       many messages there are, so we can prepare the Merkle trees. */
    for (std::size_t round = 0; round < this->num_interaction_rounds_; ++round)
    {
        const domain_to_oracles_map mapping = this->oracles_in_round_by_domain(round);
        const round_parameters<FieldT> round_params = this->get_round_parameters(round);

        /* Don't double instantiate holographic MTs in the prover */
        if (this->is_holographic_ && round == 0)
        {
            /** The prover will already have the indexed MTs instantiated here,
             *  the verifier won't. So if we are in the prover,
             *  the following check will pass*/
            if (this->Merkle_trees_.size() > 0)
            {
                continue;
            }
        }
        for (auto &kv : mapping)
        {
            // Make the Merkle tree ZK if at least one oracle in this merkle tree is zk.
            bool make_zk = false;
            for (size_t i = 0; i < kv.second.size(); i++)
            {
                const size_t id = kv.second[i].id();
                const oracle_registration registration = this->oracle_registrations_[id];
                if (registration.make_zk())
                {
                    make_zk = true;
                    break;
                }
            }
            /* Make this a per-oracle setting as well */
            const std::size_t size = this->domains_[kv.first.id()].num_elements() / round_params.quotient_map_size_;
            const merkle_tree<FieldT, MT_root_hash> MT(
                size,
                this->parameters_.leafhasher_,
                this->parameters_.compression_hasher,
                this->digest_len_bytes_,
                make_zk,
                this->parameters_.security_parameter);
            this->Merkle_trees_.emplace_back(MT);
        }
    }
}

template<typename FieldT, typename MT_root_hash>
void bcs_protocol<FieldT, MT_root_hash>::set_round_parameters(const round_parameters<FieldT> &params) {
    /* Rounds are 0 indexed. */
    std::size_t cur_round = this->num_interaction_rounds();
    if (cur_round == this->round_params_.size() - 1) {
        throw std::logic_error("Already set round parameters for this round");
    }
    /* Set blank round params until one round before current param. */
    while (this->round_params_.size() < cur_round) {
        this->round_params_.emplace_back(round_parameters<FieldT>());
    }
    this->round_params_.emplace_back(params);
}

template<typename FieldT, typename MT_root_hash>
std::vector<size_t> bcs_protocol<FieldT, MT_root_hash>::get_MT_depths() const
{
    std::vector<size_t> depths;
    for (size_t i = 0; i < this->Merkle_trees_.size(); i++)
    {
        depths.emplace_back(this->Merkle_trees_[i].depth());
    }
    return depths;
}

template<typename FieldT, typename MT_root_hash>
std::vector<bool> bcs_protocol<FieldT, MT_root_hash>::get_MT_zk_flags() const
{
    std::vector<bool> make_zk_flags;
    for (size_t i = 0; i < this->Merkle_trees_.size(); i++)
    {
        make_zk_flags.emplace_back(this->Merkle_trees_[i].zk());
    }
    return make_zk_flags;
}

template<typename FieldT, typename MT_root_hash>
std::vector<round_parameters<FieldT>> bcs_protocol<FieldT, MT_root_hash>::get_all_round_params() const
{
    std::vector<round_parameters<FieldT>> all_round_params;
    for (size_t i = 0; i < this->Merkle_trees_.size(); i++)
    {
        all_round_params.emplace_back(this->get_round_parameters(i));
    }
    return all_round_params;
}

template<typename FieldT, typename MT_root_hash>
round_parameters<FieldT> bcs_protocol<FieldT, MT_root_hash>::get_round_parameters(const std::size_t round) const {
    if (round >= this->round_params_.size()) {
        return round_parameters<FieldT>();
    }
    return this->round_params_[round];
}

template<typename FieldT, typename MT_root_hash>
std::size_t bcs_protocol<FieldT, MT_root_hash>::obtain_random_query_position(const random_query_position_handle &position)
{
    /* Obtains the next query position using the latest hashchain state. */
    /* TODO: Simplify this. */
    /* TODO: Make IOP infra obtain all query positions at once */
    const std::size_t subspace_size = this->domains_[
        this->random_query_position_registrations_[position.id()].domain().id()].num_elements();
    const std::size_t result =
        this->hashchain_->squeeze_query_positions(1, subspace_size)[0];
    return result;
}

template<typename FieldT, typename MT_root_hash>
void bcs_protocol<FieldT, MT_root_hash>::run_hashchain_for_round(
    const std::size_t round,
    const std::vector<MT_root_hash> round_MT_roots,
    const std::vector<std::vector<FieldT> > prover_messages)
{
    /* Assume the Merkle tree is already created. */
    for (auto MT_root : round_MT_roots)
    {
        this->hashchain_->absorb(MT_root);
    }

    /* Add the prover message hash as a "root" and update the pseudorandom state */
    this->absorb_prover_messages(round, prover_messages);
    this->squeeze_verifier_random_messages(round);
}

template<typename FieldT, typename MT_root_hash>
void bcs_protocol<FieldT, MT_root_hash>::absorb_prover_messages(
    const size_t round,
    const std::vector<std::vector<FieldT> > &all_prover_messages)
{
    /* Hash explicitly sent prover messages */
    const std::size_t min_message_id =
        (round == 0 ? 0 : this->num_prover_messages_at_end_of_round_[round - 1]);
    const std::size_t max_message_id = this->num_prover_messages_at_end_of_round_[round];

    std::vector<FieldT> message_concat = { FieldT::zero() };
    for (std::size_t message_id = min_message_id; message_id < max_message_id; ++message_id)
    {
        message_concat.insert(message_concat.end(),
                              all_prover_messages[message_id].begin(),
                              all_prover_messages[message_id].end());

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
    this->hashchain_->absorb(message_concat);
}

template<typename FieldT, typename MT_root_hash>
void bcs_protocol<FieldT, MT_root_hash>::squeeze_verifier_random_messages(
    const size_t ended_round)
{
    /* Squeeze verifier randomness */
    size_t start_verifier_random_message_index =
        this->num_verifier_random_messages_at_end_of_round_[ended_round];
    size_t end_verifier_random_message_index = (ended_round == this->num_interaction_rounds_ - 1) ? 0 :
        this->num_verifier_random_messages_at_end_of_round_[ended_round + 1];

    for (std::size_t i = start_verifier_random_message_index; i < end_verifier_random_message_index; i++)
    {
        const std::size_t message_length = this->verifier_random_message_registrations_[i].size();
        const std::vector<FieldT> result = this->hashchain_->squeeze(message_length);
        this->verifier_random_messages_.insert(std::pair<size_t, std::vector<FieldT>>(i, result));
    }
}

/* Serializes the provided oracle's evaluated contents into multiple rows,
   and appends these rows to all_oracles_evaluated_contents.
   It is serialized into multiple rows according to the quotient map relation in round params.
   The default round parameters result in oracle_evaluated_contents being appended as a single row. */
template<typename FieldT, typename MT_root_hash>
void bcs_protocol<FieldT, MT_root_hash>::serialize_leaf_data_by_round_params(
    const std::vector<FieldT> &oracle_evaluated_contents,
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
        all_evaluated_contents.emplace_back(std::move(vector_for_ith_element_of_each_coset));
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

template<typename FieldT>
std::size_t query_position_to_merkle_tree_position(const std::size_t query_position,
                                                   const std::size_t num_leaves,
                                                   const round_parameters<FieldT> &round_params)
{
    if (round_params.quotient_map_size_ == 1) {
        return query_position;
    }
    if (round_params.quotient_map_type_ == affine_subspace_type) {
        return query_position / round_params.quotient_map_size_;
    } else if (round_params.quotient_map_type_ == multiplicative_coset_type) {
        return query_position % num_leaves;
    }
    return 0;
}

template<typename FieldT, typename MT_root_hash>
void print_detailed_transcript_data(
    const bool holographic,
    const bcs_transformation_transcript<FieldT, MT_root_hash> &transcript,
    const bcs_transformation_parameters<FieldT, MT_root_hash> &params,
    const bcs_protocol<FieldT, MT_root_hash> bcs)
{
    /* Calculate round by round details */

    const std::vector<size_t> MT_depths = bcs.get_MT_depths();
    const std::vector<bool> make_zk = bcs.get_MT_zk_flags();
    const std::vector<round_parameters<FieldT>> round_params = bcs.get_all_round_params();

    const size_t digest_len_bytes = 2 * (params.security_parameter / 8);
    const size_t field_size = (libff::log_of_field_size_helper<FieldT>(FieldT::zero()) + 7) / 8;
    std::vector<size_t> two_to_one_hashes_by_round;
    std::vector<size_t> leaf_hashes_by_round;
    std::vector<size_t> zk_hashes_by_round;
    std::vector<size_t> IOP_size_by_round;
    std::vector<size_t> BCS_size_by_round;
    size_t total_prover_message_size;

    for (size_t round = 0; round < MT_depths.size(); round++)
    {
        const size_t MT_size = 1ull << MT_depths[round];
        merkle_tree<FieldT, MT_root_hash> MT(
            MT_size,
            params.leafhasher_,
            params.compression_hasher,
            digest_len_bytes,
            false,
            params.security_parameter);

        /** We have to merge the query positions that correspond to the same
            leaf after applying the round parameters */
        std::vector<size_t> query_positions;
        for (size_t i = 0; i < transcript.query_positions_[round].size(); i++)
        {
            size_t query_position = transcript.query_positions_[round][i];
            size_t MT_position = query_position_to_merkle_tree_position(
                query_position, MT_size, round_params[round]);
            if(std::find(query_positions.begin(), query_positions.end(), MT_position) == query_positions.end()) {
                /* query positions does not contain x */
                query_positions.emplace_back(MT_position);
            }
        }
        size_t num_two_to_one_hashes_in_round =
            MT.count_hashes_to_verify_set_membership_proof(
            query_positions);
        two_to_one_hashes_by_round.emplace_back(num_two_to_one_hashes_in_round);
        const size_t num_values_per_leaf = transcript.query_responses_[round][0].size();
        const size_t num_leaves = transcript.query_responses_[round].size();
        leaf_hashes_by_round.emplace_back(num_values_per_leaf * num_leaves);

        if (make_zk[round]) {
            zk_hashes_by_round.emplace_back(num_leaves);
        } else {
            zk_hashes_by_round.emplace_back(0);
        }

        // TODO: Should we change sizeof(FieldT) to FieldT::num_bits / 8
        IOP_size_by_round.emplace_back(
            leaf_hashes_by_round[round] * field_size);
        /* MT root + membership proof size (includes zk hash) */
        BCS_size_by_round.emplace_back(
            transcript.MT_set_membership_proofs_[round].size_in_bytes()
            + digest_len_bytes);
    }
    size_t num_prover_messages = 0;
    for (size_t i = 0; i < transcript.prover_messages_.size(); i++)
    {
        num_prover_messages += transcript.prover_messages_[i].size();
    }
    total_prover_message_size = num_prover_messages * field_size;
    if (holographic)
    {
        BCS_size_by_round[0] -= digest_len_bytes;
    }

    /* Print summary of argument size first */
    printf("\n");

    libff::print_indent(); printf("* Argument size in bytes (IOP): %zu\n", transcript.IOP_size_in_bytes());
    libff::print_indent(); printf("* Argument size in bytes (BCS): %zu\n", transcript.BCS_size_in_bytes());
    libff::print_indent(); printf("* Argument size in bytes (total): %zu\n", transcript.size_in_bytes());

    printf("\nIf we were to remove pruning of authentication paths in BCS,\n"
            "the argument would have the following sizes:\n");
    libff::print_indent(); printf("* Argument size in bytes (BCS, no pruning): %zu\n", transcript.BCS_size_in_bytes_without_pruning());
    libff::print_indent(); printf("* Argument size in bytes (total, no pruning): %zu\n", transcript.size_in_bytes_without_pruning());

    printf("\n");
    printf("total prover messages size: %lu\n", total_prover_message_size);
    const size_t total_two_to_one_hashes = std::accumulate(
        two_to_one_hashes_by_round.begin(), two_to_one_hashes_by_round.end(), 0);
    const size_t total_leaves_hashed = std::accumulate(
        leaf_hashes_by_round.begin(), leaf_hashes_by_round.end(), 0);
    const size_t total_zk_hashes = std::accumulate(
        zk_hashes_by_round.begin(), zk_hashes_by_round.end(), 0);
    const size_t total_hashes = total_two_to_one_hashes + total_leaves_hashed + total_zk_hashes;
    printf("total two to one hashes: %lu\n", total_two_to_one_hashes);
    printf("total leaves hashed: %lu\n", total_leaves_hashed);
    printf("total hashes: %lu\n", total_hashes);
    printf("\n");

    printf("Transcript info by round\n");
    printf("Per round IOP sizes don't include prover messages\n");
    for (size_t round = 0; round < MT_depths.size(); round++)
    {
        printf("\nround %lu\n", round);
        printf("MT_depth %lu\n", MT_depths[round]);
        printf("IOP size: %lu bytes\n", IOP_size_by_round[round]);
        printf("BCS size: %lu bytes\n", BCS_size_by_round[round]);
        printf("number of two to one hashes: %lu\n", two_to_one_hashes_by_round[round]);
        printf("number of leaves hashed: %lu\n", leaf_hashes_by_round[round]);
        if (make_zk[round])
        {
            printf("number of zk hashes: %lu\n", zk_hashes_by_round[round]);
        }

        std::vector<oracle_registration> oracles = bcs.get_oracle_registrations_by_round(round);
        printf("oracles in round: ");
        for (auto &oracle_itr : oracles)
        {
            printf("%s, ", oracle_itr.name().c_str());
        }
        printf("\n");
    }
    printf("\n\n");
}

} // namespace libiop
