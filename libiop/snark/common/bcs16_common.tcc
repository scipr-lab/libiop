#include <algorithm>
#include <numeric>
#include <set>

#include "libiop/algebra/fft.hpp"
#include "libiop/common/profiling.hpp"

namespace libiop {

template<typename FieldT>
std::size_t bcs16_transformation_transcript<FieldT>::IOP_size_in_bytes() const
{
    const std::size_t prover_messages_length =
        std::accumulate(this->prover_messages_.begin(),
                        this->prover_messages_.end(),
                        0,
                        [] (const std::size_t av,
                            const std::vector<FieldT> &msg) {
                            return av + msg.size();
                        });
    const std::size_t prover_messages_size =
        sizeof(FieldT) * prover_messages_length;

    const std::size_t query_responses_length =
        std::accumulate(this->query_responses_.begin(),
                        this->query_responses_.end(),
                        0,
                        [] (const std::size_t av,
                            const std::vector<std::vector<FieldT> > &resp) {
                            return av + (resp.empty() ? 0 : resp.size() * resp[0].size());
                        });
    const std::size_t query_responses_size =
        sizeof(FieldT) * query_responses_length;

    return (prover_messages_size +
            query_responses_size);
}


template<typename FieldT>
std::size_t bcs16_transformation_transcript<FieldT>::BCS_size_in_bytes() const
{
    const std::size_t MT_roots_size =
        std::accumulate(this->MT_roots_.begin(),
                        this->MT_roots_.end(),
                        0,
                        [] (const std::size_t av, const hash_digest &h) {
                            return av + h.size();
                        });

    const std::size_t MT_multi_membership_proofs_size =
        std::accumulate(this->MT_multi_membership_proofs_.begin(),
                        this->MT_multi_membership_proofs_.end(),
                        0,
                        [] (const std::size_t av,
                            const merkle_tree_multi_membership_proof &pi) {
                            return av + pi.size_in_bytes();
                        });

    return (MT_roots_size +
            MT_multi_membership_proofs_size);
}

template<typename FieldT>
std::size_t bcs16_transformation_transcript<FieldT>::size_in_bytes() const
{
    return (this->IOP_size_in_bytes() +
            this->BCS_size_in_bytes());
}

template<typename FieldT>
std::size_t bcs16_transformation_transcript<FieldT>::BCS_size_in_bytes_without_pruning() const
{
    const std::size_t MT_roots_size =
        std::accumulate(this->MT_roots_.begin(),
                        this->MT_roots_.end(),
                        0,
                        [] (const std::size_t av, const hash_digest &h) {
                            return av + h.size();
                        });

    const std::size_t digest_size_bytes = this->MT_roots_[0].size();

    return (MT_roots_size + digest_size_bytes * total_depth_without_pruning);
}

template<typename FieldT>
std::size_t bcs16_transformation_transcript<FieldT>::size_in_bytes_without_pruning() const
{
    return (this->IOP_size_in_bytes() +
            this->BCS_size_in_bytes_without_pruning());
}

template<typename FieldT>
bcs16_protocol<FieldT>::bcs16_protocol(const bcs16_transformation_parameters<FieldT> &parameters) :
    iop_protocol<FieldT>(),
    parameters_(parameters)
{
    this->digest_len_bytes_ = 2 * (this->parameters_.security_parameter / 8);
    printf("\nBCS parameters\n");
    print_indent(); printf("* digest_len (bytes) = %zu\n", this->digest_len_bytes_);
    print_indent(); printf("* digest_len (bits) = %zu\n", 8 * this->digest_len_bytes_);
}

template<typename FieldT>
void bcs16_protocol<FieldT>::seal_interaction_registrations()
{
    iop_protocol<FieldT>::seal_interaction_registrations();

    /* Now that all the interactions have been registered, we know how
       many messages there are, so we can prepare the Merkle trees. */
    for (std::size_t round = 0; round < this->num_interaction_rounds_; ++round)
    {
        const domain_to_oracles_map mapping = this->oracles_in_round(round);
        const round_parameters<FieldT> round_params = this->get_round_parameters(round);

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
            const merkle_tree<FieldT> MT(size,
                                         this->parameters_.field_hasher,
                                         this->parameters_.zk_hasher,
                                         this->parameters_.compression_hasher,
                                         this->digest_len_bytes_,
                                         make_zk,
                                         this->parameters_.security_parameter);
            this->Merkle_trees_.emplace_back(MT);
        }
    }
}

template<typename FieldT>
void bcs16_protocol<FieldT>::set_round_parameters(const round_parameters<FieldT> &params) {
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

template<typename FieldT>
round_parameters<FieldT> bcs16_protocol<FieldT>::get_round_parameters(const std::size_t round) {
    if (round >= this->round_params_.size()) {
        return round_parameters<FieldT>();
    }
    return this->round_params_[round];
}

template<typename FieldT>
std::size_t bcs16_protocol<FieldT>::obtain_random_query_position(const random_query_position_handle &position)
{
    /* Always use the last pseudorandom state */
    const std::size_t subspace_size = this->domains_[
        this->random_query_position_registrations_[position.id()].domain().id()].num_elements();
    const std::size_t result =
        this->parameters_.integer_randomness_extractor(*this->pseudorandom_state_.rbegin(),
                                                       position.id(),
                                                       subspace_size);
    return result;
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

} // namespace libiop
