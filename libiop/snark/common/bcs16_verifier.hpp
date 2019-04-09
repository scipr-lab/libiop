/**@file
 *****************************************************************************
 BCS16 transformation verifier, for converting an IOP verifier into a SNARK verifier
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_SNARK_COMMON_BCS16_VERIFIER_HPP_
#define LIBIOP_SNARK_COMMON_BCS16_VERIFIER_HPP_

#include <set>

#include "libiop/common/profiling.hpp"
#include "libiop/snark/common/bcs16_common.hpp"

namespace libiop {

template<typename FieldT>
class bcs16_verifier : public bcs16_protocol<FieldT> {
protected:
    bcs16_transformation_transcript<FieldT> transcript_;
    std::map<std::pair<std::size_t, std::size_t>, FieldT> oracle_id_and_pos_idx_to_value_;
    bool transcript_is_valid_;
public:
    bcs16_verifier(const bcs16_transformation_parameters<FieldT> &parameters,
                  const bcs16_transformation_transcript<FieldT> &transcript);

    /* This performs the actual verification: go through the rounds,
       computing the "random" verifier messages and checking the authentication
       paths. */
    virtual void seal_interaction_registrations();

    virtual void signal_prover_round_done(); /* throws an exception */

    virtual std::vector<FieldT> obtain_verifier_random_message(const verifier_random_message_handle &random_message);
    virtual FieldT get_oracle_evaluation_at_point(
        const oracle_handle_ptr &handle,
        const std::size_t evaluation_position,
        const bool record=false);
    virtual std::vector<FieldT> receive_prover_message(const prover_message_handle &message);

    bool transcript_is_valid() const;
protected:
    void parse_query_responses_from_transcript();
    std::vector<std::vector<FieldT>> query_responses_to_MT_leaf_responses(
        std::vector<std::size_t> &query_positions,
        std::vector<std::vector<FieldT> > &query_responses,
        const std::size_t round,
        const std::size_t num_leaves);
};

} // namespace libiop

#include "libiop/snark/common/bcs16_verifier.tcc"

#endif // LIBIOP_SNARK_COMMON_BCS16_VERIFIER_HPP_
