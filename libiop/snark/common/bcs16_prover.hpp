/**@file
 *****************************************************************************
 BCS16 prover, for converting an IOP prover into a SNARK prover
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_SNARK_COMMON_BCS16_PROVER_HPP_
#define LIBIOP_SNARK_COMMON_BCS16_PROVER_HPP_

#include <set>

#include "libiop/common/profiling.hpp"
#include "libiop/snark/common/bcs16_common.hpp"

namespace libiop {

template<typename FieldT>
class bcs16_prover : public bcs16_protocol<FieldT> {
protected:
    std::size_t MTs_processed_ = 0;
    void serialize_leaf_data_by_round_params(std::vector<FieldT> &oracle_evaluated_contents,
                                             std::vector<std::vector<FieldT>> &all_oracles_evaluated_contents,
                                             const domain_handle &evaluation_domain,
                                             const round_parameters<FieldT> &round_params);
public:
    bcs16_prover(const bcs16_transformation_parameters<FieldT> &parameters);
    /*
      The overloaded method for signal_prover_round_done performs
      hashing of all oracles and prover messages submitted in the
      current round.
     */
    virtual void signal_prover_round_done();

    /*
      We also overload verifier's randomness extraction functions
      below to use the extracted values.
    */
    virtual std::vector<FieldT> obtain_verifier_random_message(const verifier_random_message_handle &random_message);

    virtual FieldT obtain_query_response(const query_handle &query);
    bcs16_transformation_transcript<FieldT> get_transcript();

    std::size_t MT_size() const;
    std::size_t state_size() const;

    void describe_sizes() const;
};

} // namespace libiop

#include "libiop/snark/common/bcs16_prover.tcc"

#endif // LIBIOP_SNARK_COMMON_BCS16_PROVER_HPP_
