/**@file
*****************************************************************************
FRI low-degree test on a dummy oracle, compiled with the BCS transform (used
for profiling FRI).
*****************************************************************************
* @author     This file is part of libiop (see AUTHORS)
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/
#ifndef LIBIOP_SNARK_FRI_SNARK_HPP_
#define LIBIOP_SNARK_FRI_SNARK_HPP_

#include <cstddef>
#include <iostream>

#include "libiop/bcs/bcs_common.hpp"
#include "libiop/bcs/bcs_prover.hpp"
#include "libiop/bcs/bcs_verifier.hpp"

namespace libiop {

template<typename FieldT>
struct FRI_snark_parameters {
    std::size_t codeword_domain_dim_;
    std::size_t security_level_;
    libiop::bcs_hash_type hash_enum_;
    std::size_t RS_extra_dimensions_;
    std::size_t localization_parameter_;
    std::vector<std::size_t> localization_parameter_array_;
    std::size_t num_interactive_repetitions_;
    std::size_t num_query_repetitions_;

    std::size_t num_oracles_;

    libff::field_type field_type_;

    void describe();
};

template<typename FieldT, typename hash_type>
using FRI_snark_proof = bcs_transformation_transcript<FieldT, hash_type>;

template<typename FieldT, typename hash_type>
void FRI_snark_print_detailed_argument_size(
    FRI_snark_parameters<FieldT> params,
    FRI_snark_proof<FieldT, hash_type> argument);

template<typename FieldT, typename hash_type>
FRI_snark_proof<FieldT, hash_type> FRI_snark_prover(const FRI_snark_parameters<FieldT> &parameters);

template<typename FieldT, typename hash_type>
bool FRI_snark_verifier(const FRI_snark_proof<FieldT, hash_type> &proof,
                        const FRI_snark_parameters<FieldT> &parameters);

} // namespace libiop

#include "libiop/snark/fri_snark.tcc"

#endif // LIBIOP_SNARK_FRI_SNARK_HPP_
