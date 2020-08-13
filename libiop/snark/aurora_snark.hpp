/**@file
 *****************************************************************************
 R1CS zk-STARK obtained by combining our IOP for R1CS and the BCS16
 transformation.
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_SNARK_AURORA_SNARK_HPP_
#define LIBIOP_SNARK_AURORA_SNARK_HPP_

#include <cstddef>
#include <iostream>

#include "libiop/protocols/aurora_iop.hpp"
#include "libiop/protocols/ldt/fri/fri_ldt.hpp"
#include "libiop/protocols/ldt/ldt_reducer.hpp"
#include "libiop/bcs/bcs_common.hpp"
#include "libiop/bcs/bcs_prover.hpp"
#include "libiop/bcs/bcs_verifier.hpp"
#include "libiop/relations/r1cs.hpp"

namespace libiop {

template<typename FieldT, typename hash_type>
class aurora_snark_parameters {
    protected:
    size_t security_parameter_;
    LDT_reducer_soundness_type LDT_reducer_soundness_type_;
    FRI_soundness_type FRI_soundness_type_;
    size_t RS_extra_dimensions_;
    bool make_zk_;
    field_subset_type domain_type_;
    size_t num_constraints_;
    size_t num_variables_;

    std::vector<size_t> FRI_localization_parameter_array_;
    size_t FRI_localization_parameter_;

    void initialize_bcs_params(const bcs_hash_type hash_enum);
    void initialize_iop_params();
    public:
    aurora_snark_parameters(const size_t security_parameter,
                            const LDT_reducer_soundness_type ldt_reducer_soundness_type,
                            const FRI_soundness_type fri_soundness_type,
                            const bcs_hash_type hash_enum,
                            const std::vector<size_t> FRI_localization_parameter_array,
                            const size_t RS_extra_dimensions,
                            const bool make_zk,
                            const field_subset_type domain_type,
                            const size_t num_constraints,
                            const size_t num_variables);

    aurora_snark_parameters(const size_t security_parameter,
                            const LDT_reducer_soundness_type ldt_reducer_soundness_type,
                            const FRI_soundness_type fri_soundness_type,
                            const bcs_hash_type hash_enum,
                            const size_t FRI_localization_parameter,
                            const size_t RS_extra_dimensions,
                            const bool make_zk,
                            const field_subset_type domain_type,
                            const size_t num_constraints,
                            const size_t num_variables);

    void reset_fri_localization_parameters(const std::vector<size_t> FRI_localization_parameter_array);
    void print() const;

    bcs_transformation_parameters<FieldT, hash_type> bcs_params_;
    aurora_iop_parameters<FieldT> iop_params_;
};

template<typename FieldT, typename hash_type>
using aurora_snark_argument = bcs_transformation_transcript<FieldT, hash_type>;

template<typename FieldT, typename hash_type>
aurora_snark_argument<FieldT, hash_type> aurora_snark_prover(
    const r1cs_constraint_system<FieldT> &constraint_system,
    const r1cs_primary_input<FieldT> &primary_input,
    const r1cs_auxiliary_input<FieldT> &auxiliary_input,
    const aurora_snark_parameters<FieldT, hash_type> &parameters);

template<typename FieldT, typename hash_type>
bool aurora_snark_verifier(const r1cs_constraint_system<FieldT> &constraint_system,
                           const r1cs_primary_input<FieldT> &primary_input,
                           const aurora_snark_argument<FieldT, hash_type> &proof,
                           const aurora_snark_parameters<FieldT, hash_type> &parameters);


} // namespace libiop

#include "libiop/snark/aurora_snark.tcc"

#endif // LIBIOP_SNARK_AURORA_SNARK_HPP_
