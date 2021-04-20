#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>
#include "libiop/algebra/field_subset/subspace.hpp"
#include "libiop/bcs/bcs_common.hpp"
#include "libiop/bcs/hashing/blake2b.hpp"

namespace libiop {

template<typename FieldT, typename MT_root_hash>
void ligero_snark_parameters<FieldT, MT_root_hash>::describe()
{
    libff::print_indent(); printf("Interleaved R1CS SNARK parameters:\n");
    libff::print_indent(); printf("Security level: %zu\n", security_level_);
    libff::print_indent(); printf("Height/width ratio: %f\n", height_width_ratio_);
    libff::print_indent(); printf("RS extra dimensions: %zu\n", RS_extra_dimensions_);
    libff::print_indent(); printf("Zero-knowledge: %d\n", make_zk_);
    libff::print_indent(); printf("Domain type = %s\n", field_subset_type_names[this->domain_type_]);
}

template<typename FieldT, typename MT_root_hash>
ligero_iop_parameters<FieldT> obtain_iop_parameters_from_ligero_snark_params(
    const ligero_snark_parameters<FieldT, MT_root_hash> &parameters,
    const std::size_t num_constraints,
    const std::size_t num_variables)
{
    ligero_iop_parameters<FieldT> iop_parameters(
        parameters.security_level_,
        parameters.LDT_reducer_soundness_type_,
        parameters.RS_extra_dimensions_,
        parameters.height_width_ratio_,
        parameters.make_zk_,
        parameters.domain_type_,
        num_constraints,
        num_variables);
    iop_parameters.print();
    return iop_parameters;
}

template<typename FieldT, typename MT_root_hash>
ligero_snark_argument<FieldT, MT_root_hash> ligero_snark_prover(
    const r1cs_constraint_system<FieldT> &constraint_system,
    const r1cs_primary_input<FieldT> &primary_input,
    const r1cs_auxiliary_input<FieldT> &auxiliary_input,
    const ligero_snark_parameters<FieldT, MT_root_hash> &parameters)
{
    libff::enter_block("Ligero SNARK prover");
    const ligero_iop_parameters<FieldT> iop_params =
        obtain_iop_parameters_from_ligero_snark_params<FieldT>(
            parameters,
            constraint_system.num_constraints(),
            constraint_system.num_variables());

    bcs_prover<FieldT, MT_root_hash> IOP(parameters.bcs_params_);
    ligero_iop<FieldT> full_protocol(IOP,
                                     constraint_system,
                                     iop_params);
    full_protocol.register_interactions();
    IOP.seal_interaction_registrations();
    full_protocol.register_queries();
    IOP.seal_query_registrations();

    full_protocol.produce_proof(primary_input, auxiliary_input);

    libff::enter_block("Obtain transcript");
    const ligero_snark_argument<FieldT, MT_root_hash> transcript = IOP.get_transcript();
    libff::leave_block("Obtain transcript");

    IOP.describe_sizes();

    libff::leave_block("Ligero SNARK prover");
    return transcript;
}

template<typename FieldT, typename MT_root_hash>
bool ligero_snark_verifier(const r1cs_constraint_system<FieldT> &constraint_system,
                           const r1cs_primary_input<FieldT> &primary_input,
                           const ligero_snark_argument<FieldT, MT_root_hash> &proof,
                           const ligero_snark_parameters<FieldT, MT_root_hash> &parameters)
{
    libff::enter_block("Ligero SNARK verifier");
    const ligero_iop_parameters<FieldT> iop_params =
        obtain_iop_parameters_from_ligero_snark_params<FieldT>(
            parameters,
            constraint_system.num_constraints(),
            constraint_system.num_variables());

    bcs_verifier<FieldT, MT_root_hash> IOP(parameters.bcs_params_, proof);

    ligero_iop<FieldT> full_protocol(IOP,
                                     constraint_system,
                                     iop_params);
    full_protocol.register_interactions();
    IOP.seal_interaction_registrations();
    full_protocol.register_queries();
    IOP.seal_query_registrations();

    libff::enter_block("Check semantic validity of IOP transcript");
    const bool IOP_transcript_valid = IOP.transcript_is_valid();
    libff::leave_block("Check semantic validity of IOP transcript");

    libff::enter_block("Check verifier predicate");
    const bool full_protocol_accepts = full_protocol.verifier_predicate(primary_input);
    libff::leave_block("Check verifier predicate");

    libff::print_indent(); printf("* IOP transcript valid: %s\n", IOP_transcript_valid ? "true" : "false");
    libff::print_indent(); printf("* Full protocol decision predicate satisfied: %s\n", full_protocol_accepts ? "true" : "false");
    const bool decision = IOP_transcript_valid && full_protocol_accepts;
    libff::leave_block("Ligero SNARK verifier");

    return decision;
}

} // namespace libiop
