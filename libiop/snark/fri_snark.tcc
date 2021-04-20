#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>
#include "libiop/algebra/field_subset/subspace.hpp"
#include "libiop/protocols/fri_iop.hpp"
#include "libiop/bcs/bcs_common.hpp"
#include "libiop/bcs/common_bcs_parameters.hpp"

namespace libiop {

template<typename FieldT>
void FRI_snark_parameters<FieldT>::describe()
{
    libff::print_indent(); printf("* FRI SNARK parameters:\n");
    libff::print_indent(); printf("* Security level: %zu\n", security_level_);
    libff::print_indent(); printf("* RS extra dimensions: %zu\n", RS_extra_dimensions_);
    libff::print_indent(); printf("* Localization parameter: %zu\n", localization_parameter_);
    libff::print_indent(); printf("* Localization parameter array: %zu\n", localization_parameter_array_);
    libff::print_indent(); printf("* Num query repetitions: %zu\n", num_query_repetitions_);
}

template<typename FieldT, typename hash_type>
const std::pair<bcs_transformation_parameters<FieldT, hash_type>,
                FRI_iop_protocol_parameters>
obtain_bcs_and_FRI_parameters_from_FRI_snark_parameters(const FRI_snark_parameters<FieldT> &parameters)
{
    bcs_transformation_parameters<FieldT, hash_type> bcs_parameters =
        default_bcs_params<FieldT, hash_type>(
            parameters.hash_enum_, parameters.security_level_, parameters.codeword_domain_dim_);

    FRI_iop_protocol_parameters FRI_parameters;
    FRI_parameters.RS_extra_dimensions_ = parameters.RS_extra_dimensions_;
    FRI_parameters.codeword_domain_dim_ = parameters.codeword_domain_dim_;
    FRI_parameters.localization_parameter_ = parameters.localization_parameter_;
    FRI_parameters.localization_parameter_array_ = parameters.localization_parameter_array_;
    FRI_parameters.num_query_repetitions_ = parameters.num_query_repetitions_;
    FRI_parameters.num_interactive_repetitions_ = parameters.num_interactive_repetitions_;
    FRI_parameters.num_oracles_ = parameters.num_oracles_;
    FRI_parameters.field_type_ = parameters.field_type_;

    return std::make_pair(bcs_parameters, FRI_parameters);
}

template<typename FieldT, typename hash_type>
FRI_snark_proof<FieldT, hash_type> FRI_snark_prover(const FRI_snark_parameters<FieldT> &parameters)
{
    libff::enter_block("FRI SNARK prover");
    const std::pair<bcs_transformation_parameters<FieldT, hash_type>,
                    FRI_iop_protocol_parameters>
        bcs_and_FRI_parameters =
        obtain_bcs_and_FRI_parameters_from_FRI_snark_parameters<FieldT, hash_type>(parameters);

    bcs_transformation_parameters<FieldT, hash_type> bcs_parameters = bcs_and_FRI_parameters.first;
    bcs_prover<FieldT, hash_type> IOP(bcs_parameters);

    FRI_iop_protocol<FieldT> full_protocol(IOP,
                                           {FieldT::zero()},
                                           bcs_and_FRI_parameters.second);
    full_protocol.register_interactions();
    IOP.seal_interaction_registrations();
    full_protocol.register_queries();
    IOP.seal_query_registrations();

    full_protocol.produce_proof();

    libff::enter_block("Obtain transcript");
    libff::enter_block("Run verifier to populate virtual oracle data structures");
    full_protocol.verifier_predicate();
    libff::leave_block("Run verifier to populate virtual oracle data structures");

    const FRI_snark_proof<FieldT, hash_type> transcript = IOP.get_transcript();
    libff::leave_block("Obtain transcript");

    IOP.describe_sizes();

    libff::leave_block("FRI SNARK prover");
    return transcript;
}

template<typename FieldT, typename hash_type>
bool FRI_snark_verifier(const FRI_snark_proof<FieldT, hash_type> &proof,
                        const FRI_snark_parameters<FieldT> &parameters)
{
    libff::enter_block("FRI SNARK verifier");
    const std::pair<bcs_transformation_parameters<FieldT, hash_type>,
                    FRI_iop_protocol_parameters>
        bcs_and_FRI_parameters =
        obtain_bcs_and_FRI_parameters_from_FRI_snark_parameters<FieldT, hash_type>(parameters);

    bcs_transformation_parameters<FieldT, hash_type> bcs_parameters = bcs_and_FRI_parameters.first;
    bcs_verifier<FieldT, hash_type> IOP(bcs_parameters, proof);

    FRI_iop_protocol<FieldT> full_protocol(IOP,
                                           {FieldT::zero()},
                                           bcs_and_FRI_parameters.second);
    full_protocol.register_interactions();
    IOP.seal_interaction_registrations();
    full_protocol.register_queries();
    IOP.seal_query_registrations();

    libff::enter_block("Check semantic validity of IOP transcript");
    const bool IOP_transcript_valid = IOP.transcript_is_valid();
    libff::leave_block("Check semantic validity of IOP transcript");

    const bool full_protocol_accepts = full_protocol.verifier_predicate();

    libff::print_indent(); printf("* IOP transcript valid: %s\n", IOP_transcript_valid ? "true" : "false");
    libff::print_indent(); printf("* Full protocol decision predicate satisfied: %s\n", full_protocol_accepts ? "true" : "false");
    const bool decision = IOP_transcript_valid && full_protocol_accepts;
    libff::leave_block("FRI SNARK verifier");

    return decision;
}

template<typename FieldT, typename hash_type>
void FRI_snark_print_detailed_argument_size(
    FRI_snark_parameters<FieldT> params,
    FRI_snark_proof<FieldT, hash_type> argument)
{
    /* TODO: Lower all this boiler plate */
    const std::pair<bcs_transformation_parameters<FieldT, hash_type>,
                    FRI_iop_protocol_parameters>
        bcs_and_FRI_parameters =
        obtain_bcs_and_FRI_parameters_from_FRI_snark_parameters<FieldT, hash_type>(params);

    /* We go through registration on the verifier to know what the domains look like */
    bcs_verifier<FieldT, hash_type> verifier(bcs_and_FRI_parameters.first, argument);
    FRI_iop_protocol<FieldT> full_protocol(verifier,
                                           {FieldT::zero()},
                                           bcs_and_FRI_parameters.second);
    full_protocol.register_interactions();
    verifier.seal_interaction_registrations();
    full_protocol.register_queries();
    verifier.seal_query_registrations();
    const bool holographic = false;

    print_detailed_transcript_data<FieldT, hash_type>(
        holographic,
        argument,
        bcs_and_FRI_parameters.first,
        verifier);
}

} // namespace libiop
