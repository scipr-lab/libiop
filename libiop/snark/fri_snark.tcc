#include "libiop/algebra/field_subset/subspace.hpp"
#include "libiop/common/common.hpp"
#include "libiop/common/profiling.hpp"
#include "libiop/protocols/fri_iop.hpp"
#include "libiop/snark/common/bcs16_common.hpp"

namespace libiop {

template<typename FieldT>
void FRI_snark_parameters<FieldT>::describe()
{
    print_indent(); printf("* FRI SNARK parameters:\n");
    print_indent(); printf("* Security level: %zu\n", security_level_);
    print_indent(); printf("* RS extra dimensions: %zu\n", RS_extra_dimensions_);
    print_indent(); printf("* Localization parameter: %zu\n", localization_parameter_);
    print_indent(); printf("* Localization parameter array: %zu\n", localization_parameter_array_);
    print_indent(); printf("* Num query repetitions: %zu\n", num_query_repetitions_);
}

template<typename FieldT>
const std::pair<bcs16_transformation_parameters<FieldT>,
                FRI_iop_protocol_parameters>
obtain_bcs16_and_FRI_parameters_from_FRI_snark_parameters(const FRI_snark_parameters<FieldT> &parameters)
{
    bcs16_transformation_parameters<FieldT> bcs_parameters;
    bcs_parameters.security_parameter = parameters.security_level_;
    bcs_parameters.field_hasher = blake2b_field_element_hash<FieldT>;
    bcs_parameters.zk_hasher = blake2b_zk_element_hash;
    bcs_parameters.compression_hasher = blake2b_two_to_one_hash;
    bcs_parameters.FieldT_randomness_extractor = blake2b_FieldT_randomness_extractor<FieldT>;
    bcs_parameters.integer_randomness_extractor = blake2b_integer_randomness_extractor;

    FRI_iop_protocol_parameters FRI_parameters;
    FRI_parameters.RS_extra_dimensions_ = parameters.RS_extra_dimensions_;
    FRI_parameters.codeword_domain_dim_ = parameters.codeword_domain_dim_;
    FRI_parameters.localization_parameter_ = parameters.localization_parameter_;
    FRI_parameters.localization_parameter_array_ = parameters.localization_parameter_array_;
    FRI_parameters.num_interactive_repetitions_ = parameters.num_interactive_repetitions_;
    FRI_parameters.num_query_repetitions_ = parameters.num_query_repetitions_;
    FRI_parameters.field_type_ = parameters.field_type_;
    FRI_parameters.num_oracles_ = parameters.num_oracles_;

    return std::make_pair(bcs_parameters, FRI_parameters);
}

template<typename FieldT>
FRI_snark_proof<FieldT> FRI_snark_prover(const std::vector<FieldT> evaluations,
                                         const FRI_snark_parameters<FieldT> &parameters)
{
    enter_block("FRI SNARK prover");
    const std::pair<bcs16_transformation_parameters<FieldT>,
                    FRI_iop_protocol_parameters>
        bcs16_and_FRI_parameters =
        obtain_bcs16_and_FRI_parameters_from_FRI_snark_parameters<FieldT>(parameters);

    bcs16_transformation_parameters<FieldT> bcs16_parameters = bcs16_and_FRI_parameters.first;
    bcs16_prover<FieldT> IOP(bcs16_parameters);

    FRI_iop_protocol<FieldT> full_protocol(IOP,
                                           evaluations,
                                           bcs16_and_FRI_parameters.second);
    full_protocol.register_interactions();
    IOP.seal_interaction_registrations();
    full_protocol.register_queries();
    IOP.seal_query_registrations();

    full_protocol.produce_proof();

    enter_block("Obtain transcript");
    enter_block("Run verifier to populate virtual oracle data structures");
    full_protocol.verifier_predicate();
    leave_block("Run verifier to populate virtual oracle data structures");

    const FRI_snark_proof<FieldT> transcript = IOP.get_transcript();
    leave_block("Obtain transcript");

    IOP.describe_sizes();

    leave_block("FRI SNARK prover");
    return transcript;
}

template<typename FieldT>
bool FRI_snark_verifier(const std::vector<FieldT> evaluations,
                        const FRI_snark_proof<FieldT> &proof,
                        const FRI_snark_parameters<FieldT> &parameters)
{
    enter_block("FRI SNARK verifier");
    const std::pair<bcs16_transformation_parameters<FieldT>,
                    FRI_iop_protocol_parameters>
        bcs16_and_FRI_parameters =
        obtain_bcs16_and_FRI_parameters_from_FRI_snark_parameters<FieldT>(parameters);

    bcs16_transformation_parameters<FieldT> bcs16_parameters = bcs16_and_FRI_parameters.first;
    bcs16_verifier<FieldT> IOP(bcs16_parameters, proof);

    FRI_iop_protocol<FieldT> full_protocol(IOP,
                                           evaluations,
                                           bcs16_and_FRI_parameters.second);
    full_protocol.register_interactions();
    IOP.seal_interaction_registrations();
    full_protocol.register_queries();
    IOP.seal_query_registrations();

    enter_block("Check semantic validity of IOP transcript");
    const bool IOP_transcript_valid = IOP.transcript_is_valid();
    leave_block("Check semantic validity of IOP transcript");

    const bool full_protocol_accepts = full_protocol.verifier_predicate();

    print_indent(); printf("* IOP transcript valid: %s\n", IOP_transcript_valid ? "true" : "false");
    print_indent(); printf("* Full protocol decision predicate satisfied: %s\n", full_protocol_accepts ? "true" : "false");
    const bool decision = IOP_transcript_valid && full_protocol_accepts;
    leave_block("FRI SNARK verifier");

    return decision;
}

} // namespace libiop
