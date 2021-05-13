#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>
#include "libiop/algebra/field_subset/subspace.hpp"
#include "libiop/bcs/bcs_common.hpp"
#include "libiop/bcs/common_bcs_parameters.hpp"

namespace libiop {


/** Initialize snark with FRI localization parameter array */
template<typename FieldT, typename hash_type>
fractal_snark_parameters<FieldT, hash_type>::fractal_snark_parameters(
    const size_t security_parameter,
    const LDT_reducer_soundness_type ldt_reducer_soundness_type,
    const FRI_soundness_type fri_soundness_type,
    const bcs_hash_type hash_enum,
    const std::vector<size_t> FRI_localization_parameter_array,
    const size_t RS_extra_dimensions,
    const bool make_zk,
    const field_subset_type domain_type,
    const std::shared_ptr<r1cs_constraint_system<FieldT>> constraint_system) :
    security_parameter_(security_parameter),
    LDT_reducer_soundness_type_(ldt_reducer_soundness_type),
    FRI_soundness_type_(fri_soundness_type),
    RS_extra_dimensions_(RS_extra_dimensions),
    make_zk_(make_zk),
    domain_type_(domain_type),
    constraint_system_(constraint_system),
    FRI_localization_parameter_array_(FRI_localization_parameter_array)
{
    this->initialize_bcs_params(hash_enum);
    this->initialize_iop_params();
}

/** Initialize snark with FRI localization parameter */
template<typename FieldT, typename hash_type>
fractal_snark_parameters<FieldT, hash_type>::fractal_snark_parameters(
    const size_t security_parameter,
    const LDT_reducer_soundness_type ldt_reducer_soundness_type,
    const FRI_soundness_type fri_soundness_type,
    const bcs_hash_type hash_enum,
    const size_t FRI_localization_parameter,
    const size_t RS_extra_dimensions,
    const bool make_zk,
    const field_subset_type domain_type,
    const std::shared_ptr<r1cs_constraint_system<FieldT>> constraint_system) :
    security_parameter_(security_parameter),
    LDT_reducer_soundness_type_(ldt_reducer_soundness_type),
    FRI_soundness_type_(fri_soundness_type),
    RS_extra_dimensions_(RS_extra_dimensions),
    make_zk_(make_zk),
    domain_type_(domain_type),
    constraint_system_(constraint_system),
    FRI_localization_parameter_(FRI_localization_parameter)
{
    this->initialize_bcs_params(hash_enum);
    this->initialize_iop_params();
}

template<typename FieldT, typename hash_type>
void fractal_snark_parameters<FieldT, hash_type>::reset_fri_localization_parameters(
    const std::vector<size_t> FRI_localization_parameter_array)
{
    this->FRI_localization_parameter_array_ = FRI_localization_parameter_array;
    this->initialize_iop_params();
}


template<typename FieldT, typename hash_type>
void fractal_snark_parameters<FieldT, hash_type>::initialize_iop_params()
{
    this->iop_params_ = fractal_iop_parameters<FieldT>(
        this->security_parameter_,
        this->bcs_params_.pow_params_.work_parameter(),
        this->RS_extra_dimensions_,
        this->make_zk_,
        this->constraint_system_);
    if (this->FRI_localization_parameter_array_.size() == 0) {
        this->iop_params_.set_ldt_parameters(this->FRI_localization_parameter_,
                                             this->FRI_soundness_type_,
                                             this->LDT_reducer_soundness_type_);
    } else {
        this->iop_params_.set_ldt_parameters(this->FRI_localization_parameter_array_,
                                             this->FRI_soundness_type_,
                                             this->LDT_reducer_soundness_type_);
        this->FRI_localization_parameter_array_ =
            this->iop_params_.FRI_params_.get_localization_parameters();
    }
}

template<typename FieldT, typename hash_type>
void fractal_snark_parameters<FieldT, hash_type>::initialize_bcs_params(const bcs_hash_type hash_enum)
{
    this->bcs_params_ = default_bcs_params<FieldT, hash_type>(hash_enum, this->security_parameter_, 
        libff::log2(this->constraint_system_->num_constraints()));
}

template<typename FieldT, typename hash_type>
void fractal_snark_parameters<FieldT, hash_type>::print() const
{
    libff::print_indent(); printf("\nFractal SNARK parameters\n");
    libff::print_indent(); printf("* security parameter (bits) = %zu\n", security_parameter_);
    libff::print_indent(); printf("* RS extra dimensions = %zu\n", RS_extra_dimensions_);
    libff::print_indent(); printf("* LDT reducer soundness type = %s\n",
        LDT_reducer_soundness_type_to_string(LDT_reducer_soundness_type_));
    libff::print_indent(); printf("* FRI soundness type = %s\n",
        FRI_soundness_type_to_string(FRI_soundness_type_));
    libff::print_indent(); printf("* zero-knowledge = %s\n", make_zk_ ? "true" : "false");
    libff::print_indent(); printf("* domain type = %s\n", field_subset_type_names[this->domain_type_]);

    this->iop_params_.print();
}

template<typename FieldT, typename hash_type>
std::pair<bcs_prover_index<FieldT, hash_type>, bcs_verifier_index<FieldT, hash_type>>
fractal_snark_indexer(
    const fractal_snark_parameters<FieldT, hash_type> &parameters)
{
    libff::enter_block("Fractal SNARK indexer");
    parameters.print();
    bcs_indexer<FieldT, hash_type> IOP(parameters.bcs_params_);
    fractal_iop<FieldT> full_protocol(IOP, parameters.iop_params_);
    IOP.seal_interaction_registrations();
    IOP.seal_query_registrations();
    full_protocol.produce_index();

    bcs_prover_index<FieldT, hash_type> prover_index = IOP.get_bcs_prover_index();
    bcs_verifier_index<FieldT, hash_type> verifier_index = IOP.get_verifier_index();
    std::pair<bcs_prover_index<FieldT, hash_type>, bcs_verifier_index<FieldT, hash_type>> index =
        std::make_pair(std::move(prover_index), verifier_index);
    libff::leave_block("Fractal SNARK indexer");
    return index;
}

template<typename FieldT, typename hash_type>
fractal_snark_argument<FieldT, hash_type> fractal_snark_prover(
    bcs_prover_index<FieldT, hash_type> &index,
    const r1cs_primary_input<FieldT> &primary_input,
    const r1cs_auxiliary_input<FieldT> &auxiliary_input,
    const fractal_snark_parameters<FieldT, hash_type> &parameters)
{
    libff::enter_block("Fractal SNARK prover");
    parameters.print();

    bcs_prover<FieldT, hash_type> IOP(parameters.bcs_params_, index);
    fractal_iop<FieldT> full_protocol(IOP, parameters.iop_params_);
    full_protocol.register_interactions();
    IOP.seal_interaction_registrations();
    full_protocol.register_queries();
    IOP.seal_query_registrations();

    full_protocol.produce_proof(primary_input, auxiliary_input, index.iop_index_);

    libff::enter_block("Obtain transcript");
    const fractal_snark_argument<FieldT, hash_type> transcript = IOP.get_transcript();
    libff::leave_block("Obtain transcript");

    IOP.describe_sizes();

    libff::leave_block("Fractal SNARK prover");
    return transcript;
}

template<typename FieldT, typename hash_type>
bool fractal_snark_verifier(
    const bcs_verifier_index<FieldT, hash_type> &index,
    const r1cs_primary_input<FieldT> &primary_input,
    const fractal_snark_argument<FieldT, hash_type> &proof,
    const fractal_snark_parameters<FieldT, hash_type> &parameters)
{
    libff::enter_block("Fractal SNARK verifier");
    parameters.print();

    bcs_verifier<FieldT, hash_type> IOP(parameters.bcs_params_, proof, index);

    libff::enter_block("Fractal IOP protocol initialization and registration");
    fractal_iop<FieldT> full_protocol(IOP, parameters.iop_params_);
    full_protocol.register_interactions();
    IOP.seal_interaction_registrations();
    full_protocol.register_queries();
    IOP.seal_query_registrations();
    libff::leave_block("Fractal IOP protocol initialization and registration");

    libff::enter_block("Check semantic validity of IOP transcript");
    const bool IOP_transcript_valid = IOP.transcript_is_valid();
    libff::leave_block("Check semantic validity of IOP transcript");

    const bool full_protocol_accepts = full_protocol.verifier_predicate(primary_input);

    libff::print_indent(); printf("* IOP transcript valid: %s\n", IOP_transcript_valid ? "true" : "false");
    libff::print_indent(); printf("* Full protocol decision predicate satisfied: %s\n", full_protocol_accepts ? "true" : "false");
    const bool decision = IOP_transcript_valid && full_protocol_accepts;
    libff::leave_block("Fractal SNARK verifier");

    return decision;
}

} // namespace libiop
