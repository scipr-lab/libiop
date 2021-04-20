#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>


#ifndef CPPDEBUG /* Ubuntu's Boost does not provide binaries compatible with libstdc++'s debug mode so we just reduce functionality here */
#include <boost/program_options.hpp>
#endif

#include "boost_profile.cpp"
#include <libff/algebra/fields/binary/gf64.hpp>
#include <libff/algebra/fields/binary/gf256.hpp>
#include <libff/algebra/fields/binary/gf192.hpp>
#include <libff/algebra/fields/binary/gf256.hpp>
#include <libff/algebra/curves/edwards/edwards_pp.hpp>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/field_utils/field_utils.hpp>

#include "libiop/snark/fractal_snark.hpp"
#include "libiop/bcs/bcs_common.hpp"
#include "libiop/protocols/fractal_hiop.hpp"
#include "libiop/protocols/ldt/fri/argument_size_optimizer.hpp"
#include "libiop/relations/examples/r1cs_examples.hpp"

#ifndef CPPDEBUG
bool process_prover_command_line(const int argc, const char** argv,
                                 options &options, bool heuristic_fri_soundness, bool optimize_localization)
{
    namespace po = boost::program_options;

    try
    {
        po::options_description desc = gen_options(options);
        desc.add_options()
             ("optimize_localization", po::value<bool>(&optimize_localization)->default_value(false))
             ("heuristic_fri_soundness", po::value<bool>(&heuristic_fri_soundness)->default_value(true));

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);


        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return false;
        }        

        po::notify(vm);
        options.hash_enum = static_cast<libiop::bcs_hash_type>(options.hash_enum_val);
    }
    catch(std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return false;
    }

    return true;
}
#endif

using namespace libiop;

template<typename FieldT, typename hash_type>
void print_argument_size(
    fractal_snark_parameters<FieldT, hash_type> params,
    bcs_verifier_index<FieldT, hash_type> index,
    fractal_snark_argument<FieldT, hash_type> argument)
{
    /* We go through registration on the verifier to know what the domains look like */
    bcs_verifier<FieldT, hash_type> verifier(params.bcs_params_, argument, index);
    fractal_iop<FieldT> full_protocol(verifier, params.iop_params_);
    full_protocol.register_interactions();
    verifier.seal_interaction_registrations();
    full_protocol.register_queries();
    verifier.seal_query_registrations();
    const bool holographic = true;
    print_detailed_transcript_data<FieldT>(
        holographic,
        argument,
        params.bcs_params_,
        verifier);
}

template<typename FieldT, typename hash_type>
void instrument_fractal_snark(options &options,
                              LDT_reducer_soundness_type ldt_reducer_soundness_type,
                              FRI_soundness_type fri_soundness_type,
                              bool &optimize_localization)
{
    // TODO: Unhard code this
    const size_t RS_extra_dimensions = 3;
    const size_t fri_localization_parameter = 2;
    field_subset_type domain_type = affine_subspace_type;
    if (options.is_multiplicative) {
        domain_type = multiplicative_coset_type;
    }

    for (std::size_t log_n = options.log_n_min; log_n <= options.log_n_max; ++log_n)
    {
        libff::print_separator();

        const std::size_t n = 1ul << log_n;
        /* k+1 needs to be a power of 2 (proof system artifact) so we just fix it to 15 here */
        size_t k = 15;
        if (domain_type == multiplicative_coset_type)
        {
            k = 0;
        }
        const std::size_t m = n - 1;
        r1cs_example<FieldT> example = generate_r1cs_example<FieldT>(n, k, m);

        fractal_snark_parameters<FieldT, hash_type> parameters(
            options.security_level,
            ldt_reducer_soundness_type,
            fri_soundness_type,
            options.hash_enum,
            fri_localization_parameter,
            RS_extra_dimensions,
            options.make_zk,
            domain_type,
            std::make_shared<r1cs_constraint_system<FieldT>>(example.constraint_system_));

        std::vector<std::size_t> localization_parameter_array;
        if (optimize_localization)
        {
            const size_t codeword_dim = parameters.iop_params_.codeword_domain().dimension();
            size_t num_query_sets = parameters.iop_params_.FRI_params_.query_repetitions();
            size_t interactive_repetitions =
                parameters.iop_params_.FRI_params_.interactive_repetitions() *
                parameters.iop_params_.LDT_reducer_params_.num_output_LDT_instances();
            const size_t hash_size = (parameters.bcs_params_.security_parameter + 3) / 4;
            std::vector<size_t> oracle_locality_vector = parameters.iop_params_.locality_vector();
            // if (parameters.iop_params_.make_zk())
            // {
            //     /* Handle the zk leaves for the SNARK. TODO: Where should this go? */
            //     oracle_locality_vector[1] += 1;
            // }

            /* TODO: Get exact max tested degree */
            const size_t max_tested_degree = 6 * parameters.iop_params_.index_domain().num_elements();
            localization_parameter_array =
                compute_argument_size_optimal_localization_parameters<FieldT>(
                    oracle_locality_vector, codeword_dim,
                    num_query_sets, interactive_repetitions,
                    max_tested_degree, hash_size);

            parameters.reset_fri_localization_parameters(localization_parameter_array);
        }

        libff::enter_block("Check satisfiability of R1CS example");
        const bool is_satisfied = example.constraint_system_.is_satisfied(
            example.primary_input_, example.auxiliary_input_);
        assert(is_satisfied);
        libff::leave_block("Check satisfiability of R1CS example");
        printf("\n");
        libff::print_indent(); printf("* R1CS number of constraints: %zu\n", example.constraint_system_.num_constraints());
        libff::print_indent(); printf("* R1CS number of variables: %zu\n", example.constraint_system_.num_variables());
        libff::print_indent(); printf("* R1CS number of variables for primary input: %zu\n", example.primary_input_.size());
        libff::print_indent(); printf("* R1CS number of variables for auxiliary input: %zu\n", example.auxiliary_input_.size());
        libff::print_indent(); printf("* R1CS size of constraint system (bytes): %zu\n", example.constraint_system_.size_in_bytes());
        libff::print_indent(); printf("* R1CS size of primary input (bytes): %zu\n", example.primary_input_.size() * sizeof(FieldT));
        libff::print_indent(); printf("* R1CS size of auxiliary input (bytes): %zu\n", example.auxiliary_input_.size() * sizeof(FieldT));
        printf("\n");

        std::pair<bcs_prover_index<FieldT, hash_type>, bcs_verifier_index<FieldT, hash_type>> index =
            fractal_snark_indexer(parameters);

        /** TODO: Print some useful data about the indexed data */

        const fractal_snark_argument<FieldT, hash_type> argument =
            fractal_snark_prover(
                index.first,
                example.primary_input_,
                example.auxiliary_input_,
                parameters);

        parameters = fractal_snark_parameters<FieldT, hash_type>(
            options.security_level,
            ldt_reducer_soundness_type,
            fri_soundness_type,
            options.hash_enum,
            fri_localization_parameter,
            RS_extra_dimensions,
            options.make_zk,
            domain_type,
            std::make_shared<r1cs_constraint_system<FieldT>>(example.constraint_system_));
        if (optimize_localization)
        {
            parameters.reset_fri_localization_parameters(localization_parameter_array);
        }

        print_argument_size(parameters, index.second, argument);

        const bool bit = fractal_snark_verifier<FieldT, hash_type>(
            index.second,
            example.primary_input_,
            argument,
            parameters);

        printf("\n\n");

        libff::print_indent(); printf("* Verifier satisfied: %s\n", bit ? "true" : "false");
    }
}

int main(int argc, const char * argv[])
{
    /* Set up R1CS */
    options default_vals;

    bool optimize_localization = false;
    bool heuristic_fri_soundness = true;

#ifdef CPPDEBUG
    /* set reasonable defaults */
    if (argc > 1)
    {
        printf("There is no argument parsing in CPPDEBUG mode.");
        exit(1);
    }
    libff::UNUSED(argv);

#else
    if (!process_prover_command_line(argc, argv, default_vals, heuristic_fri_soundness, optimize_localization))
    {
        return 1;
    }
#endif
    /** TODO: eventually get a string from program options, and then have a from string methods in protocols */
    LDT_reducer_soundness_type ldt_reducer_soundness_type = LDT_reducer_soundness_type::proven;
    if (default_vals.heuristic_ldt_reducer_soundness)
    {
        ldt_reducer_soundness_type = LDT_reducer_soundness_type::optimistic_heuristic;
    }
    FRI_soundness_type fri_soundness_type = FRI_soundness_type::proven;
    if (heuristic_fri_soundness) {
        fri_soundness_type = FRI_soundness_type::heuristic;
    }
    libff::start_profiling();

    printf("Selected parameters:\n");
    printf("- log_n_min = %zu\n", default_vals.log_n_min);
    printf("- log_n_max = %zu\n", default_vals.log_n_max);
    printf("- security_level = %zu\n", default_vals.security_level);
    printf("- LDT_reducer_soundness_type = %s\n", LDT_reducer_soundness_type_to_string(ldt_reducer_soundness_type));
    printf("- FRI_soundness_type = %s\n", FRI_soundness_type_to_string(fri_soundness_type));
    printf("- is_multiplicative = %s\n", default_vals.is_multiplicative ? "true" : "false");
    printf("- field_size = %zu\n", default_vals.field_size);
    printf("- make_zk = %s\n", default_vals.make_zk ? "true" : "false");
    printf("- hash_enum = %s\n", bcs_hash_type_names[default_vals.hash_enum]);
    
    if (default_vals.is_multiplicative) {
        switch (default_vals.field_size) {
            case 181:
                libff::edwards_pp::init_public_params();
                instrument_fractal_snark<libff::edwards_Fr, binary_hash_digest>(
                    default_vals, ldt_reducer_soundness_type, fri_soundness_type, optimize_localization);
                break;
            case 256:
                libff::alt_bn128_pp::init_public_params();
                if (default_vals.hash_enum == libiop::blake2b_type)
                {
                    instrument_fractal_snark<libff::alt_bn128_Fr, binary_hash_digest>(
                        default_vals, ldt_reducer_soundness_type, fri_soundness_type, optimize_localization);
                }
                else
                {
                    instrument_fractal_snark<libff::alt_bn128_Fr, libff::alt_bn128_Fr>(
                        default_vals, ldt_reducer_soundness_type, fri_soundness_type, optimize_localization);
                }
                break;
            default:
                throw std::invalid_argument("Field size not supported.");
        }

    } else {
        switch (default_vals.field_size)
        {
            case 64:
                instrument_fractal_snark<libff::gf64, binary_hash_digest>(
                    default_vals, ldt_reducer_soundness_type, fri_soundness_type, optimize_localization);
                break;
            case 128:
                instrument_fractal_snark<libff::gf128, binary_hash_digest>(
                    default_vals, ldt_reducer_soundness_type, fri_soundness_type, optimize_localization);
                break;
            case 192:
                instrument_fractal_snark<libff::gf192, binary_hash_digest>(
                    default_vals, ldt_reducer_soundness_type, fri_soundness_type, optimize_localization);
                break;
            case 256:
                instrument_fractal_snark<libff::gf256, binary_hash_digest>(
                    default_vals, ldt_reducer_soundness_type, fri_soundness_type, optimize_localization);
                break;
            default:
                throw std::invalid_argument("Field size not supported.");
        }
    }
}
