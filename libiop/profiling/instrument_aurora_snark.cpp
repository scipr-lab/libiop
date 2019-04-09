#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>


#ifndef CPPDEBUG /* Ubuntu's Boost does not provide binaries compatible with libstdc++'s debug mode so we just reduce functionality here */
#include <boost/program_options.hpp>
#endif

#include "libiop/algebra/fields/gf64.hpp"
#include "libiop/algebra/fields/gf128.hpp"
#include "libiop/algebra/fields/gf192.hpp"
#include "libiop/algebra/fields/gf256.hpp"
#include "libiop/algebra/fields/utils.hpp"

#include "libiop/snark/aurora_snark.hpp"
#include "libiop/protocols/aurora_iop.hpp"
#include "libiop/protocols/ldt/fri/optimizer.hpp"
#include "libiop/relations/examples/r1cs_examples.hpp"
#include <libff/algebra/curves/edwards/edwards_pp.hpp>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>

#ifndef CPPDEBUG
bool process_prover_command_line(const int argc, const char** argv,
                                 std::size_t &log_n_min,
                                 std::size_t &log_n_max,
                                 std::size_t &security_level,
                                 std::size_t &field_size,
                                 bool &heuristic_ldt_reducer_soundness,
                                 bool &heuristic_fri_soundness,
                                 bool &make_zk,
                                 bool &is_multiplicative,
                                 bool &optimize_localization)
{
    namespace po = boost::program_options;

    try
    {
        po::options_description desc("Usage");
        desc.add_options()
            ("help", "print this help message")
            ("log_n_min", po::value<std::size_t>(&log_n_min)->default_value(8))
            ("log_n_max", po::value<std::size_t>(&log_n_max)->default_value(20))
            ("security_level", po::value<std::size_t>(&security_level)->default_value(128))
            ("field_size", po::value<std::size_t>(&field_size)->default_value(192))
            ("heuristic_ldt_reducer_soundness", po::value<bool>(&heuristic_ldt_reducer_soundness)->default_value(true))
            ("heuristic_fri_soundness", po::value<bool>(&heuristic_fri_soundness)->default_value(true))
            ("make_zk", po::value<bool>(&make_zk)->default_value(false))
            ("is_multiplicative", po::value<bool>(&is_multiplicative)->default_value(false))
            ("optimize_localization", po::value<bool>(&optimize_localization)->default_value(false));

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return false;
        }

        po::notify(vm);
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

template<typename FieldT>
void instrument_aurora_snark(const std::size_t log_n_min,
                             const std::size_t log_n_max,
                             std::size_t security_level,
                             LDT_reducer_soundness_type ldt_reducer_soundness_type,
                             FRI_soundness_type fri_soundness_type,
                             const bool make_zk,
                             const bool is_multiplicative,
                             bool optimize_localization)
{
    // TODO: Unhard code this
    const size_t RS_extra_dimensions = 3 + (make_zk ? 0 : 2);
    const size_t fri_localization_parameter = 2;
    field_subset_type domain_type = affine_subspace_type;
    if (is_multiplicative) {
        domain_type = multiplicative_coset_type;
    }

    for (std::size_t log_n = log_n_min; log_n <= log_n_max; ++log_n)
    {
        print_separator();

        const std::size_t n = 1ul << log_n;
        /* k+1 needs to be a power of 2 (proof system artifact) so we just fix it to 15 here */
        const std::size_t k = 15;
        const std::size_t m = (n <= 116 ? n : n - 101);
        r1cs_example<FieldT> example = generate_r1cs_example<FieldT>(n, k, m);


        aurora_snark_parameters<FieldT> parameters(security_level,
                                                   ldt_reducer_soundness_type,
                                                   fri_soundness_type,
                                                   fri_localization_parameter,
                                                   RS_extra_dimensions,
                                                   make_zk,
                                                   domain_type,
                                                   example.constraint_system_.num_constraints(),
                                                   example.constraint_system_.num_variables());

        std::vector<std::size_t> localization_parameter_array;
        if (optimize_localization)
        {
            std::size_t num_query_sets = parameters.iop_params_.FRI_params_.query_repetitions();
            std::size_t field_size = log_of_field_size_helper<FieldT>(FieldT(0));
            size_t locality = parameters.iop_params_.locality();
            localization_parameter_array = brute_force_optimal_localization_parameters(
                log_n, num_query_sets, locality, field_size, RS_extra_dimensions);

            parameters.reset_fri_localization_parameters(localization_parameter_array);
        }

        enter_block("Check satisfiability of R1CS example");
        const bool is_satisfied = example.constraint_system_.is_satisfied(
            example.primary_input_, example.auxiliary_input_);
        assert(is_satisfied);
        leave_block("Check satisfiability of R1CS example");
        printf("\n");
        print_indent(); printf("* R1CS number of constraints: %zu\n", example.constraint_system_.num_constraints());
        print_indent(); printf("* R1CS number of variables: %zu\n", example.constraint_system_.num_variables());
        print_indent(); printf("* R1CS number of variables for primary input: %zu\n", example.primary_input_.size());
        print_indent(); printf("* R1CS number of variables for auxiliary input: %zu\n", example.auxiliary_input_.size());
        print_indent(); printf("* R1CS size of constraint system (bytes): %zu\n", example.constraint_system_.size_in_bytes());
        print_indent(); printf("* R1CS size of primary input (bytes): %zu\n", example.primary_input_.size() * sizeof(FieldT));
        print_indent(); printf("* R1CS size of auxiliary input (bytes): %zu\n", example.auxiliary_input_.size() * sizeof(FieldT));
        printf("\n");
        const aurora_snark_argument<FieldT> proof = aurora_snark_prover<FieldT>(
            example.constraint_system_,
            example.primary_input_,
            example.auxiliary_input_,
            parameters);

        printf("\n");

        print_indent(); printf("* Argument size in bytes (IOP): %zu\n", proof.IOP_size_in_bytes());
        print_indent(); printf("* Argument size in bytes (BCS): %zu\n", proof.BCS_size_in_bytes());
        print_indent(); printf("* Argument size in bytes (total): %zu\n", proof.size_in_bytes());

        printf("\nIf we were to remove pruning of authentication paths in BCS,\n"
               "the argument would have the following sizes:\n");
        print_indent(); printf("* Argument size in bytes (BCS, no pruning): %zu\n", proof.BCS_size_in_bytes_without_pruning());
        print_indent(); printf("* Argument size in bytes (total, no pruning): %zu\n", proof.size_in_bytes_without_pruning());

        printf("\n");

        const bool bit = aurora_snark_verifier<FieldT>(
            example.constraint_system_,
            example.primary_input_,
            proof,
            parameters);

        printf("\n\n");

        print_indent(); printf("* Verifier satisfied: %s\n", bit ? "true" : "false");
    }
}

int main(int argc, const char * argv[])
{
    /* Set up R1CS */
    std::size_t log_n_min;
    std::size_t log_n_max;
    std::size_t security_level;
    std::size_t field_size;
    bool heuristic_ldt_reducer_soundness;
    bool heuristic_fri_soundness;
    bool make_zk;
    bool is_multiplicative;
    bool optimize_localization;

#ifdef CPPDEBUG
    /* set reasonable defaults */
    if (argc > 1)
    {
        printf("There is no argument parsing in CPPDEBUG mode.");
        exit(1);
    }
    libiop::UNUSED(argv);

    log_n_min = 8;
    log_n_max = 20;
    security_level = 128;
    field_size = 64;
    heuristic_ldt_reducer_soundness = true;
    heuristic_fri_soundness = true;
    make_zk = false;
    is_multiplicative = false;
    optimize_localization = false;
#else
    if (!process_prover_command_line(argc, argv, log_n_min, log_n_max, security_level, field_size,
        heuristic_ldt_reducer_soundness, heuristic_fri_soundness, make_zk, is_multiplicative, optimize_localization))
    {
        return 1;
    }
#endif
    /** TODO: eventually get a string from program options, and then have a from string methods in protocols */
    LDT_reducer_soundness_type ldt_reducer_soundness_type = LDT_reducer_soundness_type::proven;
    if (heuristic_ldt_reducer_soundness)
    {
        ldt_reducer_soundness_type = LDT_reducer_soundness_type::optimistic_heuristic;
    }
    FRI_soundness_type fri_soundness_type = FRI_soundness_type::proven;
    if (heuristic_fri_soundness) {
        fri_soundness_type = FRI_soundness_type::heuristic;
    }
    start_profiling();

    printf("Selected parameters:\n");
    printf("- log_n_min = %zu\n", log_n_min);
    printf("- log_n_max = %zu\n", log_n_max);
    printf("- security_level = %zu\n", security_level);
    printf("- LDT_reducer_soundness_type = %s\n", LDT_reducer_soundness_type_to_string(ldt_reducer_soundness_type));
    printf("- FRI_soundness_type = %s\n", FRI_soundness_type_to_string(fri_soundness_type));
    printf("- is_multiplicative = %s\n", is_multiplicative ? "true" : "false");
    printf("- field_size = %zu\n", field_size);
    printf("- make_zk = %s\n", make_zk ? "true" : "false");

    if (is_multiplicative) {
        switch (field_size) {
            case 181:
                edwards_pp::init_public_params();
                instrument_aurora_snark<edwards_Fr>(log_n_min, log_n_max, security_level,
                    ldt_reducer_soundness_type, fri_soundness_type, make_zk, is_multiplicative, optimize_localization);
                break;
            case 256:
                libff::alt_bn128_pp::init_public_params();
                instrument_aurora_snark<libff::alt_bn128_Fr>(log_n_min, log_n_max, security_level,
                    ldt_reducer_soundness_type, fri_soundness_type, make_zk, is_multiplicative, optimize_localization);
                break;
            default:
                throw std::invalid_argument("Field size not supported.");
        }

    } else {
        switch (field_size)
        {
            case 64:
                instrument_aurora_snark<gf64>(log_n_min, log_n_max, security_level,
                    ldt_reducer_soundness_type, fri_soundness_type, make_zk, is_multiplicative, optimize_localization);
                break;
            case 128:
                instrument_aurora_snark<gf128>(log_n_min, log_n_max, security_level,
                    ldt_reducer_soundness_type, fri_soundness_type, make_zk, is_multiplicative, optimize_localization);
                break;
            case 192:
                instrument_aurora_snark<gf192>(log_n_min, log_n_max, security_level,
                    ldt_reducer_soundness_type, fri_soundness_type, make_zk, is_multiplicative, optimize_localization);
                break;
            case 256:
                instrument_aurora_snark<gf256>(log_n_min, log_n_max, security_level,
                    ldt_reducer_soundness_type, fri_soundness_type, make_zk, is_multiplicative, optimize_localization);
                break;
            default:
                throw std::invalid_argument("Field size not supported.");
        }
    }
}
