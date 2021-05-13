#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

#ifndef CPPDEBUG /* Ubuntu's Boost does not provide binaries compatible with libstdc++'s debug mode so we just reduce functionality here */
#include <boost/program_options.hpp>
#endif

#include "boost_profile.cpp"
#include <libff/algebra/curves/edwards/edwards_pp.hpp>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>

#include <libff/algebra/field_utils/field_utils.hpp>
#include <libff/algebra/fields/binary/gf64.hpp>
#include <libff/algebra/fields/binary/gf256.hpp>
#include <libff/algebra/fields/binary/gf192.hpp>
#include <libff/algebra/fields/binary/gf256.hpp>

#include "libiop/algebra/field_subset/subgroup.hpp"
#include "libiop/algebra/fft.hpp"

#include <libff/common/utils.hpp>
#include "libiop/iop/iop.hpp"
#include "libiop/protocols/ldt/fri/fri_ldt.hpp"
#include "libiop/snark/fri_snark.hpp"

#ifndef CPPDEBUG
bool process_prover_command_line(const int argc, const char** argv, options &options,
                                 std::size_t localization_parameter,
                                 std::size_t num_localization_steps,
                                 std::size_t num_oracles,
                                 std::size_t num_interactive_repetitions,
                                 std::size_t num_query_repetitions)
{
    namespace po = boost::program_options;

    try
    {
        po::options_description desc = gen_options(options);
        desc.add_options()
            ("localization_parameter", po::value<std::size_t>(&localization_parameter)->default_value(2), "Only used when num_localization_steps is 0")
            ("num_localization_steps", po::value<std::size_t>(&num_localization_steps)->default_value(0))
            ("num_oracles", po::value<std::size_t>(&num_oracles)->default_value(1))
            ("num_interactive_repetitions", po::value<std::size_t>(&num_interactive_repetitions)->default_value(1))
            ("num_query_repetitions", po::value<std::size_t>(&num_query_repetitions)->default_value(64));

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

template<typename FieldT, typename hash_type>
void instrument_FRI(options &options,
                    std::size_t localization_parameter,
                    std::size_t num_localization_steps,
                    std::size_t num_oracles,
                    std::size_t num_interactive_repetitions,
                    std::size_t num_query_repetitions)
{
    for (std::size_t log_n = options.log_n_min; log_n <= options.log_n_max; ++log_n)
    {
        libff::print_separator();
        const std::size_t poly_degree_bound = 1ull << log_n;
        const std::size_t RS_extra_dimensions = 2; /* \rho = 2^{-RS_extra_dimensions} */
        const std::size_t codeword_domain_dim = log_n + RS_extra_dimensions;


        std::vector<std::size_t> localization_parameter_array;
        if (num_localization_steps != 0)
        {
            std::size_t remaining = codeword_domain_dim - RS_extra_dimensions - 1;
            std::size_t vals = remaining / num_localization_steps;
            localization_parameter_array = std::vector<std::size_t>(num_localization_steps, vals);
            localization_parameter_array.insert(localization_parameter_array.begin(), 1);
        }

        std::cout << "Codeword domain dimension: " << codeword_domain_dim << "\n"
                << "RS_extra_dimensions: " << RS_extra_dimensions << "\n"
                << "poly_degree_bound: " << poly_degree_bound << "\n"
                << "\n";

        const polynomial<FieldT> poly = polynomial<FieldT>::random_polynomial(poly_degree_bound);

        /* Set up the protocol blueprint */
        iop_protocol<FieldT> IOP;

        const std::size_t codeword_domain_size = 1ull << codeword_domain_dim;

        FRI_snark_parameters<FieldT> params;
        params.codeword_domain_dim_ = codeword_domain_dim;
        params.security_level_ = options.security_level;
        params.hash_enum_ = options.hash_enum;
        params.RS_extra_dimensions_ = RS_extra_dimensions;
        params.localization_parameter_array_ = localization_parameter_array;
        params.localization_parameter_ = localization_parameter;
        params.num_interactive_repetitions_ = num_interactive_repetitions;
        params.num_query_repetitions_ = num_query_repetitions;
        params.field_type_ = libff::get_field_type<FieldT>(FieldT::zero());
        params.num_oracles_ = num_oracles;

        const FRI_snark_proof<FieldT, hash_type> proof = FRI_snark_prover<FieldT, hash_type>(params);
        printf("\n");
        libff::print_indent(); printf("* Argument size in bytes (IOP): %zu\n", proof.IOP_size_in_bytes());
        libff::print_indent(); printf("* Argument size in bytes (BCS): %zu\n", proof.BCS_size_in_bytes());
        libff::print_indent(); printf("* Argument size in bytes (total): %zu\n", proof.size_in_bytes());

        printf("\nIf we were to remove the pruning of BCS merkle tree paths feature,\n"
               "the argument would have the following sizes:\n");
        libff::print_indent(); printf("* Argument size in bytes (BCS, no pruning): %zu\n", proof.BCS_size_in_bytes_without_pruning());
        libff::print_indent(); printf("* Argument size in bytes (total, no pruning): %zu\n", proof.size_in_bytes_without_pruning());
        printf("\n");

        const bool bit = FRI_snark_verifier<FieldT, hash_type>(proof, params);

        libff::print_indent(); printf("* Verifier satisfied: %s\n", bit ? "true" : "false");
    }
}

int main(int argc, const char * argv[])
{
    options default_vals;

    std::size_t localization_parameter = 2;
    std::size_t num_localization_steps = 0;
    std::size_t num_oracles = 1;
    std::size_t num_interactive_repetitions = 1;
    std::size_t num_query_repetitions = 10;

#ifdef CPPDEBUG
    /* set reasonable defaults */
    if (argc > 1)
    {
        printf("There is no argument parsing in CPPDEBUG mode.");
        exit(1);
    }
    libff::UNUSED(argv);

#else
    if (!process_prover_command_line(argc, argv, default_vals, localization_parameter,
                                     num_interactive_repetitions, num_query_repetitions,
                                     num_localization_steps, num_oracles))
    {
        return 1;
    }

#endif

    libff::start_profiling();

    printf("Selected parameters:\n");
    printf("* log_n_min = %zu\n", default_vals.log_n_min);
    printf("* log_n_max = %zu\n", default_vals.log_n_max);
    printf("* security_level = %zu\n", default_vals.security_level);

    if (default_vals.is_multiplicative) {
        switch (default_vals.field_size) {
            case 181:
                libff::edwards_pp::init_public_params();
                instrument_FRI<libff::edwards_Fr, binary_hash_digest>(default_vals, localization_parameter,
                                     num_interactive_repetitions, num_query_repetitions,
                                     num_localization_steps, num_oracles);
                break;
            case 256:
                libff::alt_bn128_pp::init_public_params();
                instrument_FRI<libff::alt_bn128_Fr, binary_hash_digest>(default_vals, localization_parameter,
                                     num_interactive_repetitions, num_query_repetitions,
                                     num_localization_steps, num_oracles);
                break;
            default:
                throw std::invalid_argument("Field size not supported.");
        }
    } else {
        switch (default_vals.field_size)
        {
            case 64:
                instrument_FRI<libff::gf64, binary_hash_digest>(default_vals, localization_parameter,
                                     num_interactive_repetitions, num_query_repetitions,
                                     num_localization_steps, num_oracles);
                break;
            case 128:
                instrument_FRI<libff::gf128, binary_hash_digest>(default_vals, localization_parameter,
                                     num_interactive_repetitions, num_query_repetitions,
                                     num_localization_steps, num_oracles);
                break;
            case 192:
                instrument_FRI<libff::gf192, binary_hash_digest>(default_vals, localization_parameter,
                                     num_interactive_repetitions, num_query_repetitions,
                                     num_localization_steps, num_oracles);
                break;
            case 256:
                instrument_FRI<libff::gf256, binary_hash_digest>(default_vals, localization_parameter,
                                     num_interactive_repetitions, num_query_repetitions,
                                     num_localization_steps, num_oracles);
                break;
            default:
                throw std::invalid_argument("Field size not supported.");
        }
    }
}
