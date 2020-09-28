#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

#ifndef CPPDEBUG /* Ubuntu's Boost does not provide binaries compatible with libstdc++'s debug mode so we just reduce functionality here */
#include <boost/program_options.hpp>
#endif

#include "snark_types.hpp"
#include "boost_profile.cpp"
#include "libff/algebra/curves/edwards/edwards_pp.hpp"
#include "libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp"

#include "libiop/algebra/fields/utils.hpp"
#include "libiop/algebra/fft.hpp"
#include "libiop/algebra/fields/gf64.hpp"
#include "libiop/algebra/fields/gf128.hpp"
#include "libiop/algebra/fields/gf192.hpp"
#include "libiop/algebra/fields/gf256.hpp"
#include "libiop/algebra/field_subset/subgroup.hpp"

#include "libiop/common/common.hpp"
#include "libiop/iop/iop.hpp"
#include "libiop/protocols/ldt/fri/fri_ldt.hpp"
#include "libiop/snark/fri_snark.hpp"

#ifndef CPPDEBUG
bool process_prover_command_line(const int argc, const char** argv, options &options)
{
    namespace po = boost::program_options;

    try
    {
        po::options_description desc = gen_options(options);
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
void instrument_FRI(options &options)
{
    for (std::size_t log_n = options.log_n_min; log_n <= options.log_n_max; ++log_n)
    {
        print_separator();
        const std::size_t poly_degree_bound = 1ull << log_n;
        const std::size_t RS_extra_dimensions = 2; /* \rho = 2^{-RS_extra_dimensions} */
        const std::size_t codeword_domain_dim = log_n + RS_extra_dimensions;


        std::vector<std::size_t> localization_parameter_array;
        if (options.num_localization_steps != 0)
        {
            std::size_t remaining = codeword_domain_dim - RS_extra_dimensions - 1;
            std::size_t vals = remaining / options.num_localization_steps;
            localization_parameter_array = std::vector<std::size_t>(options.num_localization_steps, vals);
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
        params.RS_extra_dimensions_ = options.RS_extra_dimensions;
        params.localization_parameter_array_ = localization_parameter_array;
        params.localization_parameter_ = options.localization_parameter;
        params.num_interactive_repetitions_ = options.num_interactive_repetitions;
        params.num_query_repetitions_ = options.num_query_repetitions;
        params.field_type_ = get_field_type<FieldT>(FieldT::zero());
        params.num_oracles_ = options.num_oracles;

        const FRI_snark_proof<FieldT, hash_type> proof = FRI_snark_prover<FieldT, hash_type>(params);
        printf("\n");
        print_indent(); printf("* Argument size in bytes (IOP): %zu\n", proof.IOP_size_in_bytes());
        print_indent(); printf("* Argument size in bytes (BCS): %zu\n", proof.BCS_size_in_bytes());
        print_indent(); printf("* Argument size in bytes (total): %zu\n", proof.size_in_bytes());

        printf("\nIf we were to remove the pruning of BCS merkle tree paths feature,\n"
               "the argument would have the following sizes:\n");
        print_indent(); printf("* Argument size in bytes (BCS, no pruning): %zu\n", proof.BCS_size_in_bytes_without_pruning());
        print_indent(); printf("* Argument size in bytes (total, no pruning): %zu\n", proof.size_in_bytes_without_pruning());
        printf("\n");

        const bool bit = FRI_snark_verifier<FieldT, hash_type>(proof, params);

        libiop::print_indent(); printf("* Verifier satisfied: %s\n", bit ? "true" : "false");
    }
}

int main(int argc, const char * argv[])
{
    options default_vals = {8, 20, 128, 181, 2, 0, 1, 1, 10, 
            (size_t) blake2b_type, 2, 0.1, true, true, true, false, false, blake2b_type};

#ifdef CPPDEBUG
    /* set reasonable defaults */
    if (argc > 1)
    {
        printf("There is no argument parsing in CPPDEBUG mode.");
        exit(1);
    }
    libiop::UNUSED(argv);

#else
    if (!process_prover_command_line(argc, argv, default_vals))
    {
        return 1;
    }
#endif

    start_profiling();

    printf("Selected parameters:\n");
    printf("* log_n_min = %zu\n", default_vals.log_n_min);
    printf("* log_n_max = %zu\n", default_vals.log_n_max);
    printf("* security_level = %zu\n", default_vals.security_level);

    if (default_vals.is_multiplicative) {
        switch (default_vals.field_size) {
            case 181:
                edwards_pp::init_public_params();
                instrument_FRI<edwards_Fr, binary_hash_digest>(default_vals);
                break;
            case 256:
                libff::alt_bn128_pp::init_public_params();
                instrument_FRI<alt_bn128_Fr, binary_hash_digest>(default_vals);
                break;
            default:
                throw std::invalid_argument("Field size not supported.");
        }
    } else {
        switch (default_vals.field_size)
        {
            case 64:
                instrument_FRI<gf64, binary_hash_digest>(default_vals);
                break;
            case 128:
                instrument_FRI<gf128, binary_hash_digest>(default_vals);
                break;
            case 192:
                instrument_FRI<gf192, binary_hash_digest>(default_vals);
                break;
            case 256:
                instrument_FRI<gf256, binary_hash_digest>(default_vals);
                break;
            default:
                throw std::invalid_argument("Field size not supported.");
        }
    }
}
