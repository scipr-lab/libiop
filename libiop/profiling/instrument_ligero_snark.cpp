#include <algorithm>
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

#include <libff/algebra/fields/binary/gf64.hpp>
#include <libff/algebra/fields/binary/gf256.hpp>
#include <libff/algebra/fields/binary/gf192.hpp>
#include <libff/algebra/fields/binary/gf256.hpp>

#include "libiop/snark/ligero_snark.hpp"
#include "libiop/bcs/bcs_common.hpp"
#include "libiop/bcs/common_bcs_parameters.hpp"
#include "libiop/relations/examples/r1cs_examples.hpp"

#ifndef CPPDEBUG
bool process_prover_command_line(const int argc, const char** argv,
                                 options &options, 
                                 float &height_width_ratio,
                                 std::size_t &RS_extra_dimensions)
{
    namespace po = boost::program_options;

    try
    {
        po::options_description desc = gen_options(options);
        desc.add_options()
             ("height_width_ratio", po::value<float>(&height_width_ratio)->default_value(0.1))
             ("RS_extra_dimensions", po::value<std::size_t>(&RS_extra_dimensions)->default_value(2));

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
void instrument_ligero_snark(options &options,
                             LDT_reducer_soundness_type ldt_reducer_soundness_type,
                             const field_subset_type domain_type,
                             float height_width_ratio, 
                             std::size_t RS_extra_dimensions)
{
    ligero_snark_parameters<FieldT, hash_type> parameters;
    parameters.security_level_ = options.security_level;
    parameters.LDT_reducer_soundness_type_ = ldt_reducer_soundness_type;
    parameters.height_width_ratio_ = height_width_ratio;
    parameters.RS_extra_dimensions_ = RS_extra_dimensions;
    parameters.make_zk_ = options.make_zk;
    parameters.domain_type_ = domain_type;
    parameters.bcs_params_ = default_bcs_params<FieldT, hash_type>(options.hash_enum, options.security_level, options.log_n_min);
    parameters.describe();

    for (std::size_t log_n = options.log_n_min; log_n <= options.log_n_max; ++log_n)
    {
        libff::print_separator();
        const std::size_t n = 1ul << log_n;
        /* k+1 needs to be a power of 2 (proof system artifact) and k <= n+2 (example generation artifact) so we just fix it to 15 here */
        const std::size_t k = 15;
        const std::size_t m = n - 1;
        r1cs_example<FieldT> example = generate_r1cs_example<FieldT>(n, k, m);
        parameters.bcs_params_ = default_bcs_params<FieldT, hash_type>(options.hash_enum, options.security_level, log_n);

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
        const ligero_snark_argument<FieldT, hash_type> proof = ligero_snark_prover<FieldT, hash_type>(
            example.constraint_system_,
            example.primary_input_,
            example.auxiliary_input_,
            parameters);

        printf("\n");

        libff::print_indent(); printf("* Argument size in bytes (IOP): %zu\n", proof.IOP_size_in_bytes());
        libff::print_indent(); printf("* Argument size in bytes (BCS): %zu\n", proof.BCS_size_in_bytes());
        libff::print_indent(); printf("* Argument size in bytes (total): %zu\n", proof.size_in_bytes());

        printf("\nIf we were to remove pruning of authentication paths in BCS,\n"
               "the argument would have the following sizes:\n");
        libff::print_indent(); printf("* Argument size in bytes (BCS, no pruning): %zu\n", proof.BCS_size_in_bytes_without_pruning());
        libff::print_indent(); printf("* Argument size in bytes (total, no pruning): %zu\n", proof.size_in_bytes_without_pruning());

        printf("\n");

        const bool bit = ligero_snark_verifier<FieldT, hash_type>(
            example.constraint_system_,
            example.primary_input_,
            proof,
            parameters);

        printf("\n\n");

        libff::print_indent(); printf("* Verifier satisfied: %s\n", bit ? "true" : "false");
    }
}

int main(int argc, const char * argv[])
{

    options default_vals;

    float height_width_ratio = 0.1;
    std::size_t RS_extra_dimensions = 2;

#ifdef CPPDEBUG
    /* set reasonable defaults */

#else
    if (!process_prover_command_line(argc, argv, default_vals, height_width_ratio, RS_extra_dimensions))
    {
        return 1;
    }
#endif

    /** TODO: eventually get a string from program options, and then have a from string method in LDT reducer */
    LDT_reducer_soundness_type ldt_reducer_soundness_type = LDT_reducer_soundness_type::proven;
    if (default_vals.heuristic_ldt_reducer_soundness)
    {
        ldt_reducer_soundness_type = LDT_reducer_soundness_type::optimistic_heuristic;
    }
    libff::start_profiling();

    printf("Selected parameters:\n");
    printf("- log_n_min = %zu\n", default_vals.log_n_min);
    printf("- log_n_max = %zu\n", default_vals.log_n_max);
    printf("- height_width_ratio = %f\n", height_width_ratio);
    printf("- RS_extra_dimensions = %zu\n", RS_extra_dimensions);
    printf("- security_level = %zu\n", default_vals.security_level);
    printf("- LDT_reducer_soundness_type = %s\n", LDT_reducer_soundness_type_to_string(ldt_reducer_soundness_type));
    printf("- field_size = %zu\n", default_vals.field_size);
    printf("- make_zk = %d\n", default_vals.make_zk);
    printf("- hash_enum = %s\n", bcs_hash_type_names[default_vals.hash_enum]);

    if (default_vals.is_multiplicative)
    {
        switch (default_vals.field_size) {
            case 181:
                libff::edwards_pp::init_public_params();
                instrument_ligero_snark<libff::edwards_Fr, binary_hash_digest>(
                                        default_vals, ldt_reducer_soundness_type, multiplicative_coset_type,
                                        height_width_ratio, RS_extra_dimensions);
                break;
            case 256:
                libff::alt_bn128_pp::init_public_params();
                if (default_vals.hash_enum == blake2b_type)
                {
                    instrument_ligero_snark<libff::alt_bn128_Fr, binary_hash_digest>(
                                            default_vals, ldt_reducer_soundness_type, multiplicative_coset_type,
                                            height_width_ratio, RS_extra_dimensions);
                } 
                else
                {
                    instrument_ligero_snark<libff::alt_bn128_Fr, libff::alt_bn128_Fr>(
                                            default_vals, ldt_reducer_soundness_type, multiplicative_coset_type,
                                            height_width_ratio, RS_extra_dimensions);
                }
                break;
            default:
                throw std::invalid_argument("Field size not supported.");
        }
    }
    else
    {
        switch (default_vals.field_size)
        {
            case 64:
                instrument_ligero_snark<libff::gf64, binary_hash_digest>(
                                        default_vals, ldt_reducer_soundness_type, affine_subspace_type,
                                        height_width_ratio, RS_extra_dimensions);
                break;
            case 128:
                instrument_ligero_snark<libff::gf128, binary_hash_digest>(
                                        default_vals, ldt_reducer_soundness_type, affine_subspace_type,
                                        height_width_ratio, RS_extra_dimensions);
                break;
            case 192:
                instrument_ligero_snark<libff::gf192, binary_hash_digest>(
                                        default_vals, ldt_reducer_soundness_type, affine_subspace_type,
                                        height_width_ratio, RS_extra_dimensions);
                break;
            case 256:
                instrument_ligero_snark<libff::gf256, binary_hash_digest>(
                                        default_vals, ldt_reducer_soundness_type, affine_subspace_type,
                                        height_width_ratio, RS_extra_dimensions);
                break;
            default:
                throw std::invalid_argument("Field size not supported.");
        }
    }
}
