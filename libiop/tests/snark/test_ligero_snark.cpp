#include <cstdint>

#include <gtest/gtest.h>

#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>

#include <libff/algebra/fields/binary/gf64.hpp>
#include "libiop/relations/examples/r1cs_examples.hpp"
#include "libiop/snark/ligero_snark.hpp"
#include "libiop/bcs/common_bcs_parameters.hpp"

namespace libiop {

TEST(InterleavedR1CSSnarkTest, SimpleTest) {
    /* Set up R1CS */
    typedef libff::gf64 FieldT;

    std::size_t num_constraints = 16;
    std::size_t constraint_dim = 4;
    std::size_t num_inputs = 8;
    std::size_t num_variables = 15;
    r1cs_example<FieldT> ex = generate_r1cs_example<FieldT>(num_constraints, num_inputs, num_variables);

    r1cs_constraint_system<FieldT> constraints = ex.constraint_system_;
    r1cs_primary_input<FieldT> primary_input = ex.primary_input_;
    r1cs_auxiliary_input<FieldT> auxiliary_input = ex.auxiliary_input_;

    EXPECT_TRUE(constraints.is_satisfied(primary_input, auxiliary_input));

    /* Actual SNARK test */
    for (std::size_t i = 0; i < 2; i++)
    {
        ligero_snark_parameters<FieldT, binary_hash_digest> parameters;
        parameters.security_level_ = 128;
        parameters.height_width_ratio_ = 0.001;
        parameters.RS_extra_dimensions_ = 2;
        parameters.make_zk_ = (i == 1);
        parameters.domain_type_ = affine_subspace_type;
        parameters.LDT_reducer_soundness_type_ = LDT_reducer_soundness_type::proven;
        parameters.bcs_params_ = default_bcs_params<FieldT, binary_hash_digest>(
            blake2b_type, parameters.security_level_, constraint_dim);

        const ligero_snark_argument<FieldT, binary_hash_digest> argument =
            ligero_snark_prover<FieldT>(constraints, primary_input, auxiliary_input, parameters);

        const bool bit = ligero_snark_verifier<FieldT, binary_hash_digest>(constraints, primary_input, argument, parameters);

        EXPECT_TRUE(bit);
    }
}

TEST(InterleavedR1CSSnarkMultiplicativeTest, SimpleTest) {
    /* Set up R1CS */
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;

    std::size_t num_constraints = 16;
    std::size_t constraint_dim = 4;
    std::size_t num_inputs = 8;
    std::size_t num_variables = 15;
    r1cs_example<FieldT> ex = generate_r1cs_example<FieldT>(num_constraints, num_inputs, num_variables);

    r1cs_constraint_system<FieldT> constraints = ex.constraint_system_;
    r1cs_primary_input<FieldT> primary_input = ex.primary_input_;
    r1cs_auxiliary_input<FieldT> auxiliary_input = ex.auxiliary_input_;

    EXPECT_TRUE(constraints.is_satisfied(primary_input, auxiliary_input));

    /* Actual SNARK test */
    ligero_snark_parameters<FieldT, binary_hash_digest> parameters;
    parameters.security_level_ = 128;
    parameters.height_width_ratio_ = 0.001;
    parameters.RS_extra_dimensions_ = 2;
    parameters.make_zk_ = true;
    parameters.domain_type_ = multiplicative_coset_type;
    parameters.LDT_reducer_soundness_type_ = LDT_reducer_soundness_type::proven;
    parameters.bcs_params_ = default_bcs_params<FieldT, binary_hash_digest>(
        blake2b_type, parameters.security_level_, constraint_dim);

    const ligero_snark_argument<FieldT, binary_hash_digest> argument =
        ligero_snark_prover<FieldT, binary_hash_digest>(constraints, primary_input, auxiliary_input, parameters);

    const bool bit = ligero_snark_verifier<FieldT, binary_hash_digest>(constraints, primary_input, argument, parameters);

    EXPECT_TRUE(bit);
}

}
