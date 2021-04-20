#include <cstdint>

#include <gtest/gtest.h>

#include "libiop/bcs/hashing/hashing.hpp"
#include "libiop/bcs/hashing/hash_enum.hpp"
#include "libiop/bcs/pow.hpp"

#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>

namespace libiop {

TEST(BinaryPoWTest, SimpleTest) {
    /* Arbitrary field */
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;
    typedef binary_hash_digest hash_type;
    const size_t security_parameter = 128;
    const size_t digets_len_bytes = 2 * security_parameter/8;

    const size_t log_work = 20;
    const size_t cost_per_hash = 1;
    pow_parameters params = pow_parameters(log_work, cost_per_hash);
    pow<FieldT, hash_type> prover = pow<FieldT, hash_type>(params, digets_len_bytes);

    // 32 char challenge
    const std::string challenge = "abcdefghijklmnopqrstuvwxyzabcdef";
    two_to_one_hash_function<hash_type> compressive_hash = 
        get_two_to_one_hash<hash_type, FieldT>(blake2b_type, security_parameter);

    const hash_type proof = prover.solve_pow(compressive_hash, challenge);
    EXPECT_TRUE(prover.verify_pow(compressive_hash, challenge, proof));
}

TEST(AlgeraicPoWTest, SimpleTest) {
    /* Set up field / pow params */
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;
    typedef FieldT hash_type;
    const size_t digets_len_bytes = 32;

    const size_t log_work = 20;
    const size_t cost_per_hash = 200;
    pow_parameters params = pow_parameters(log_work, cost_per_hash);
    pow<FieldT, hash_type> prover = pow<FieldT, hash_type>(params, digets_len_bytes);

    const FieldT challenge = FieldT::random_element();
    two_to_one_hash_function<hash_type> compressive_hash = 
        get_two_to_one_hash<FieldT, hash_type>(high_alpha_poseidon_type, 128);

    const hash_type proof = prover.solve_pow(compressive_hash, challenge);
    EXPECT_TRUE(prover.verify_pow(compressive_hash, challenge, proof));
}

}
