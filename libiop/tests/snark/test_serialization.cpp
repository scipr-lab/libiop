#include <cstdint>
#include <type_traits>
#include <sstream>

#include <gtest/gtest.h>

#include <libff/algebra/fields/binary/gf64.hpp>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>

#include "libiop/algebra/polynomials/polynomial.hpp"
#include "libiop/iop/iop.hpp"
#include "libiop/tests/bcs/dummy_bcs_protocol.hpp"
#include "libiop/bcs/bcs_prover.hpp"
#include "libiop/bcs/bcs_indexer.hpp"
#include "libiop/bcs/bcs_verifier.hpp"
#include "libiop/bcs/hashing/blake2b.hpp"
#include "libiop/bcs/hashing/dummy_algebraic_hash.hpp"
#include "libiop/snark/aurora_snark.hpp"
#include "libiop/relations/examples/r1cs_examples.hpp"


namespace libiop {

TEST(VectorSerialization, BCSTest) {
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;
    std::vector<FieldT> start = {FieldT::zero(), FieldT::one()};
    std::ostringstream s1;
    serialize_Field_Elem_vec<FieldT>(s1, start);
    std::cout << s1.str() << std::endl;
    std::istringstream s2(s1.str());
    std::vector<FieldT> end;
    deserialize_Field_Elem_vec<FieldT>(s2, end);
    assert(start[0] == end[0]);
    assert(start[1] == end[1]);    
}

TEST(TranscriptSerializationOnSnark, SimpleTest) {
    /* Set up R1CS */
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;
    typedef FieldT hash_type;

    const size_t num_constraints = 1 << 5;
    const size_t num_inputs = (1 << 2) - 1;
    const size_t num_variables = (1 << 4) - 1;
    const size_t security_parameter = 128;
    const size_t RS_extra_dimensions = 2;
    const size_t FRI_localization_parameter = 3;
    const LDT_reducer_soundness_type ldt_reducer_soundness_type = LDT_reducer_soundness_type::optimistic_heuristic;
    const FRI_soundness_type fri_soundness_type = FRI_soundness_type::heuristic;
    const field_subset_type domain_type = multiplicative_coset_type;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);
    EXPECT_TRUE(r1cs_params.constraint_system_.is_satisfied(
        r1cs_params.primary_input_, r1cs_params.auxiliary_input_));

    /* Actual SNARK test */
    for (std::size_t i = 0; i < 2; i++) {
        const bool make_zk = false;
        aurora_snark_parameters<FieldT, hash_type> params(
            security_parameter,
            ldt_reducer_soundness_type,
            fri_soundness_type,
            high_alpha_poseidon_type,
            FRI_localization_parameter,
            RS_extra_dimensions,
            make_zk,
            domain_type,
            num_constraints,
            num_variables);
        const aurora_snark_argument<FieldT, hash_type> argument = aurora_snark_prover<FieldT>(
            r1cs_params.constraint_system_,
            r1cs_params.primary_input_,
            r1cs_params.auxiliary_input_,
            params);
        
        // Test serialization
        std::ostringstream s1;
        argument.serialize(s1);
    std::cout << s1.str() << std::endl;
        std::istringstream s2(s1.str());
        aurora_snark_argument<FieldT, hash_type> deserialized_argument;
        deserialized_argument.deserialize(s2);

        const bool bit = aurora_snark_verifier<FieldT>(
            r1cs_params.constraint_system_,
            r1cs_params.primary_input_,
            deserialized_argument,
            params);

        EXPECT_TRUE(bit) << "failed on make_zk = " << i << " test";
    }
}

}
