#include <cstdint>
#include <gtest/gtest.h>
#include <vector>
#include <type_traits>
#include <random>

#include <libff/algebra/fields/binary/gf64.hpp>
#include "libiop/algebra/utils.hpp"
#include "libiop/bcs/hashing/blake2b.hpp"
#include "libiop/bcs/hashing/dummy_algebraic_hash.hpp"
#include "libiop/bcs/hashing/hashing.hpp"
#include "libiop/bcs/merkle_tree.hpp"

namespace libiop {

using std::size_t;

template< bool B, class T = void >
using enable_if_t = typename libff::enable_if<B,T>::type;

// Binary hash type
template<typename FieldT, typename hash_type,
    enable_if_t<!std::is_same<FieldT, hash_type>::value, int> = 42>
merkle_tree<FieldT, hash_type> new_MT(
    const size_t size, const size_t digest_len_bytes, const bool make_zk,
    const size_t security_parameter, const size_t cap_size=2)
{
    return merkle_tree<FieldT, binary_hash_digest>(
        size,
        std::make_shared<blake2b_leafhash<FieldT>>(security_parameter),
        blake2b_two_to_one_hash,
        blake2b_many_to_one_hash,
        digest_len_bytes,
        make_zk,
        security_parameter,
        cap_size);
}

template<typename FieldT, typename hash_type,
    enable_if_t<std::is_same<FieldT, hash_type>::value, int> = 42>
merkle_tree<FieldT, hash_type> new_MT(
    const size_t size, const size_t digest_len_bytes, const bool make_zk,
    const size_t security_parameter, const size_t cap_size=2)
{
    return merkle_tree<FieldT, FieldT>(
        size,
        std::make_shared<dummy_algebraic_leafhash<FieldT>>(),
        dummy_algebraic_two_to_one_hash<FieldT>,
        dummy_algebraic_cap_hash<FieldT>,
        digest_len_bytes,
        make_zk,
        security_parameter,
        cap_size);
}

/** Constructs a merkle tree with leaf size 2. Generates and verifies membership proofs for
 *  each leaf individually, and makes sure reversing the contents of each leaf causes the
 *  verification to fail (unless the leaf contents are symmetric). */
template<typename FieldT, typename hash_type>
void run_simple_MT_test(const size_t size, const size_t digest_len_bytes, const bool make_zk,
                        const size_t security_parameter, const size_t cap_size) {
    merkle_tree<FieldT, hash_type> tree =
        new_MT<FieldT, hash_type>(size, digest_len_bytes, make_zk, security_parameter, cap_size);
    const std::vector<FieldT> vec1 = random_vector<FieldT>(size);
    const std::vector<FieldT> vec2 = random_vector<FieldT>(size);

    tree.construct({ vec1, vec2 });

    const hash_type root = tree.get_root();

    for (size_t i = 0; i < size; ++i)
    {
        /* membership proof for the set {i} */
        const std::vector<size_t> set = {i};

        const merkle_tree_set_membership_proof<hash_type> ap = tree.get_set_membership_proof(set);
        const std::vector<std::vector<FieldT>> contents =
            { {vec1[i] , vec2[i]} };

        const bool is_valid = tree.validate_set_membership_proof(
            root,
            set,
            contents,
            ap);
        EXPECT_TRUE(is_valid);

        const std::vector<std::vector<FieldT>> reverse_contents =
            { {vec2[i] , vec1[i]} };

        const bool reverse_is_valid = tree.validate_set_membership_proof(
            root,
            set,
            reverse_contents,
            ap);
        if (vec1[i] == vec2[i])
        {
            EXPECT_EQ(contents, reverse_contents);
            EXPECT_TRUE(reverse_is_valid);
        }
        else
        {
            EXPECT_NE(contents, reverse_contents);
            EXPECT_FALSE(reverse_is_valid);
        }
    }
}

TEST(MerkleTreeTest, SimpleTest) {
    typedef libff::gf64 FieldT;

    const size_t size = 16;
    const std::vector<size_t> cap_sizes = {2, 4, 8, 16}; // Test all possible cap sizes.
    const size_t digest_len_bytes = 256/8;
    const size_t security_parameter = 128;

    for (size_t cap_size : cap_sizes)
    {
        run_simple_MT_test<FieldT, binary_hash_digest>(size, digest_len_bytes, false,
                                                       security_parameter, cap_size);
        run_simple_MT_test<FieldT, FieldT>(size, digest_len_bytes, false,
                                           security_parameter, cap_size);
    }
}

TEST(MerkleTreeZKTest, SimpleTest) {
    typedef libff::gf64 FieldT;

    const size_t size_small = 16;
    const size_t size_large = 1ull << 18; /* The goal is to test batch randomness logic */
    const size_t cap_size = 4;
    const size_t digest_len_bytes = 256/8;
    const size_t security_parameter = 128;
    run_simple_MT_test<FieldT, binary_hash_digest>(size_small, digest_len_bytes, true,
                                                   security_parameter, cap_size);
    run_simple_MT_test<FieldT, FieldT>(size_large, digest_len_bytes, true,
                                       security_parameter, cap_size);
}

/** Constructs a merkle tree with 8 leaves each of size 2, and cap size 4. Generates and verifies
 *  membership proofs for every possible subset of leaves. */
void run_fixed_multi_test(const bool make_zk) {
    typedef libff::gf64 FieldT;

    // The size is fixed because large values would quickly cause the program run out of memory.
    const size_t size = 8;
    const size_t cap_size = 4;
    const size_t security_parameter = 128;
    const size_t digest_len_bytes = 256/8;
    const bool algebraic_hash = false;

    merkle_tree<FieldT, binary_hash_digest> tree = new_MT<FieldT, binary_hash_digest>(
        size,
        digest_len_bytes,
        make_zk,
        security_parameter,
        cap_size);

    const std::vector<FieldT> vec1 = random_vector<FieldT>(size);
    const std::vector<FieldT> vec2 = random_vector<FieldT>(size);

    tree.construct({ vec1, vec2 });

    const binary_hash_digest root = tree.get_root();

    std::vector<std::vector<FieldT>> leaves;
    for (size_t i = 0; i < size; ++i)
    {
        std::vector<FieldT> leaf({ vec1[i], vec2[i] });
        leaves.emplace_back(leaf);
    }

    /* This code generates every possible subset. `subset` is a binary string that encodes for each
       element, whether it is in this subset. */
    for (size_t subset = 0; subset < (1ull<<size); ++subset)
    {
        std::vector<size_t> subset_elements;
        std::vector<std::vector<FieldT>> subset_leaves;
        for (size_t k = 0; k < size; ++k)
        {
            if (subset & (1ull<<k))
            {
                subset_elements.emplace_back(k);
                subset_leaves.emplace_back(leaves[k]);
            }
        }

        const merkle_tree_set_membership_proof<binary_hash_digest> mp =
                tree.get_set_membership_proof(subset_elements);

        const bool is_valid = tree.validate_set_membership_proof(root,
                                                                 subset_elements,
                                                                 subset_leaves,
                                                                 mp);
        EXPECT_TRUE(is_valid);
    }
}

TEST(MerkleTreeTest, FixedMultiTest) {
    const bool make_zk = false;
    run_fixed_multi_test(make_zk);
}

TEST(MerkleTreeZKTest, FixedMultiTest) {
    const bool make_zk = true;
    run_fixed_multi_test(make_zk);
}

/** Constructs a merkle tree with leaf size 2. Generates and verifies membership proofs for some
 *  randomly generated sorted subset of leaves of specified size, with no duplicates. Queries with
 *  unsorted, duplicated lists of leaves currently only work when it is not zero knowledge. */
void run_random_multi_test(const size_t size, const size_t digest_len_bytes, const bool make_zk,
                          const size_t security_parameter, const size_t cap_size,
                          const size_t subset_size) {
    typedef libff::gf64 FieldT;

    const bool algebraic_hash = false;
    const size_t num_iterations = 1; // The number of randomly generated subsets to test.

    merkle_tree<FieldT, binary_hash_digest> tree = new_MT<FieldT, binary_hash_digest>(
        size,
        digest_len_bytes,
        make_zk,
        security_parameter,
        cap_size);

    const std::vector<FieldT> vec1 = random_vector<FieldT>(size);
    const std::vector<FieldT> vec2 = random_vector<FieldT>(size);

    tree.construct({ vec1, vec2 });

    const binary_hash_digest root = tree.get_root();

    std::vector<std::vector<FieldT>> leaves;
    leaves.reserve(size);
    std::vector<size_t> shuffled_leaf_indices;
    shuffled_leaf_indices.reserve(size);
    for (size_t i = 0; i < size; ++i)
    {
        std::vector<FieldT> leaf({ vec1[i], vec2[i] });
        leaves.emplace_back(leaf);
        shuffled_leaf_indices.emplace_back(i);
    }

    for (size_t i = 0; i < num_iterations; i++)
    {
        std::vector<size_t> subset_elements;
        std::vector<std::vector<FieldT>> subset_leaves;
        /* The commented-out code generates subsets that are unsorted and may be repeats.
           They are not used because the code currently cannot handle these cases if it is
           zero knowledge. */
        // for (size_t j = 0; j < subset_size; j++)
        // {
        //     size_t k = randombytes_uniform(size);
        //     subset_elements.emplace_back(k);
        //     subset_leaves.emplace_back(leaves[k]);
        // }

        // Generate a random sorted subset of indices at the beginning of shuffled_leaf_indices.
        std::shuffle(shuffled_leaf_indices.begin(), shuffled_leaf_indices.end(),
                     std::default_random_engine(i));
        std::sort(shuffled_leaf_indices.begin(), shuffled_leaf_indices.begin() + subset_size);
        for (size_t j = 0; j < subset_size; j++)
        {
            size_t k = shuffled_leaf_indices[j];
            subset_elements.emplace_back(k);
            subset_leaves.emplace_back(leaves[k]);
        }

        const merkle_tree_set_membership_proof<binary_hash_digest> mp =
                tree.get_set_membership_proof(subset_elements);

        const bool is_valid = tree.validate_set_membership_proof(root,
                                                                 subset_elements,
                                                                 subset_leaves,
                                                                 mp);
        EXPECT_TRUE(is_valid);
    }
}

TEST(MerkleTreeTest, RandomMultiTest) {
    const size_t security_parameter = 128;
    const size_t digest_len_bytes = 256/8;
    const bool make_zk = false;
    // Test a small and a large tree.
    run_random_multi_test(16, digest_len_bytes, make_zk, security_parameter, 4, 5);
    run_random_multi_test(1ull << 16, digest_len_bytes, make_zk, security_parameter, 256, 100);
}

TEST(MerkleTreeZKTest, RandomMultiTest) {
    const size_t security_parameter = 128;
    const size_t digest_len_bytes = 256/8;
    const bool make_zk = true;
    // Test a small and a large tree.
    run_random_multi_test(16, digest_len_bytes, make_zk, security_parameter, 4, 5);
    run_random_multi_test(1ull << 16, digest_len_bytes, make_zk, security_parameter, 256, 100);
}

/** Verify that count_internal_hash_complexity_to_verify_set_membership is correct for a fixed tree
 *  size and query set, for various cap sizes. */
TEST(MerkleTreeHashCountTest, SimpleTest)
{
    typedef libff::gf64 FieldT;
    bool make_zk = false;
    const size_t num_leaves = 8;
    const size_t security_parameter = 128;
    const size_t hash_length = 32;
    const bool algebraic_hash = false;

    const std::vector<size_t> cap_sizes = {2, 4, 8};
    const std::vector<size_t> expected_num_hashes = {12, 10, 8};

    const std::vector<size_t> positions = {1, 3, 6, 7};

    for (size_t i = 0; i < cap_sizes.size(); i++)
    {
        merkle_tree<FieldT, binary_hash_digest> tree = new_MT<FieldT, binary_hash_digest>(
            num_leaves,
            hash_length,
            make_zk,
            security_parameter,
            cap_sizes[i]);

        size_t actual_num_hashes = tree.count_internal_hash_complexity_to_verify_set_membership(
            positions);
        ASSERT_EQ(expected_num_hashes[i], actual_num_hashes);
    }
}

}
