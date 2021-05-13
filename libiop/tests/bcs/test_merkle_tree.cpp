#include <cstdint>
#include <gtest/gtest.h>
#include <vector>
#include <type_traits>

#include <libff/algebra/fields/binary/gf64.hpp>
#include "libiop/algebra/utils.hpp"
#include "libiop/bcs/hashing/blake2b.hpp"
#include "libiop/bcs/hashing/dummy_algebraic_hash.hpp"
#include "libiop/bcs/hashing/hashing.hpp"
#include "libiop/bcs/merkle_tree.hpp"

namespace libiop {

template< bool B, class T = void >
using enable_if_t = typename libff::enable_if<B,T>::type;
// Binary hash type
template<typename FieldT, typename hash_type,
    enable_if_t<!std::is_same<FieldT, hash_type>::value, int> = 42>
merkle_tree<FieldT, hash_type> new_MT(
    const std::size_t size, const std::size_t digest_len_bytes, const bool make_zk,
    const std::size_t security_parameter)
{
    return merkle_tree<FieldT, binary_hash_digest>(
        size,
        std::make_shared<blake2b_leafhash<FieldT>>(security_parameter),
        blake2b_two_to_one_hash,
        digest_len_bytes,
        make_zk,
        security_parameter);
}
template<typename FieldT, typename hash_type,
    enable_if_t<std::is_same<FieldT, hash_type>::value, int> = 42>
merkle_tree<FieldT, hash_type> new_MT(
    const std::size_t size, const std::size_t digest_len_bytes, const bool make_zk,
    const std::size_t security_parameter)
{
    return merkle_tree<FieldT, FieldT>(
        size,
        std::make_shared<dummy_algebraic_leafhash<FieldT>>(),
        dummy_algebraic_two_to_one_hash<FieldT>,
        digest_len_bytes,
        make_zk,
        security_parameter);
}

template<typename FieldT, typename hash_type>
void run_simple_MT_test(const std::size_t size, const std::size_t digest_len_bytes, const bool make_zk,
                        const std::size_t security_parameter) {
    merkle_tree<FieldT, hash_type> tree =
        new_MT<FieldT, hash_type>(size, digest_len_bytes, make_zk, security_parameter);
    const std::vector<FieldT> vec1 = random_vector<FieldT>(size);
    const std::vector<FieldT> vec2 = random_vector<FieldT>(size);

    tree.construct({ vec1, vec2 });

    const hash_type root = tree.get_root();

    for (std::size_t i = 0; i < size; ++i)
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

    const std::size_t size = 16;
    const std::size_t digest_len_bytes = 256/8;
    const std::size_t security_parameter = 128;
    run_simple_MT_test<FieldT, binary_hash_digest>(size, digest_len_bytes, false, security_parameter);
    run_simple_MT_test<FieldT, FieldT>(size, digest_len_bytes, false, security_parameter);
}

TEST(MerkleTreeZKTest, SimpleTest) {
    typedef libff::gf64 FieldT;

    const std::size_t size_small = 16;
    const std::size_t size_large = 1ull << 18; /* The goal is to test batch randomness logic */
    const std::size_t digest_len_bytes = 256/8;
    const std::size_t security_parameter = 128;
    run_simple_MT_test<FieldT, binary_hash_digest>(size_small, digest_len_bytes, true, security_parameter);
    run_simple_MT_test<FieldT, FieldT>(size_large, digest_len_bytes, true, security_parameter);
}

void run_multi_test(const bool make_zk) {
    typedef libff::gf64 FieldT;

    const std::size_t size = 8;
    const std::size_t security_parameter = 128;
    const std::size_t digest_len_bytes = 256/8;
    const bool algebraic_hash = false;

    merkle_tree<FieldT, binary_hash_digest> tree = new_MT<FieldT, binary_hash_digest>(
        size,
        digest_len_bytes,
        make_zk,
        security_parameter);

    const std::vector<FieldT> vec1 = random_vector<FieldT>(size);
    const std::vector<FieldT> vec2 = random_vector<FieldT>(size);

    tree.construct({ vec1, vec2 });

    const binary_hash_digest root = tree.get_root();

    std::vector<std::vector<FieldT>> leafs;
    for (std::size_t i = 0; i < size; ++i)
    {
        std::vector<FieldT> leaf({ vec1[i], vec2[i] });
        leafs.emplace_back(leaf);
    }

    for (std::size_t subset = 0; subset < (1ull<<size); ++subset)
    {
        std::vector<std::size_t> subset_elements;
        std::vector<std::vector<FieldT>> subset_leafs;
        for (std::size_t k = 0; k < size; ++k)
        {
            if (subset & (1ull<<k))
            {
                subset_elements.emplace_back(k);
                subset_leafs.emplace_back(leafs[k]);
            }
        }

        const merkle_tree_set_membership_proof<binary_hash_digest> mp = tree.get_set_membership_proof(subset_elements);

        const bool is_valid = tree.validate_set_membership_proof(root,
                                                                 subset_elements,
                                                                 subset_leafs,
                                                                 mp);
        EXPECT_TRUE(is_valid);
    }
}

TEST(MerkleTreeTest, MultiTest) {
    const bool make_zk = false;
    run_multi_test(make_zk);
}

TEST(MerkleTreeZKTest, MultiTest) {
    const bool make_zk = true;
    run_multi_test(make_zk);
}

TEST(MerkleTreeTwoToOneHashTest, SimpleTest)
{
    typedef libff::gf64 FieldT;
    bool make_zk = false;
    const size_t num_leaves = 8;
    const size_t security_parameter = 128;
    const size_t hash_length = 32;
    const bool algebraic_hash = false;

    merkle_tree<FieldT, binary_hash_digest> tree = new_MT<FieldT, binary_hash_digest>(
        num_leaves,
        hash_length,
        make_zk,
        security_parameter);

    std::vector<size_t> positions = {1, 3, 6, 7};
    size_t expected_num_hashes = 6;
    size_t actual_num_hashes = tree.count_hashes_to_verify_set_membership_proof(positions);
    ASSERT_EQ(expected_num_hashes, actual_num_hashes);
}

}
