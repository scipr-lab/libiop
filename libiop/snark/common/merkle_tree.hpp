/**@file
 *****************************************************************************
 Merkle tree interfaces.

 Includes support for zero knowledge merkle trees, and multi-membership proofs.
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_SNARK_COMMON_MERKLE_TREE_HPP_
#define LIBIOP_SNARK_COMMON_MERKLE_TREE_HPP_

#include <cstddef>
#include <numeric>
#include <vector>

#include "libiop/snark/common/hashing.hpp"

namespace libiop {

/* Authentication path for a single position */
struct merkle_tree_membership_proof {
    std::vector<hash_digest> auxiliary_hashes;
    hash_digest randomness_hash;

    /* TODO: factor out to .cpp */
    std::size_t size_in_bytes() const
    {
        return std::accumulate(this->auxiliary_hashes.begin(),
                               this->auxiliary_hashes.end(),
                               0,
                               [] (const std::size_t av, const hash_digest &h) { return av + h.size(); });
    }
};

/* Authentication paths for multiple positions */
struct merkle_tree_multi_membership_proof {
    std::vector<hash_digest> auxiliary_hashes;
    std::vector<hash_digest> randomness_hashes;

    /* TODO: Write a test for this */
    std::size_t size_in_bytes() const
    {
        return std::accumulate(this->auxiliary_hashes.begin(),
                               this->auxiliary_hashes.end(),
                               0,
                               [] (const std::size_t av, const hash_digest &h) { return av + h.size(); }) +
               std::accumulate(this->randomness_hashes.begin(),
                               this->randomness_hashes.end(),
                               0,
                               [] (const std::size_t av, const hash_digest &h) { return av + h.size(); });
    }
};

template<typename FieldT>
class merkle_tree {
protected:
    bool constructed_;
    std::vector<hash_digest> inner_nodes_;

    std::size_t num_leaves_;
    field_element_hash_function<FieldT> leaf_hasher_;
    zk_element_hash_function zk_leaf_hasher_;
    two_to_one_hash_function node_hasher_;
    std::size_t digest_len_bytes_;
    bool make_zk_;
    std::size_t num_zk_bytes_;

    /* Each element will be hashed (individually) to produce a random hash digest. */
    std::vector<std::vector<uint8_t>> zk_leaf_randomness_elements_;

public:
    /* Create a merkle tree with the given configuration.
    If make_zk is true, 2 * security parameter random bytes will be appended to each leaf
    before hashing, to prevent a low entropy leaf value from being inferred
    from its hash. */
    merkle_tree(const std::size_t num_leaves,
                const field_element_hash_function<FieldT> &leaf_hasher,
                const zk_element_hash_function &zk_leaf_hasher,
                const two_to_one_hash_function &node_hasher,
                const std::size_t digest_len_bytes,
                const bool make_zk,
                const std::size_t security_parameter);

    /** This treats each leaf as a column.
     * e.g. The ith leaf is the vector formed by leaf_contents[j][i] for all j */
    void construct(const std::vector<std::vector<FieldT> > &leaf_contents);

    hash_digest get_root() const;

    merkle_tree_membership_proof get_membership_proof(const std::size_t position) const;
    bool validate_membership_proof(const hash_digest &root,
                                   const std::size_t position,
                                   const hash_digest& contents_hash,
                                   const merkle_tree_membership_proof &proof);

    merkle_tree_multi_membership_proof get_multi_membership_proof(
        const std::vector<std::size_t> &positions) const;
    bool validate_multi_membership_proof(
        const hash_digest &root,
        const std::vector<std::size_t> &positions,
        const std::vector<hash_digest> &contents_hashes,
        const merkle_tree_multi_membership_proof &proof);

    std::size_t num_leaves() const;
    std::size_t depth() const;
    std::size_t num_total_bytes() const;
};

} // namespace libiop

#include "libiop/snark/common/merkle_tree.tcc"

#endif // LIBIOP_SNARK_COMMON_MERKLE_TREE_HPP_
