/**@file
 *****************************************************************************
 Merkle tree interfaces.

 Includes support for zero knowledge merkle trees, and set membership-proofs.
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_SNARK_COMMON_MERKLE_TREE_HPP_
#define LIBIOP_SNARK_COMMON_MERKLE_TREE_HPP_

#include <cstddef>
#include <numeric>
#include <vector>
#include <bits/stdc++.h>

#include "libiop/algebra/field_subset/field_subset.hpp"
#include "libiop/bcs/hashing/hashing.hpp"

namespace libiop {

/* Authentication paths for a set of positions */
template<typename hash_digest_type>
struct merkle_tree_set_membership_proof {
    std::vector<hash_digest_type> auxiliary_hashes;
    std::vector<zk_salt_type> randomness_hashes;

    /* TODO: Write a test for this */
    std::size_t size_in_bytes() const
    {
        return std::accumulate(this->auxiliary_hashes.begin(),
                               this->auxiliary_hashes.end(),
                               0,
                               [] (const std::size_t av, const hash_digest_type &h) { 
                                   return av + get_hash_size<hash_digest_type>(h); }) +
               std::accumulate(this->randomness_hashes.begin(),
                               this->randomness_hashes.end(),
                               0,
                               [] (const std::size_t av, const zk_salt_type &h) {
                                    return av + get_hash_size<zk_salt_type>(h); });
    }
};

template<typename FieldT, typename hash_digest_type>
class merkle_tree {
protected:
    bool constructed_;
    /* inner_nodes_ is a vector of the (num_leaves - 1) nodes in the tree, with the root at
       index 0, left child at 1, right child at 2, etc. If cap_size_ is greater than 2, the
       first (log_2(cap_size_) - 1) layers under the root are empty to make the math easier. */
    std::vector<hash_digest_type> inner_nodes_;

    std::size_t num_leaves_;
    std::shared_ptr<leafhash<FieldT, hash_digest_type>> leaf_hasher_;
    two_to_one_hash_function<hash_digest_type> node_hasher_;
    std::size_t digest_len_bytes_;
    bool make_zk_;
    std::size_t num_zk_bytes_;
    /* The top log_2(cap_size_) layers are hashed with a single computation to improve efficiency.
       The root along with its cap_size_ direct children are referred to as the "cap," and the
       operation that transforms these children to the root is the cap hash.
       See https://github.com/scipr-lab/libiop/issues/41. */
    cap_hash_function<hash_digest_type> cap_hasher_;
    /* cap_size_ is the number of direct children the root has. It must be a power of 2 and at
       least 2. For example if cap_size == 4, the root has 4 children, and in inner_nodes_ the
       indices 1 and 2 are unused. */
    std::size_t cap_size_;

    /* Each element will be hashed (individually) to produce a random hash digest. */
    std::vector<zk_salt_type> zk_leaf_randomness_elements_;
    void sample_leaf_randomness();
    void compute_inner_nodes();

    /* Helper functions for dealing with the tree strucutre. Correctness not guaranteed
       when out of bounds. */
    std::size_t parent_of(const std::size_t node_index) const;
    std::size_t left_child_of(const std::size_t node_index) const;
    std::size_t right_child_of(const std::size_t node_index) const;
    bool is_in_cap(const std::size_t node_index) const;
    std::size_t cap_children_start() const; // Inclusive.
    std::size_t cap_children_end() const; // Exclusive.
public:
    /* Create a merkle tree with the given configuration.
    If make_zk is true, 2 * security parameter random bytes will be appended to each leaf
    before hashing, to prevent a low entropy leaf value from being inferred from its hash.
    cap_size is the number of children of the root and must be a power of 2. */
    merkle_tree(const std::size_t num_leaves,
                const std::shared_ptr<leafhash<FieldT, hash_digest_type>> &leaf_hasher,
                const two_to_one_hash_function<hash_digest_type> &node_hasher,
                const cap_hash_function<hash_digest_type> &cap_hasher,
                const std::size_t digest_len_bytes,
                const bool make_zk,
                const std::size_t security_parameter,
                const std::size_t cap_size=2);

    /** This treats each leaf as a column.
     * e.g. The ith leaf is the vector formed by leaf_contents[j][i] for all j */
    void construct(const std::vector<std::shared_ptr<std::vector<FieldT>>> &leaf_contents);
    // TODO: Remove this overload in favor of only using the former
    void construct(const std::vector<std::vector<FieldT> > &leaf_contents);
    /** Leaf contents is a table with `r` rows
     *  (`r` typically being the number of oracles)
     *  and (MT_num_leaves * coset_serialization_size) columns.
     *  Each MT leaf is the serialization of a table with `r` rows,
     *  and coset_serialization_size columns.
     *
     *  This is done here rather than the BCS layer to avoid needing to copy the data,
     *  as this will take a significant amount of memory.
     */
    void construct_with_leaves_serialized_by_cosets(
        const std::vector<std::shared_ptr<std::vector<FieldT>>> &leaf_contents,
        size_t coset_serialization_size);

    /** Takes in a set of query positions to input oracles to a domain of size:
     *  `num_leaves * coset_serialization_size`,
     *  and the associated evaluations for each query position.
     *
     *  This function then serializes these evaluations into leaf entries.
     *  The rows of a leaf entry are the same as in the eva
    */
    std::vector<std::vector<FieldT>> serialize_leaf_values_by_coset(
        const std::vector<size_t> &query_positions,
        const std::vector<std::vector<FieldT> > &query_responses,
        const size_t coset_serialization_size) const;

    hash_digest_type get_root() const;

    /* These two functions do not currently work if the given positions aren't sorted or
       have duplicates, AND the tree is set to be zero knowledge. */
    merkle_tree_set_membership_proof<hash_digest_type> get_set_membership_proof(
        const std::vector<std::size_t> &positions) const;
    bool validate_set_membership_proof(
        const hash_digest_type &root,
        const std::vector<std::size_t> &positions,
        const std::vector<std::vector<FieldT>> &leaf_contents,
        const merkle_tree_set_membership_proof<hash_digest_type> &proof);

    /** Returns a number that is proportional to the hashing runtime of verifying a set membership
     *  proof. Each two-to-one hash is counted as 2 units, and each input of the cap hash is 1 unit.
     *  Leaf hashes are not counted. */
    size_t count_internal_hash_complexity_to_verify_set_membership(
        const std::vector<std::size_t> &positions) const;

    std::size_t num_leaves() const;
    std::size_t depth() const;
    bool zk() const;
    std::size_t num_total_bytes() const;
};

} // namespace libiop

#include "libiop/bcs/merkle_tree.tcc"

#endif // LIBIOP_SNARK_COMMON_MERKLE_TREE_HPP_
