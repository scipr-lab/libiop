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
    std::vector<hash_digest_type> inner_nodes_;

    std::size_t num_leaves_;
    std::shared_ptr<leafhash<FieldT, hash_digest_type>> leaf_hasher_;
    two_to_one_hash_function<hash_digest_type> node_hasher_;
    std::size_t digest_len_bytes_;
    bool make_zk_;
    std::size_t num_zk_bytes_;
    cap_hash_function<hash_digest_type> cap_hasher_;
    std::size_t cap_size_;

    /* Each element will be hashed (individually) to produce a random hash digest. */
    std::vector<zk_salt_type> zk_leaf_randomness_elements_;
    void sample_leaf_randomness();
    void compute_inner_nodes();
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

    merkle_tree_set_membership_proof<hash_digest_type> get_set_membership_proof(
        const std::vector<std::size_t> &positions) const;
    bool validate_set_membership_proof(
        const hash_digest_type &root,
        const std::vector<std::size_t> &positions,
        const std::vector<std::vector<FieldT>> &leaf_contents,
        const merkle_tree_set_membership_proof<hash_digest_type> &proof);

    /* Returns number of two to one hashes */
    size_t count_hashes_to_verify_set_membership_proof(
        const std::vector<std::size_t> &positions) const;

    std::size_t num_leaves() const;
    std::size_t depth() const;
    bool zk() const;
    std::size_t num_total_bytes() const;
};

} // namespace libiop

#include "libiop/bcs/merkle_tree.tcc"

#endif // LIBIOP_SNARK_COMMON_MERKLE_TREE_HPP_
