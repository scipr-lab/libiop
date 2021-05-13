/**@file
 *****************************************************************************
 Blake2b implementation of relevant hash functions
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_SNARK_COMMON_HASHING_BLAKE2B_HPP_
#define LIBIOP_SNARK_COMMON_HASHING_BLAKE2B_HPP_

#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>
#include <libff/algebra/field_utils/field_utils.hpp>
#include "libiop/bcs/hashing/hashing.hpp"

namespace libiop {

/** blake2b hash-chain. */
template<typename FieldT, typename MT_root_type>
class blake2b_hashchain : public hashchain<FieldT, MT_root_type>
{
    protected:
        binary_hash_digest internal_state_;
        const size_t security_parameter_;
        size_t digest_len_bytes_;
        size_t squeeze_index_ = 0;
    public:
        blake2b_hashchain(size_t security_parameter);
        void absorb(const MT_root_type new_input);
        /* internally does absorb(hash(new_input)) */
        void absorb(const std::vector<FieldT> &new_input);
        std::vector<FieldT> squeeze(const size_t num_elements);
        std::vector<size_t> squeeze_query_positions(
            const size_t num_positions, const size_t range_of_positions);

        MT_root_type squeeze_root_type();

        /* Needed for C++ polymorphism */
        std::shared_ptr<hashchain<FieldT, MT_root_type>> new_hashchain();
    protected:
        void absorb_hash_digest(const binary_hash_digest new_input);
        void absorb_internal(const typename libff::enable_if<std::is_same<MT_root_type, binary_hash_digest>::value, MT_root_type>::type new_input);
        void absorb_internal(const typename libff::enable_if<std::is_same<MT_root_type, FieldT>::value, MT_root_type>::type new_input);
};

template<typename FieldT>
class blake2b_leafhash : public leafhash<FieldT, binary_hash_digest>
{
    protected:
    size_t digest_len_bytes_;
    public:
    blake2b_leafhash(size_t security_parameter);
    binary_hash_digest hash(const std::vector<FieldT> &leaf);
    binary_hash_digest zk_hash(const std::vector<FieldT> &leaf,
        const zk_salt_type &zk_salt);
};

template<typename FieldT>
binary_hash_digest blake2b_field_element_hash(const std::vector<FieldT> &data,
                                       const std::size_t digest_len_bytes);

template<typename FieldT>
std::vector<FieldT> blake2b_FieldT_randomness_extractor(const binary_hash_digest &root,
                                                        const std::size_t index,
                                                        const std::size_t num_elements);

/* Returns a random integer of size less than upper_bound using input root, and key index */
std::size_t blake2b_integer_randomness_extractor(const binary_hash_digest &root,
                                                 const std::size_t index,
                                                 const std::size_t upper_bound);

binary_hash_digest blake2b_zk_element_hash(const std::vector<uint8_t> &first,
                                    const std::size_t digest_len_bytes);

binary_hash_digest blake2b_two_to_one_hash(const binary_hash_digest &first,
                                    const binary_hash_digest &second,
                                    const std::size_t digest_len_bytes);

} // namespace libiop

#include "libiop/bcs/hashing/blake2b.tcc"

#endif // LIBIOP_SNARK_COMMON_HASHING_BLAKE2B_HPP_
