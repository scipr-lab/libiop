/**@file
 *****************************************************************************
 Dummy algebraic hash function implementation
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_SNARK_COMMON_HASHING_DUMMY_ALGEBRAIC_HASH_HPP_
#define LIBIOP_SNARK_COMMON_HASHING_DUMMY_ALGEBRAIC_HASH_HPP_

#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>
#include <ostream>

#include "libiop/bcs/hashing/hashing.hpp"

namespace libiop {

/** This is an insecure hashchain,
 *  just meant to test if everything executes correctly with an algebraic hashchain */
template<typename FieldT, typename MT_root_type>
class dummy_algebraic_hashchain : public hashchain<FieldT, MT_root_type>
{
    protected:
        FieldT internal_state_;
        size_t squeeze_index_ = 0;
    public:
        dummy_algebraic_hashchain();
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
        void absorb_internal(const typename libff::enable_if<std::is_same<MT_root_type, binary_hash_digest>::value, MT_root_type>::type new_input);
        void absorb_internal(const typename libff::enable_if<std::is_same<MT_root_type, FieldT>::value, MT_root_type>::type new_input);

        MT_root_type squeeze_root_type_internal(const typename libff::enable_if<std::is_same<MT_root_type, binary_hash_digest>::value, MT_root_type>::type dummy);
        MT_root_type squeeze_root_type_internal(const typename libff::enable_if<std::is_same<MT_root_type, FieldT>::value, MT_root_type>::type dummy);
};

template<typename FieldT>
class dummy_algebraic_leafhash : public leafhash<FieldT, FieldT>
{
    public:
    dummy_algebraic_leafhash() {};
    FieldT hash(const std::vector<FieldT> &leaf);
    FieldT zk_hash(const std::vector<FieldT> &leaf,
        const zk_salt_type &zk_salt);
};

template<typename FieldT>
FieldT dummy_algebraic_two_to_one_hash(
    const FieldT &first,
    const FieldT &second,
    const std::size_t digest_len_bytes);

} // namespace libiop

#include "libiop/bcs/hashing/dummy_algebraic_hash.tcc"

#endif // LIBIOP_SNARK_COMMON_HASHING_DUMMY_ALGEBRAIC_HASH_HPP_
