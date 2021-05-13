/**@file
 *****************************************************************************
 Hash function declarations and implementations.
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_SNARK_COMMON_HASHING_HASHING_HPP_
#define LIBIOP_SNARK_COMMON_HASHING_HASHING_HPP_

#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>
#include <libff/algebra/field_utils/field_utils.hpp>

namespace libiop {

// TODO: Document that this is marginally memory inefficient
typedef std::string binary_hash_digest;
typedef std::string zk_salt_type;



/* An abstract class for stateful hash-chains */
template<typename FieldT, typename MT_root_type>
class hashchain
{
    public:
    virtual void absorb(const MT_root_type new_input) = 0;
    virtual void absorb(const std::vector<FieldT> &new_input) = 0;
    virtual std::vector<FieldT> squeeze(size_t num_elements) = 0;
    virtual std::vector<size_t> squeeze_query_positions(
        size_t num_positions, size_t range_of_positions) = 0;

    virtual MT_root_type squeeze_root_type() = 0;
    /* Needed for C++ polymorphism*/
    virtual std::shared_ptr<hashchain<FieldT, MT_root_type>> new_hashchain() = 0;
};

/* An abstract class for leaf hashes */
template<typename FieldT, typename leaf_hash_type>
class leafhash
{
    public:
    virtual leaf_hash_type hash(const std::vector<FieldT> &leaf) = 0;
    virtual leaf_hash_type zk_hash(const std::vector<FieldT> &leaf,
        const zk_salt_type &zk_salt) = 0;
};

template<typename hash_type>
using two_to_one_hash_function = std::function<hash_type(const hash_type&, const hash_type&, const std::size_t)>;

/* Sizeof algebraic hash */
template<typename hash_type>
size_t get_hash_size(const typename libff::enable_if<!std::is_same<hash_type, binary_hash_digest>::value, hash_type>::type h)
{
    const size_t field_size =
        (libff::log_of_field_size_helper<hash_type>(hash_type::zero()) + 7) / 8;
    return field_size;
}

/* Sizeof binary hash */
template<typename hash_type>
size_t get_hash_size(const typename libff::enable_if<std::is_same<hash_type, binary_hash_digest>::value, hash_type>::type h)
{
    return h.size();
}


template<typename FieldT>
class hash_circuit_description
{
public:
    /* leaf hash and 2 to 1 hash complexities */
    static size_t arity_m_hash_complexity(
        const size_t m);
    static size_t hash_chain_complexity(
        const size_t sponge_state_size,
        const size_t input_size);
};

} // namespace libiop

#endif // LIBIOP_SNARK_COMMON_HASHING_HASHING_HPP_
