#include "sodium/crypto_generichash_blake2b.h"
#include <libff/algebra/field_utils/field_utils.hpp>
#include "libiop/bcs/hashing/hashing.hpp"
#include <cstring>
#include <sstream>
#include <ostream>
#include <stdexcept>

namespace libiop {

/* TODO: Move dummy hashchain to its own file */
template<typename FieldT, typename hash_data_type>
dummy_algebraic_hashchain<FieldT, hash_data_type>::dummy_algebraic_hashchain()
{
    this->internal_state_ = FieldT::zero();
}

template<typename FieldT, typename hash_data_type>
std::shared_ptr<hashchain<FieldT, hash_data_type>>
    dummy_algebraic_hashchain<FieldT, hash_data_type>::new_hashchain()
{
    return std::make_shared<dummy_algebraic_hashchain<FieldT, hash_data_type>>();
}

template<typename FieldT, typename hash_data_type>
void dummy_algebraic_hashchain<FieldT, hash_data_type>::absorb(const hash_data_type new_input)
{
    this->absorb_internal(new_input);
}

template<typename FieldT, typename hash_data_type>
void dummy_algebraic_hashchain<FieldT, hash_data_type>::absorb_internal(
    const typename libff::enable_if<std::is_same<hash_data_type, FieldT>::value, hash_data_type>::type new_input)
{
    this->internal_state += new_input;
}

template<typename FieldT, typename hash_data_type>
void dummy_algebraic_hashchain<FieldT, hash_data_type>::absorb_internal(
    const typename libff::enable_if<std::is_same<hash_data_type, binary_hash_digest>::value, hash_data_type>::type new_input)
{
    std::stringstream ss(new_input);
    int64_t x = 0;
    while(ss >> x)           // get an int64 value from string stream iss
    {
        this->internal_state_ += FieldT(x);
    }
}

template<typename FieldT, typename hash_data_type>
void dummy_algebraic_hashchain<FieldT, hash_data_type>::absorb(const std::vector<FieldT> &new_input)
{
    for (size_t i = 0; i < new_input.size(); i++)
    {
        this->internal_state_ += new_input[i];
    }
}

template<typename FieldT, typename hash_data_type>
std::vector<FieldT> dummy_algebraic_hashchain<FieldT, hash_data_type>::squeeze(
    const size_t num_elements)
{
    std::vector<FieldT> squeezed_outputs;
    for (size_t i = 0; i < num_elements; i++)
    {
        this->squeeze_index_++;
        FieldT x = FieldT(this->squeeze_index_) + this->internal_state_;
        squeezed_outputs.emplace_back(x);
    }
    return squeezed_outputs;
}

template<typename FieldT, typename hash_data_type>
std::vector<size_t> dummy_algebraic_hashchain<FieldT, hash_data_type>::squeeze_query_positions(
        const size_t num_positions, const size_t range_of_positions)
{
    std::vector<size_t> query_pos;
    for (size_t i = 0; i < num_positions; i++)
    {
        this->squeeze_index_++;
        query_pos.emplace_back(
            this->squeeze_index_ % range_of_positions
        );
    }
    return query_pos;
}

template<typename FieldT, typename hash_data_type>
hash_data_type dummy_algebraic_hashchain<FieldT, hash_data_type>::squeeze_root_type()
{
    hash_data_type dummy;
    return this->squeeze_root_type_internal(dummy);
}

template<typename FieldT, typename hash_data_type>
hash_data_type dummy_algebraic_hashchain<FieldT, hash_data_type>::squeeze_root_type_internal(
    const typename libff::enable_if<std::is_same<hash_data_type, FieldT>::value, hash_data_type>::type dummy)
{
    return this->squeeze(1)[0];
}

template<typename FieldT, typename hash_data_type>
hash_data_type dummy_algebraic_hashchain<FieldT, hash_data_type>::squeeze_root_type_internal(
    const typename libff::enable_if<std::is_same<hash_data_type, binary_hash_digest>::value, hash_data_type>::type dummy)
{
    FieldT h = this->squeeze(1)[0];
    std::stringstream ss;
    // TODO: Fix libff import issue here
    // ss << h;
    binary_hash_digest s = ss.str();
    return s;
}

template<typename FieldT>
FieldT dummy_algebraic_leafhash<FieldT>::hash(const std::vector<FieldT> &leaf)
{
    FieldT sum = FieldT::zero();
    for (size_t i = 0; i < leaf.size(); i++)
    {
        sum += FieldT(i) * leaf[i];
    }
    return sum;
}

template<typename FieldT>
FieldT dummy_algebraic_leafhash<FieldT>::zk_hash(
    const std::vector<FieldT> &leaf,
    const zk_salt_type &zk_salt)
{
    FieldT leaf_hash = this->hash(leaf);
    // bad method to incorporate zk salt
    std::vector<uint8_t> salt_in_uint8(zk_salt.begin(), zk_salt.end());
    for (size_t i = 0; i < salt_in_uint8.size(); i++)
    {
        leaf_hash += FieldT((long)salt_in_uint8[i]);
    }
    return leaf_hash;
}

template<typename FieldT>
FieldT dummy_algebraic_two_to_one_hash(
    const FieldT &first,
    const FieldT &second,
    const std::size_t digest_len_bytes)
{
    // Need it to be non commutative for tests
    return FieldT(2) * first + second;
}

}
