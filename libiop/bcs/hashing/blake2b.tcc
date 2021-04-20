#include "sodium/crypto_generichash_blake2b.h"
#include <libff/algebra/field_utils/field_utils.hpp>
#include "libiop/bcs/hashing/hashing.hpp"
#include <cstring>
#include <sstream>
#include <stdexcept>

namespace libiop {

template<typename FieldT, typename hash_data_type>
blake2b_hashchain<FieldT, hash_data_type>::blake2b_hashchain(size_t security_parameter) :
    security_parameter_(security_parameter)
{
    /* 2*security_parameter bits, rounded up to next byte */
    this->digest_len_bytes_ = ((2*security_parameter) + 7) / 8;
    /* TODO: Should we personalize this? */
    this->internal_state_ = binary_hash_digest(this->digest_len_bytes_, ' ');
}

template<typename FieldT, typename hash_data_type>
std::shared_ptr<hashchain<FieldT, hash_data_type>>
    blake2b_hashchain<FieldT, hash_data_type>::new_hashchain()
{
    return std::make_shared<blake2b_hashchain<FieldT, hash_data_type>>
        (this->security_parameter_);
}

template<typename FieldT, typename hash_data_type>
void blake2b_hashchain<FieldT, hash_data_type>::absorb(const hash_data_type new_input)
{
    this->absorb_internal(new_input);
}

template<typename FieldT, typename hash_data_type>
void blake2b_hashchain<FieldT, hash_data_type>::absorb_internal(
    const typename libff::enable_if<std::is_same<hash_data_type, FieldT>::value, hash_data_type>::type new_input)
{
    std::vector<FieldT> vec_new_input;
    vec_new_input.emplace_back(new_input);
    this->absorb(vec_new_input);
}

template<typename FieldT, typename hash_data_type>
void blake2b_hashchain<FieldT, hash_data_type>::absorb_internal(
    const typename libff::enable_if<std::is_same<hash_data_type, binary_hash_digest>::value, hash_data_type>::type new_input)
{
    this->absorb_hash_digest(new_input);
}

template<typename FieldT, typename hash_data_type>
void blake2b_hashchain<FieldT, hash_data_type>::absorb_hash_digest(
    const binary_hash_digest new_input)
{
    binary_hash_digest hash_input = this->internal_state_ + new_input;

    /* see https://download.libsodium.org/doc/hashing/generic_hashing.html */
    const int status = crypto_generichash_blake2b((unsigned char*)&this->internal_state_[0],
                                                  this->digest_len_bytes_,
                                                  (unsigned char*)&hash_input[0],
                                                  this->digest_len_bytes_,
                                                  NULL, 0);
    if (status != 0)
    {
        throw std::runtime_error("Got non-zero status from crypto_generichash_blake2b. (Is digest_len_bytes correct?)");
    }
}

template<typename FieldT, typename hash_data_type>
void blake2b_hashchain<FieldT, hash_data_type>::absorb(const std::vector<FieldT> &new_input)
{
    const binary_hash_digest new_input_hash =
        blake2b_field_element_hash<FieldT>(new_input, this->digest_len_bytes_);
    this->absorb_hash_digest(new_input_hash);
}

template<typename FieldT, typename hash_data_type>
std::vector<FieldT> blake2b_hashchain<FieldT, hash_data_type>::squeeze(
    const size_t num_elements)
{
    this->squeeze_index_++;
    return blake2b_FieldT_randomness_extractor<FieldT>(
        this->internal_state_,
        this->squeeze_index_,
        num_elements);
}

template<typename FieldT, typename hash_data_type>
std::vector<size_t> blake2b_hashchain<FieldT, hash_data_type>::squeeze_query_positions(
        const size_t num_positions, const size_t range_of_positions)
{
    std::vector<size_t> query_pos;
    for (size_t i = 0; i < num_positions; i++)
    {
        this->squeeze_index_++;
        query_pos.emplace_back(
            blake2b_integer_randomness_extractor(
                this->internal_state_,
                this->squeeze_index_,
                range_of_positions)
        );
    }
    return query_pos;
}

template<typename FieldT, typename hash_data_type>
hash_data_type blake2b_hashchain<FieldT, hash_data_type>::squeeze_root_type()
{
    std::vector<FieldT> x = this->squeeze(1);
    return blake2b_field_element_hash<FieldT>(x, this->digest_len_bytes_);
}


template<typename FieldT>
blake2b_leafhash<FieldT>::blake2b_leafhash(size_t security_parameter)
{
    /* 2*security_parameter bits, rounded up to next byte */
    this->digest_len_bytes_ = ((2*security_parameter) + 7) / 8;
}

template<typename FieldT>
binary_hash_digest blake2b_leafhash<FieldT>::hash(const std::vector<FieldT> &leaf)
{
    return blake2b_field_element_hash<FieldT>(leaf, this->digest_len_bytes_);
}

template<typename FieldT>
binary_hash_digest blake2b_leafhash<FieldT>::zk_hash(
    const std::vector<FieldT> &leaf,
    const zk_salt_type &zk_salt)
{
    /* TODO: This is inefficient (requires 2 hash calls),
             this should be done with 1. */
    binary_hash_digest leaf_hash = blake2b_field_element_hash<FieldT>(
        leaf, this->digest_len_bytes_);
    return blake2b_two_to_one_hash(leaf_hash, zk_salt, this->digest_len_bytes_);
}

// TODO: Consider how this interacts with field elems being in montgomery form
// don't we need to make them in canonical form first?
template<typename FieldT>
binary_hash_digest blake2b_field_element_hash(const std::vector<FieldT> &data,
                                       const std::size_t digest_len_bytes)
{

    binary_hash_digest result(digest_len_bytes, 'X');

    /* see https://download.libsodium.org/doc/hashing/generic_hashing.html */
    const int status = crypto_generichash_blake2b((unsigned char*)&result[0],
                                                  digest_len_bytes,
                                                  (result.empty() ? NULL : (unsigned char*)&data[0]),
                                                  sizeof(FieldT) * data.size(),
                                                  NULL, 0);
    if (status != 0)
    {
        throw std::runtime_error("Got non-zero status from crypto_generichash_blake2b. (Is digest_len_bytes correct?)");
    }


    return result;
}

template<typename FieldT>
FieldT blake2b_FieldT_rejection_sample(
    typename libff::enable_if<libff::is_additive<FieldT>::value, FieldT>::type _,
    unsigned char* root_plus_index,
    size_t root_plus_index_size,
    const size_t key,
    const size_t key_increment)
{
    /* No need for rejection sampling, since our binary fields are word-aligned */
    FieldT el;
    const int status = crypto_generichash_blake2b((unsigned char*)&el,
                                                   sizeof(el),
                                                   root_plus_index,
                                                   root_plus_index_size,
                                                   (unsigned char*)&key, sizeof(key));
    if (status != 0)
    {
        throw std::runtime_error("Got non-zero status from crypto_generichash_blake2b. (Is digest_len_bytes correct?)");
    }
    return el;
}

template<typename FieldT>
FieldT blake2b_FieldT_rejection_sample(
    typename libff::enable_if<libff::is_multiplicative<FieldT>::value, FieldT>::type _,
    unsigned char* root_plus_index,
    size_t root_plus_index_size,
    const size_t key,
    const size_t key_increment)
{
    FieldT el;
    bool valid = false;
    size_t cur_key = key;
    const size_t bits_per_limb = 8 * sizeof(mp_limb_t);
    const size_t num_limbs = sizeof(el.mont_repr) / sizeof(mp_limb_t);
    while (!valid)
    {
        /* crypto generichash is keyed */
        const int status = crypto_generichash_blake2b((unsigned char*)&el.mont_repr,
                                                      sizeof(el.mont_repr),
                                                      root_plus_index,
                                                      root_plus_index_size,
                                                      (unsigned char*)&cur_key,
                                                      sizeof(cur_key));
        if (status != 0)
        {
            throw std::runtime_error("Got non-zero status from crypto_generichash_blake2b. (Is digest_len_bytes correct?)");
        }
        /* clear all bits higher than MSB of modulus */
        size_t bitno = sizeof(el.mont_repr) * 8 - 1;
        while (FieldT::mod.test_bit(bitno) == false)
        {
            const std::size_t part = bitno / bits_per_limb;
            const std::size_t bit = bitno - (bits_per_limb*part);

            el.mont_repr.data[part] &= ~(1ul<<bit);
            bitno--;
        }
        /* if el.data is < modulus its valid, otherwise repeat (rejection sampling) */
        if (mpn_cmp(el.mont_repr.data, FieldT::mod.data, num_limbs) < 0)
        {
            valid = true;
        }

        cur_key += key_increment;
    }
    return el;
}

template<typename FieldT>
std::vector<FieldT> blake2b_FieldT_randomness_extractor(const binary_hash_digest &root,
                                                        const std::size_t index,
                                                        const std::size_t num_elements)
{
    const std::size_t root_plus_index_size = root.size() + sizeof(index);
    unsigned char* root_plus_index = (unsigned char*)(malloc(root_plus_index_size));
    memcpy(root_plus_index, &root[0], root.size());
    memcpy(root_plus_index + root.size(), &index, sizeof(index));

    std::vector<FieldT> result;
    result.reserve(num_elements);

    for (std::size_t i = 0; i < num_elements; ++i)
    {
        FieldT el = blake2b_FieldT_rejection_sample<FieldT>(
            FieldT::zero(),
            root_plus_index, root_plus_index_size,
            i, num_elements);

        result.emplace_back(el);
    }

    free(root_plus_index);

    return result;
}

}
