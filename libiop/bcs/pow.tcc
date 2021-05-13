#include <algorithm>
#include <stdexcept>

#include <libff/common/profiling.hpp>
#include "libiop/common/cpp17_bits.hpp"
#include <libff/common/utils.hpp>

#include <sodium/randombytes.h>
#include <libff/algebra/field_utils/field_utils.hpp>

namespace libiop {

pow_parameters::pow_parameters(
    const size_t work_parameter,
    const size_t cost_per_hash) :
    work_parameter_(work_parameter),
    cost_per_hash_(cost_per_hash)
{
}

size_t pow_parameters::pow_bitlen() const
{
    // For now we round the hash cost to 1 << floor(libff::log2(cost))
    // This makes the proof of work condition very simple, at the expense of some extra prover work
    // for the same desired security.
    size_t log_hash_cost = libff::log2(this->cost_per_hash_);
    if ((1 << log_hash_cost) > this->cost_per_hash_)
    {
        log_hash_cost -= 1;
    }
    return this->work_parameter_ - log_hash_cost;
}

// For now we don't implement the optimization for non-power-of-2 hash_costs,
// so this is always set to 0
size_t pow_parameters::pow_upperbound() const
{
    return 0;
}

size_t pow_parameters::work_parameter() const
{
    return this->work_parameter_;
}


void pow_parameters::print() const
{
    printf("\nProof of work parameters\n");
    libff::print_indent(); printf("* log of target work amount = %zu\n", this->work_parameter_);
    libff::print_indent(); printf("* assumed cost per hash = %zu\n", this->cost_per_hash_);
    libff::print_indent(); printf("* expected number of hashes = %d\n", 1 << this->pow_bitlen());
}

template<typename FieldT, typename hash_digest_type>
pow<FieldT, hash_digest_type>::pow(
    const pow_parameters params,
    const size_t digest_len_bytes) :
    parameters_(params),
    digest_len_bytes_(digest_len_bytes)
{
}


template<typename FieldT, typename hash_digest_type>
hash_digest_type pow<FieldT, hash_digest_type>::solve_pow(
    const two_to_one_hash_function<hash_digest_type> &node_hasher, 
    const hash_digest_type &challenge) const
{
    return this->solve_pow_internal(node_hasher, challenge);
}

template<typename FieldT, typename hash_digest_type>
hash_digest_type pow<FieldT, hash_digest_type>::solve_pow_internal(
    const two_to_one_hash_function<hash_digest_type> &node_hasher, 
    const typename libff::enable_if<std::is_same<hash_digest_type, FieldT>::value, hash_digest_type>::type challenge) const
{
    FieldT pow = FieldT::zero();
    while (this->verify_pow(node_hasher, challenge, pow) == false)
    {
        pow += FieldT::one();
    }
    return pow;
}

template<typename FieldT, typename hash_digest_type>
hash_digest_type pow<FieldT, hash_digest_type>::solve_pow_internal(
    const two_to_one_hash_function<hash_digest_type> &node_hasher, 
    const typename libff::enable_if<std::is_same<hash_digest_type, binary_hash_digest>::value, hash_digest_type>::type challenge) const
{
    binary_hash_digest pow;
    pow.assign(challenge);

    size_t num_words = pow.length() / sizeof(size_t);
    size_t pow_int = 0;
    while (this->verify_pow(node_hasher, challenge, pow) == false)
    {
        std::memcpy(&pow[(num_words - 1)*sizeof(size_t)], &pow_int, sizeof(size_t));
        pow_int += 1;
    }
    return pow;
}

template<typename FieldT, typename hash_digest_type>
bool pow<FieldT, hash_digest_type>::verify_pow(
    const two_to_one_hash_function<hash_digest_type> &node_hasher, 
    const hash_digest_type &challenge,
    const hash_digest_type &pow) const
{
    hash_digest_type hash = node_hasher(challenge, pow, this->digest_len_bytes_);
    return this->verify_pow_internal(hash);
}

// Function to sanity check the PoW's.
void print_string_in_hex(const std::string& input)
{
    static const char hex_digits[] = "0123456789ABCDEF";

    std::string output;
    output.reserve(input.length() * 2);
    for (unsigned char c : input)
    {
        output.push_back(hex_digits[c >> 4]);
        output.push_back(hex_digits[c & 15]);
    }
    std::cout << output << "\n";
}

template<typename FieldT, typename hash_digest_type>
bool pow<FieldT, hash_digest_type>::verify_pow_internal(
    const typename libff::enable_if<std::is_same<hash_digest_type, FieldT>::value, hash_digest_type>::type &hash) const
{
    size_t least_significant_word = libff::get_word_of_field_elem<FieldT>(hash, 0);
    size_t relevant_bits = least_significant_word & ((1 << this->parameters_.pow_bitlen()) - 1);
    if (relevant_bits <= this->parameters_.pow_upperbound())
    {
        return true;
    }
    return false;
}

template<typename FieldT, typename hash_digest_type>
bool pow<FieldT, hash_digest_type>::verify_pow_internal(
    const typename libff::enable_if<std::is_same<hash_digest_type, binary_hash_digest>::value, hash_digest_type>::type &hash) const
{
    size_t num_words = hash.length() / sizeof(size_t);
    size_t least_significant_word;
    std::memcpy(&least_significant_word, &hash[(num_words - 1)*sizeof(size_t)], sizeof(size_t));
    size_t relevant_bits = least_significant_word & ((1 << this->parameters_.pow_bitlen()) - 1);
    if (relevant_bits <= this->parameters_.pow_upperbound())
    {    
        // printf("%d\n", (1 << this->parameters_.pow_bitlen()));
        // printf("%zu\n", least_significant_word);
        // printf("%\n", relevant_bits);
        // print_string_in_hex(hash);
        return true;
    }
    return false;
}

} // libiop
