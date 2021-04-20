/**@file
 *****************************************************************************
 Proof of Work implementation for usage in BCS
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_SNARK_COMMON_BCS16_POW_HPP_
#define LIBIOP_SNARK_COMMON_BCS16_POW_HPP_

#include <algorithm>
#include <cstddef>
#include <map>
#include <vector>

#include "libiop/bcs/hashing/hashing.hpp"
#include "libiop/bcs/hashing/hash_enum.hpp"

namespace libiop {

class pow_parameters {
protected:
    std::size_t work_parameter_; /* The prover will do an expected 2^work_parameter units of work */
    std::size_t cost_per_hash_; /* How many units of work is one hash */
public:
    pow_parameters() {};
    pow_parameters(
        const size_t work_parameter,
        const size_t cost_per_hash);

    // The proof of work is satisfied if H(x) & pow_bitlen <= pow_upperbound
    size_t pow_bitlen() const;
    size_t pow_upperbound() const;
    size_t work_parameter() const;

    void print() const;
};

template<typename FieldT, typename hash_digest_type>
class pow {
protected:
    pow_parameters parameters_;
    size_t digest_len_bytes_;
public:
    pow() {};
    pow(const pow_parameters params,
        const size_t digest_len_bytes);

    // The proof of work is satisfied if H(x) & pow_bitlen < pow_upperbound
    // For binary hashes, this is done by interpreting H(x) as 4 words, each word being written little-endian.
    // This property is satisfied on the final little-endian word.
    hash_digest_type solve_pow(
        const two_to_one_hash_function<hash_digest_type> &node_hasher, 
        const hash_digest_type &challenge) const;

    bool verify_pow(
        const two_to_one_hash_function<hash_digest_type> &node_hasher, 
        const hash_digest_type &challenge,
        const hash_digest_type &pow) const;
protected:
    hash_digest_type solve_pow_internal(
        const two_to_one_hash_function<hash_digest_type> &node_hasher, 
        const typename libff::enable_if<std::is_same<hash_digest_type, binary_hash_digest>::value, hash_digest_type>::type challenge) const;
    hash_digest_type solve_pow_internal(
        const two_to_one_hash_function<hash_digest_type> &node_hasher, 
        const typename libff::enable_if<std::is_same<hash_digest_type, FieldT>::value, hash_digest_type>::type challenge) const;

    bool verify_pow_internal(
        const typename libff::enable_if<std::is_same<hash_digest_type, FieldT>::value, hash_digest_type>::type &hash) const;

    bool verify_pow_internal(
        const typename libff::enable_if<std::is_same<hash_digest_type, binary_hash_digest>::value, hash_digest_type>::type &hash) const;
};

} // namespace libiop

#include "libiop/bcs/pow.tcc"

#endif // LIBIOP_SNARK_COMMON_BCS16_POW_HPP_
