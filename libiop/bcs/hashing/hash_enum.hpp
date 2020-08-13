/**@file
 *****************************************************************************
    Enum for all supported hash types to give to the SNARK interface
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_SNARK_COMMON_HASHING_HASH_ENUM_HPP_
#define LIBIOP_SNARK_COMMON_HASHING_HASH_ENUM_HPP_

#include "libiop/bcs/hashing/hashing.hpp"
#include "libiop/bcs/hashing/algebraic_sponge.hpp"
#include "libiop/bcs/hashing/poseidon.hpp"
#include "libiop/bcs/hashing/blake2b.hpp"

#include <memory>

namespace libiop {

/** One different enum element for each qualitatively different parameterization of a hash */
enum bcs_hash_type {
    blake2b_type = 1,
    /* Poseidon with MDS matrices and round selection identical to Starkware's choices */
    starkware_poseidon_type = 2,
    high_alpha_poseidon_type = 3
};

static const char* bcs_hash_type_names[] = {"", "blake2b", "poseidon with Starkware's parameterization", "poseidon with high alpha"};

template<typename FieldT, typename MT_root_type>
std::shared_ptr<hashchain<FieldT, MT_root_type>> get_hashchain(bcs_hash_type hash_type, size_t security_parameter);

template<typename FieldT, typename leaf_hash_type>
std::shared_ptr<leafhash<FieldT, leaf_hash_type>> get_leafhash(
    const bcs_hash_type hash_type, 
    const size_t security_parameter, 
    const size_t leaf_size);

template<typename hash_type, typename FieldT>
two_to_one_hash_function<hash_type> get_two_to_one_hash(const bcs_hash_type hash_enum, const size_t security_parameter);

}
#include "libiop/bcs/hashing/hash_enum.tcc"

#endif // LIBIOP_SNARK_COMMON_HASHING_HASH_ENUM_HPP_