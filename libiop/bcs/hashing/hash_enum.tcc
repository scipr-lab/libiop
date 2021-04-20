#include "libiop/bcs/hashing/hashing.hpp"
#include "libiop/bcs/hashing/algebraic_sponge.hpp"
#include "libiop/bcs/hashing/poseidon.hpp"
#include "libiop/bcs/hashing/blake2b.hpp"

#include <cstring>
#include <sstream>
#include <stdexcept>

namespace libiop {

template<typename FieldT>
poseidon_params<FieldT> get_poseidon_parameters(const bcs_hash_type hash_enum)
{
    if (hash_enum == starkware_poseidon_type)
    {
        return default_128_bit_altbn_poseidon_params<FieldT>();
    }
    else if (hash_enum == high_alpha_poseidon_type)
    {
        return high_alpha_128_bit_altbn_poseidon_params<FieldT>();
    }
    throw std::invalid_argument("Not a poseidon hash type");
}

/* Algebraic hashchain case */
template<typename FieldT, typename MT_root_type>
std::shared_ptr<hashchain<FieldT, MT_root_type>> get_hashchain_internal(
    const typename libff::enable_if<std::is_same<MT_root_type, FieldT>::value, FieldT>::type _, 
    const bcs_hash_type hash_enum,
    const size_t security_parameter)
{
    if (hash_enum == starkware_poseidon_type || hash_enum == high_alpha_poseidon_type)
    {
        if (security_parameter != 128)
        {
            throw std::invalid_argument("Poseidon only supported for 128 bit soundness.");
        }
        poseidon_params<FieldT> params = get_poseidon_parameters<FieldT>(hash_enum);

        std::shared_ptr<algebraic_sponge<FieldT>> permutation = std::make_shared<poseidon<FieldT>>(params);
        return std::make_shared<algebraic_hashchain<FieldT, MT_root_type>>(
            permutation, 
            security_parameter - 1);
    }
    throw std::invalid_argument("bcs_hash_type unknown (algebraic hashchain)");
}

/* Binary_hash_digest hashchain */
/* Algebraic hashchain case */
template<typename FieldT, typename MT_root_type>
std::shared_ptr<hashchain<FieldT, MT_root_type>> get_hashchain_internal(
    const typename libff::enable_if<!std::is_same<MT_root_type, FieldT>::value, FieldT>::type _, 
    const bcs_hash_type hash_enum,
    const size_t security_parameter)
{
    if (hash_enum == blake2b_type)
    {
        return std::make_shared<blake2b_hashchain<FieldT, MT_root_type>>(security_parameter);
    }
    throw std::invalid_argument("bcs_hash_type unknown");
}


template<typename FieldT, typename MT_root_type>
std::shared_ptr<hashchain<FieldT, MT_root_type>> get_hashchain(bcs_hash_type hash_enum, size_t security_parameter)
{
    return get_hashchain_internal<FieldT, MT_root_type>(FieldT::zero(), hash_enum, security_parameter);
}

/* Algebraic leafhash case */
template<typename FieldT, typename leaf_hash_type>
std::shared_ptr<leafhash<FieldT, leaf_hash_type>> get_leafhash_internal(
    const typename libff::enable_if<std::is_same<leaf_hash_type, FieldT>::value, FieldT>::type _, 
    const bcs_hash_type hash_enum,
    const size_t security_parameter, 
    const size_t leaf_size)
{
    if (hash_enum == starkware_poseidon_type || hash_enum == high_alpha_poseidon_type)
    {
        if (security_parameter != 128)
        {
            throw std::invalid_argument("Poseidon only supported for 128 bit soundness.");
        }
        poseidon_params<FieldT> params = get_poseidon_parameters<FieldT>(hash_enum);
        /* security parameter is -1 b/c */
        std::shared_ptr<poseidon<FieldT>> permutation = std::make_shared<poseidon<FieldT>>(params);
        std::shared_ptr<leafhash<FieldT, leaf_hash_type>> leafhasher = std::make_shared<algebraic_leafhash<FieldT>>(
            permutation, 
            security_parameter - 1);
        return leafhasher;
    }
    throw std::invalid_argument("bcs_hash_type unknown (algebraic leaf hash)");
}

/* Binary_hash_digest leafhash */
template<typename FieldT, typename leaf_hash_type>
std::shared_ptr<leafhash<FieldT, leaf_hash_type>> get_leafhash_internal(
    const typename libff::enable_if<!std::is_same<leaf_hash_type, FieldT>::value, FieldT>::type _, 
    const bcs_hash_type hash_enum,
    const size_t security_parameter, 
    const size_t leaf_size)
{
    if (hash_enum == blake2b_type)
    {
        return std::make_shared<blake2b_leafhash<FieldT>>(security_parameter);
    }
    throw std::invalid_argument("bcs_hash_type unknown");
}

template<typename FieldT, typename leaf_hash_type>
std::shared_ptr<leafhash<FieldT, leaf_hash_type>> get_leafhash(
    const bcs_hash_type hash_enum, const size_t security_parameter, const size_t leaf_size)
{
    return get_leafhash_internal<FieldT, leaf_hash_type>(FieldT::zero(), hash_enum, security_parameter, leaf_size);
}

/* binary hash digest 2->1 hash */
template<typename hash_type, typename FieldT>
two_to_one_hash_function<hash_type> get_two_to_one_hash_internal(
    const typename libff::enable_if<!std::is_same<hash_type, FieldT>::value, FieldT>::type _, 
    const bcs_hash_type hash_enum, 
    const size_t security_parameter)
{
    if (hash_enum == blake2b_type)
    {
        return blake2b_two_to_one_hash;
    }
    throw std::invalid_argument("bcs_hash_type unknown");
}

/* algebraic 2->1 hash */
template<typename hash_type, typename FieldT>
two_to_one_hash_function<FieldT> get_two_to_one_hash_internal(
    const typename libff::enable_if<std::is_same<hash_type, FieldT>::value, FieldT>::type _, 
    const bcs_hash_type hash_enum, 
    const size_t security_parameter)
{
    if (hash_enum == starkware_poseidon_type || hash_enum == high_alpha_poseidon_type)
    {
        if (security_parameter != 128)
        {
            throw std::invalid_argument("Poseidon only supported for 128 bit soundness.");
        }
        poseidon_params<FieldT> params = get_poseidon_parameters<FieldT>(hash_enum);
        /* security parameter is -1 b/c */
        std::shared_ptr<algebraic_sponge<FieldT>> permutation = std::make_shared<poseidon<FieldT>>(params);
        /* We explicitly place this on heap with no destructor,
           as this reference has to live after the function terminates */
        std::shared_ptr<algebraic_two_to_one_hash<FieldT>> hash_class = 
            std::make_shared<algebraic_two_to_one_hash<FieldT>>(permutation, security_parameter - 1);
        std::function<FieldT(const FieldT&, const FieldT&, const std::size_t)> f = [permutation, hash_class](const FieldT& left, const FieldT& right, const std::size_t unused) -> FieldT 
        {
            return hash_class->hash(left, right);
        };
        return f;
    }
    throw std::invalid_argument("bcs_hash_type unknown (algebraic two to one hash)");
}

template<typename hash_type, typename FieldT>
two_to_one_hash_function<hash_type> get_two_to_one_hash(const bcs_hash_type hash_enum, const size_t security_parameter)
{
    return get_two_to_one_hash_internal<hash_type, FieldT>(FieldT::zero(), hash_enum, security_parameter);
}

}
