/**@file
 *****************************************************************************
 Hash function declarations and implementations.
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_SNARK_COMMON_HASHING_ALGEBRAIC_SPONGE_HPP_
#define LIBIOP_SNARK_COMMON_HASHING_ALGEBRAIC_SPONGE_HPP_

#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include "libiop/bcs/hashing/hashing.hpp"
#include <libff/algebra/field_utils/field_utils.hpp>

namespace libiop {

/** An abstract class for algebraic sponges
 *  TODO: Benchmark if state size should be fixed at compile time.
 *  TODO: Benchmark if capacity should be fixed at compile time, like libFq.
 *  One could conceivably expand rate, I think that will have to be a separate class.
 *  (This performance difference does matter)
 *
 *  Note that this is causing the hash output digest size to be fixed at compile time.
 * */
template<typename FieldT>
class algebraic_sponge
{
    protected:
    std::vector<FieldT> state_;

    bool currently_absorbing = false;
    /** This is to track what element in the rate to output next */
    size_t next_unsqueezed_elem_ = 0;
    /* This is set based on your hash of choice. */
    virtual void apply_permutation() {};
    void absorb_internal(const std::vector<FieldT> &new_input, size_t begin_index);
    void squeeze_internal(std::vector<FieldT> &output, size_t begin_index);

    algebraic_sponge<FieldT>(const size_t rate, const size_t capacity);
    public:
    const size_t rate_;
    const size_t capacity_;

    virtual double achieved_security_parameter() const { return 0.0; };
    virtual void print() const {};
    // void absorb(const FieldT[state_size] new_input);
    void absorb(const std::vector<FieldT> &new_input);   
    std::vector<FieldT> squeeze_vector(size_t num_elements);

    /* Only use in two to one hash, not black-box sponge functionality. */
    void initialize_element_of_state(const FieldT elem, const size_t index);

    /* Needed for C++ polymorphism*/
    virtual void reset() { printf("top-level reset\n"); };
    virtual std::shared_ptr<algebraic_sponge<FieldT>> new_sponge() { return NULL; };
    
    virtual ~algebraic_sponge() {};
};

template<typename FieldT>
FieldT string_to_field_elem(typename libff::enable_if<libff::is_multiplicative<FieldT>::value, FieldT>::type dummy_field_elem,
    const zk_salt_type &zk_salt);
template<typename FieldT>
FieldT string_to_field_elem(typename libff::enable_if<libff::is_additive<FieldT>::value, FieldT>::type dummy_field_elem,
    const zk_salt_type &zk_salt);

template<typename FieldT, typename MT_root_type>
class algebraic_hashchain : public hashchain<FieldT, MT_root_type>
{
    protected:
    std::shared_ptr<algebraic_sponge<FieldT>> sponge_;
    size_t security_parameter_;
    public:
    algebraic_hashchain(
        std::shared_ptr<algebraic_sponge<FieldT>> sponge,
        size_t security_parameter);

    void absorb(const MT_root_type new_input);
    void absorb(const std::vector<FieldT> &new_input);
    std::vector<FieldT> squeeze(size_t num_elements);
    std::vector<size_t> squeeze_query_positions(
        size_t num_positions, size_t range_of_positions);

    // TODO: implement for binary MT root type
    MT_root_type squeeze_root_type();

    /* Needed for C++ polymorphism*/
    std::shared_ptr<hashchain<FieldT, MT_root_type>> new_hashchain() {
        return std::make_shared<algebraic_hashchain<FieldT, MT_root_type>>(
            this->sponge_->new_sponge(), security_parameter_);
    };
    protected:
    void absorb_internal(const typename libff::enable_if<std::is_same<MT_root_type, binary_hash_digest>::value, MT_root_type>::type new_input);
    void absorb_internal(const typename libff::enable_if<std::is_same<MT_root_type, FieldT>::value, MT_root_type>::type new_input);
};

template<typename FieldT>
class algebraic_leafhash : public leafhash<FieldT, FieldT>
{
    protected:
    std::shared_ptr<algebraic_sponge<FieldT>> sponge_;

    public:
    algebraic_leafhash(
        std::shared_ptr<algebraic_sponge<FieldT>> sponge,
        size_t security_parameter);
    FieldT hash(const std::vector<FieldT> &leaf);
    FieldT zk_hash(
        const std::vector<FieldT> &leaf,
        const zk_salt_type &zk_salt);
};

template<typename FieldT>
class algebraic_two_to_one_hash
{
    protected:

    public:
    std::shared_ptr<algebraic_sponge<FieldT>> sponge_;
    algebraic_two_to_one_hash(
        std::shared_ptr<algebraic_sponge<FieldT>> sponge,
        size_t security_parameter);
    FieldT hash(const FieldT &left, const FieldT &right);
};

// This is intended to be compatible with STARKWARE's MDS generation
// template<typename FieldT>
// std::vector<std::vector<FieldT>> generate_mds_matrix(std::string name, size_t state_size);

} // namespace libiop

#include "libiop/bcs/hashing/algebraic_sponge.tcc"

#endif // LIBIOP_SNARK_COMMON_HASHING_ALGEBRAIC_SPONGE_HPP_
