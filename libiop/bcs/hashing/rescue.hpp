/**@file
 *****************************************************************************
    An implementation of Poseidon with the exponentiation S-Box
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_SNARK_COMMON_HASHING_POSEIDON_HPP_
#define LIBIOP_SNARK_COMMON_HASHING_POSEIDON_HPP_

#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>
#include "libiop/bcs/hashing/algebraic_sponge.hpp"

namespace libiop {

template<typename FieldT>
class rescue_params
{
    public:
    const size_t rounds_;
    const size_t alpha_;
    const std::vector<size_t> one_over_alpha_;
    const size_t state_size_;
    const size_t rate_;
    const size_t capacity_;
    // The ARK matrix is the "Add-Round-Key" matrix
    std::vector<std::vector<FieldT>> ark_matrix_;
    // The MDS matrix is any maximum distance separating matrix
    // We support special near-MDS matrices for state sizes 3 and 4 as well,
    // since this provides much better constraint complexity / evaluation time
    // with no down-sides on large fields.
    const bool supported_near_mds_ = false;
    std::vector<std::vector<FieldT>> mds_matrix_;

    rescue_params() {};
    /* We assume that these matrices satisfy the required properties. */
    rescue_params(size_t full_rounds, size_t partial_rounds, size_t alpha, size_t rate,
        std::vector<std::vector<FieldT>> &ark_matrix,
        bool supported_near_mds,
        std::vector<std::vector<FieldT>> &mds_matrix);

    double achieved_soundness() const;
    void print() const;
};

template<typename FieldT>
class rescue : public algebraic_sponge<FieldT>
{
    protected:
    const poseidon_params<FieldT> params_;
    const FieldT zero_singleton_;
    const FieldT a_;
    std::vector<FieldT> scratch_state_;

    /* Really should be generated at compile time */
    FieldT raise_to_alpha(const FieldT x) const;
    void apply_mix_layer();
    void apply_full_round(size_t round_id);
    void apply_ideal_partial_round(size_t round_id);
    void apply_partial_round(size_t round_id);
    void apply_permutation();
    public:
    poseidon(const poseidon_params<FieldT> params);

    double achieved_security_parameter() const;
    void print() const;
    /* TODO: Figure out right way to make this compile */
    // FieldT[state_size - capacity_] squeeze_rate();
    void reset();
    /* Needed for C++ polymorphism*/
    std::shared_ptr<algebraic_sponge<FieldT>> new_sponge() { 
        return std::make_shared<poseidon<FieldT>>(this->params_); 
    };
    ~poseidon() = default;
};

} // namespace libiop

#include "libiop/bcs/hashing/rescue.tcc"

#endif // LIBIOP_SNARK_COMMON_HASHING_POSEIDON_HPP_
