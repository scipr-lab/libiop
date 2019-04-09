/**@file
 *****************************************************************************
 Classes for multiplicative subgroups / cosets of order 2^n
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_ALGEBRA_FIELD_SUBSET_SUBGROUP_HPP_
#define LIBIOP_ALGEBRA_FIELD_SUBSET_SUBGROUP_HPP_

#include <cstddef>
#include <vector>
#include <cstdint>
#include <libfqfft/evaluation_domain/get_evaluation_domain.hpp>
#include <libfqfft/evaluation_domain/domains/basic_radix2_domain.hpp>

#include "libiop/algebra/fields/utils.hpp"

namespace libiop {

template<typename FieldT>
class multiplicative_subgroup_base {
protected:
    std::shared_ptr<std::vector<FieldT>> elems_;

    FieldT g_;
    u_long order_;

    std::shared_ptr<libfqfft::basic_radix2_domain<FieldT>> FFT_eval_domain_;

public:
    multiplicative_subgroup_base() = default;
    multiplicative_subgroup_base(const multiplicative_subgroup_base<FieldT> &other) = default;
    multiplicative_subgroup_base(multiplicative_subgroup_base<FieldT> &&other) = default;
    multiplicative_subgroup_base<FieldT>& operator=(const multiplicative_subgroup_base<FieldT> &other) = default;
    multiplicative_subgroup_base<FieldT>& operator=(multiplicative_subgroup_base<FieldT> &&other) = default;

    multiplicative_subgroup_base(FieldT order); // FIXME: Do we need a basis? Don't think so because we just need the dimension-th roots of unity as a basis for H (?)
    multiplicative_subgroup_base(std::size_t order);

    FieldT generator() const;
    u_long order() const; // FIXME: Does it work as u_long? Or should it be a FieldT/bigint?

    std::size_t dimension() const;
    std::size_t num_elements() const;

    virtual std::vector<FieldT> all_elements() const = 0;
    virtual FieldT element_by_index(const std::size_t index) const = 0;
    std::size_t reindex_by_subgroup(const std::size_t reindex_subgroup_dim, const std::size_t index) const;
    std::size_t coset_index(const std::size_t position, const std::size_t coset_size) const;
    std::size_t intra_coset_index(const std::size_t position, const std::size_t coset_size) const;
    size_t position_by_coset_indices(
        const size_t coset_index, const size_t intra_coset_index, const size_t coset_size) const;

    libfqfft::basic_radix2_domain<FieldT> FFT_eval_domain() const;

protected:
    void construct_internal(typename libiop::enable_if<is_multiplicative<FieldT>::value, FieldT>::type order);
    void construct_internal(typename libiop::enable_if<is_additive<FieldT>::value, FieldT>::type order);
};

template<typename FieldT>
class multiplicative_subgroup : public multiplicative_subgroup_base<FieldT> {
    public:
    using multiplicative_subgroup_base<FieldT>::multiplicative_subgroup_base;
    std::vector<FieldT> all_elements() const;
    FieldT element_by_index(const std::size_t index) const;
};

template<typename FieldT>
class multiplicative_coset : public multiplicative_subgroup_base<FieldT> {
protected:
    FieldT shift_;

public:
    multiplicative_coset() = default;
    multiplicative_coset(const multiplicative_coset<FieldT> &other) = default;
    multiplicative_coset(multiplicative_coset<FieldT> &&other) = default;
    multiplicative_coset<FieldT>& operator=(const multiplicative_coset<FieldT> &other) = default;
    multiplicative_coset<FieldT>& operator=(multiplicative_coset<FieldT> &&other) = default;

    multiplicative_coset(FieldT order);
    multiplicative_coset(std::size_t order);
    multiplicative_coset(FieldT order, FieldT shift);
    multiplicative_coset(std::size_t order, FieldT shift);

    std::vector<FieldT> all_elements() const;
    FieldT element_by_index(const std::size_t index) const;

    FieldT shift() const;
};



}

#include "libiop/algebra/field_subset/subgroup.tcc"

#endif // LIBIOP_ALGEBRA_FIELD_SUBSET_SUBGROUP_HPP_
