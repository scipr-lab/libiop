#include <algorithm>
#include <cstdint>
#include <iostream>

#include "libiop/common/common.hpp"

namespace libiop {

template<typename FieldT>
multiplicative_subgroup_base<FieldT>::multiplicative_subgroup_base(FieldT order)
{
    this->construct_internal(order);
}

template<typename FieldT>
multiplicative_subgroup_base<FieldT>::multiplicative_subgroup_base(std::size_t order)
{
    FieldT f(order);
    this->construct_internal(f);
}

template<typename FieldT>
void multiplicative_subgroup_base<FieldT>::construct_internal(typename libiop::enable_if<is_multiplicative<FieldT>::value, FieldT>::type order)
{
    FieldT F_order = FieldT(FieldT::mod) - 1;

    size_t order_as_size_t = order.as_bigint().data[0];
    if(!is_power_of_2(order_as_size_t))
    {
        // TODO: Support creating subgroups of other orders dividing the field size.
        throw std::invalid_argument("The order of the subgroup must be a power of two.");
    }

    this->g_ = (FieldT::multiplicative_generator)^((F_order * order.inverse()).as_bigint());

    this->elems_ = std::make_shared<std::vector<FieldT> >();
    this->order_ = order.as_ulong();

    if (is_power_of_2(this->order_))
    {
        this->FFT_eval_domain_ = std::make_shared<libfqfft::basic_radix2_domain<FieldT>>(this->order_);
    }
}

template<typename FieldT>
void multiplicative_subgroup_base<FieldT>::construct_internal(typename libiop::enable_if<is_additive<FieldT>::value, FieldT>::type order)
{
    throw std::invalid_argument("cannot create multiplicative subgroup of this field type");
}

template<typename FieldT>
FieldT multiplicative_subgroup_base<FieldT>::generator() const
{
    return this->g_;
}

template<typename FieldT>
u_long multiplicative_subgroup_base<FieldT>::order() const
{
    return this->order_;
}

template<typename FieldT>
std::size_t multiplicative_subgroup_base<FieldT>::dimension() const
{
    std::size_t num_elems = this->num_elements();

    /* dimension for a subgroup is only defined when the number of elements
       is a power of 2 */
    assert(is_power_of_2(num_elems));

    return log2(num_elems);
}

template<typename FieldT>
std::size_t multiplicative_subgroup_base<FieldT>::num_elements() const
{
    return (std::size_t) this->order_;
}

/** Given an index which assumes the first elements of this subgroup are the elements of
 *  another subgroup with dimension reindex_subgroup_dim,
 *  this returns the actual index into this subgroup. */
template<typename FieldT>
std::size_t multiplicative_subgroup_base<FieldT>::reindex_by_subgroup(const std::size_t reindex_subgroup_dim,
                                                            const std::size_t index) const
{
    /** Let this subgroup be G, and the subgroup we're re-indexing by be S.
     *  Since its a subgroup, the 0th element of S is at index 0 in G, the first element of S is at
     *  index |G|/|S|, the second at 2*|G|/|S|, etc.
     *  Thus for an index i that correspond S, the index in G is i*|G|/|S|
    */
   const std::size_t order_s = 1ull << reindex_subgroup_dim;
   const std::size_t order_g_over_s = 1ull << (this->dimension() - reindex_subgroup_dim);
   if (index < order_s) {
       return index * order_g_over_s;
   }
   /** Let i now be the index of this element in G \ S */
   const std::size_t i = index - order_s;
   /** Let x be the number of elements in G \ S, for every element in S. Then x = (|G|/|S| - 1).
    *  At index i in G \ S, the number of elements in S that appear before the index in G to which
    *  i corresponds to, is floor(i / x) + 1.
    *  The +1 is because index 0 of G is S_0, so the position is offset by at least one.
    *  The floor(i / x) term is because after x elements in G \ S, there is one more element from S
    *  that will have appeared in G. */
   const std::size_t x = order_g_over_s - 1;
   return i + (i / x) + 1;
}

template<typename FieldT>
std::size_t multiplicative_subgroup_base<FieldT>::coset_index(
    const std::size_t position, const std::size_t coset_size) const
{
    const size_t num_cosets = this->num_elements() / coset_size;
    return position % num_cosets;
}

template<typename FieldT>
std::size_t multiplicative_subgroup_base<FieldT>::intra_coset_index(
    const std::size_t position, const std::size_t coset_size) const
{
    const size_t num_cosets = this->num_elements() / coset_size;
    return position / num_cosets;
}

template<typename FieldT>
std::size_t multiplicative_subgroup_base<FieldT>::position_by_coset_indices(
        const size_t coset_index, const size_t intra_coset_index, const size_t coset_size) const
{
    const size_t num_cosets = this->num_elements() / coset_size;
    return coset_index + intra_coset_index * num_cosets;
}

template<typename FieldT>
libfqfft::basic_radix2_domain<FieldT> multiplicative_subgroup_base<FieldT>::FFT_eval_domain() const
{
    assert(is_power_of_2(this->order_));
    return *(this->FFT_eval_domain_);
}

template<typename FieldT>
std::vector<FieldT> multiplicative_subgroup<FieldT>::all_elements() const
{
    if (this->elems_->empty()) {
        std::vector<FieldT> elems;
        elems.reserve(this->order_);
        FieldT el = FieldT::one();
        for (size_t i = 0; i < this->order_; i++) {
            elems.emplace_back(el);
            el *= this->g_;
        }
        this->elems_->swap(elems);
    }
    return *this->elems_;
}

template<typename FieldT>
FieldT multiplicative_subgroup<FieldT>::element_by_index(const std::size_t index) const
{
    if (this->elems_->empty()) {
        return libiop::power(this->g_, index);
    } else {
        return this->elems_->operator[](index);
    }
}

template<typename FieldT>
multiplicative_coset<FieldT>::multiplicative_coset(FieldT order)
{
    this->construct_internal(order);

    this->shift_ = FieldT::one();
}

template<typename FieldT>
multiplicative_coset<FieldT>::multiplicative_coset(std::size_t order) :
    multiplicative_coset<FieldT>(FieldT(order))
{
}

template<typename FieldT>
multiplicative_coset<FieldT>::multiplicative_coset(FieldT order, FieldT shift)
{
    this->construct_internal(order);

    this->shift_ = shift;
}

template<typename FieldT>
multiplicative_coset<FieldT>::multiplicative_coset(std::size_t order, FieldT shift) :
    multiplicative_coset<FieldT>(FieldT(order), shift)
{
}

template<typename FieldT>
std::vector<FieldT> multiplicative_coset<FieldT>::all_elements() const
{
    if (this->elems_->empty()) {
        std::vector<FieldT> elems;
        elems.reserve(this->order_);
        FieldT el = this->shift_;
        for (size_t i = 0; i < this->order_; i++) {
            elems.emplace_back(el);
            el *= this->g_;
        }
        this->elems_->swap(elems);
    }
    return *this->elems_;
}

template<typename FieldT>
FieldT multiplicative_coset<FieldT>::element_by_index(const std::size_t index) const
{
    if (this->elems_->empty()) {
        return this->shift_ * libiop::power(this->g_, index);
    } else {
        return this->elems_->operator[](index);
    }
}


template<typename FieldT>
FieldT multiplicative_coset<FieldT>::shift() const
{
    return this->shift_;
}

}
