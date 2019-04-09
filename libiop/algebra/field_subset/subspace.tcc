#include <stdexcept>

#include "libiop/algebra/utils.hpp"

namespace libiop {

/* linear_subspace<FieldT> */

template<typename FieldT>
linear_subspace<FieldT>::linear_subspace(const std::vector<FieldT> &basis) :
    basis_(basis)
{
}

template<typename FieldT>
linear_subspace<FieldT>::linear_subspace(std::vector<FieldT> &&basis) :
    basis_(std::move(basis))
{
}

template<typename FieldT>
std::size_t linear_subspace<FieldT>::dimension() const
{
    return this->basis_.size();
}

template<typename FieldT>
std::size_t linear_subspace<FieldT>::num_elements() const
{
    return (1ull << this->basis_.size());
}

template<typename FieldT>
const std::vector<FieldT>& linear_subspace<FieldT>::basis() const
{
    return this->basis_;
}

template<typename FieldT>
const FieldT& linear_subspace<FieldT>::basis_element(const std::size_t i) const
{
    return this->basis_[i];
}

template<typename FieldT>
std::vector<FieldT> linear_subspace<FieldT>::all_elements() const
{
    return all_subset_sums<FieldT>(this->basis_);
}

template<typename FieldT>
FieldT linear_subspace<FieldT>::element_by_index(const std::size_t index) const
{
    if (index >= this->num_elements())
    {
        throw std::invalid_argument("element index out of bounds");
    }

    FieldT result = FieldT(0);
    for (std::size_t i = 0; i < this->basis_.size(); ++i)
    {
        if (index & (1ull<<i))
        {
            result += this->basis_[i];
        }
    }

    return result;
}

template<typename FieldT>
std::size_t linear_subspace<FieldT>::coset_index(const std::size_t position, const std::size_t coset_size) const
{
    return position / coset_size;
}

/** Assumes that the coset uses the first few basis vectors as this subspace. */
template<typename FieldT>
std::size_t linear_subspace<FieldT>::intra_coset_index(const std::size_t position, const std::size_t coset_size) const
{
    return position % coset_size;
}

template<typename FieldT>
std::size_t linear_subspace<FieldT>::position_by_coset_indices(
        const size_t coset_index, const size_t intra_coset_index, const size_t coset_size) const
{
    return coset_index * coset_size + intra_coset_index;
}

template<typename FieldT>
linear_subspace<FieldT> linear_subspace<FieldT>::standard_basis(const std::size_t dimension)
{
    //assert(dimension <= extension_degree_helper(FieldT::zero()));

    std::vector<FieldT> basis_elements;
    basis_elements.reserve(dimension);

    for (std::size_t i = 0; i < dimension; ++i)
    {
        basis_elements.emplace_back(FieldT(1ull<<i));
    }

    return linear_subspace<FieldT>(std::move(basis_elements));
}

template<typename FieldT>
linear_subspace<FieldT> linear_subspace<FieldT>::random_linear_subspace(const std::size_t dimension)
{
    const std::vector<FieldT> basis_elements = random_FieldT_vector<FieldT>(dimension);
    return linear_subspace<FieldT>(std::move(basis_elements));
}

/* affine_subspace<FieldT> */

template<typename FieldT>
affine_subspace<FieldT>::affine_subspace(const std::vector<FieldT> &basis,
                                         const FieldT &offset) :
    linear_subspace<FieldT>(basis), offset_(offset)
{
}

template<typename FieldT>
affine_subspace<FieldT>::affine_subspace(
    const linear_subspace<FieldT> &base_space,
    const FieldT &offset) :
    linear_subspace<FieldT>(base_space),
    offset_(offset)
{
}

template<typename FieldT>
affine_subspace<FieldT>::affine_subspace(
    linear_subspace<FieldT> &&base_space,
    const FieldT &offset) :
    linear_subspace<FieldT>(std::move(base_space)),
    offset_(offset)
{
}

template<typename FieldT>
const FieldT& affine_subspace<FieldT>::offset() const
{
    return this->offset_;
}

template<typename FieldT>
std::vector<FieldT> affine_subspace<FieldT>::all_elements() const
{
    return all_subset_sums<FieldT>(this->basis_, this->offset_);
}

template<typename FieldT>
FieldT affine_subspace<FieldT>::element_by_index(const std::size_t index) const
{
    return (this->offset_ + (linear_subspace<FieldT>::element_by_index(index)));
}

template<typename FieldT>
linear_subspace<FieldT> standard_basis(const std::size_t dimension)
{
    std::vector<FieldT> basis;
    basis.reserve(dimension);
    for (size_t i = 0; i < dimension; ++i)
    {
        basis.emplace_back(FieldT(1ull<<i));
    }

    return linear_subspace<FieldT>(std::move(basis));
}

template<typename FieldT>
affine_subspace<FieldT> affine_subspace<FieldT>::shifted_standard_basis(
    const std::size_t dimension,
    const FieldT& offset)
{
    const linear_subspace<FieldT> basis =
        linear_subspace<FieldT>::standard_basis(dimension);
    return affine_subspace<FieldT>(std::move(basis), offset);
}

template<typename FieldT>
affine_subspace<FieldT> affine_subspace<FieldT>::random_affine_subspace(const std::size_t dimension)
{
    const linear_subspace<FieldT> basis =
        linear_subspace<FieldT>::standard_basis(dimension);
    const FieldT offset = FieldT::random_element();
    return affine_subspace<FieldT>(std::move(basis), offset);
}

} // namespace libiop
