/**@file
 *****************************************************************************
 Classes for vanishing polynomials for field subsets.
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_ALGEBRA_POLYNOMIALS_VANISHING_POLYNOMIAL_HPP_
#define LIBIOP_ALGEBRA_POLYNOMIALS_VANISHING_POLYNOMIAL_HPP_

#include <cstddef>
#include <vector>

#include "libiop/algebra/exponentiation.hpp"
#include "libiop/algebra/field_subset/field_subset.hpp"
#include "libiop/algebra/field_subset/subspace.hpp"
#include "libiop/algebra/field_subset/subgroup.hpp"
#include "libiop/algebra/polynomials/polynomial.hpp"
#include "libiop/algebra/polynomials/linearized_polynomial.hpp"

namespace libiop {

template<typename FieldT>
class vanishing_polynomial {
private:
    const field_subset_type type_;
    std::size_t vp_degree_;
    // subspace type
    linearized_polynomial<FieldT> linearized_polynomial_;
    // multiplicative coset type
    FieldT vp_offset_; /* offset^|H| for cosets, 1 for subgroups */

public:
    explicit vanishing_polynomial(const field_subset<FieldT> &S);
    explicit vanishing_polynomial(const affine_subspace<FieldT> &S);
    explicit vanishing_polynomial(const multiplicative_coset<FieldT> &S);

    FieldT evaluation_at_point(const FieldT &evalpoint) const;
    std::vector<FieldT> evaluations_over_field_subset(const field_subset<FieldT> &S) const;
    std::vector<FieldT> evaluations_over_subspace(const affine_subspace<FieldT> &S) const;
    std::vector<FieldT> evaluations_over_coset(const multiplicative_coset<FieldT> &S) const;

    polynomial<FieldT> operator*(const polynomial<FieldT> &p) const;

    std::size_t degree() const;
    FieldT constant_coefficient() const;
    linearized_polynomial<FieldT> get_linearized_polynomial() const;
    polynomial<FieldT> get_polynomial() const;
    field_subset_type type() const;
};

template<typename FieldT>
std::pair<polynomial<FieldT>,
          polynomial<FieldT> >
polynomial_over_vanishing_polynomial(const polynomial<FieldT> &f,
                                     const vanishing_polynomial<FieldT> &Z);

template<typename FieldT>
linearized_polynomial<FieldT> vanishing_polynomial_from_subspace(const affine_subspace<FieldT> &S);

} // namespace libiop

#include "libiop/algebra/polynomials/vanishing_polynomial.tcc"

#endif // LIBIOP_ALGEBRA_POLYNOMIALS_VANISHING_POLYNOMIAL_HPP_
