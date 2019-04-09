#include <algorithm>

#include "libiop/algebra/exponentiation.hpp"
#include "libiop/algebra/field_subset/field_subset.hpp"
#include "libiop/algebra/polynomials/polynomial.hpp"
#include "libiop/algebra/polynomials/linearized_polynomial.hpp"
#include "libiop/algebra/fft.hpp"

namespace libiop {

/* vanishing_polynomial<FieldT>(field_subset<FieldT>) */
template<typename FieldT>
vanishing_polynomial<FieldT>::vanishing_polynomial(const field_subset<FieldT> &S) :
    type_(S.type())
{
    this->vp_degree_ = S.num_elements();
    if (this->type_ == affine_subspace_type) {
        this->linearized_polynomial_ = vanishing_polynomial_from_subspace(S.subspace());
    } else if (this->type_ == multiplicative_coset_type) {
        this->vp_offset_ = libiop::power(S.coset().shift(), this->vp_degree_);
    } else {
        throw std::invalid_argument("field_subset type unsupported.");
    }
}

template<typename FieldT>
vanishing_polynomial<FieldT>::vanishing_polynomial(const affine_subspace<FieldT> &S) :
    type_(affine_subspace_type)
{
    this->vp_degree_ = S.num_elements();
    this->linearized_polynomial_ = vanishing_polynomial_from_subspace(S);
}

template<typename FieldT>
vanishing_polynomial<FieldT>::vanishing_polynomial(const multiplicative_coset<FieldT> &S) :
    type_(multiplicative_coset_type)
{
    this->vp_degree_ = S.num_elements();
    this->vp_offset_ = libiop::power(S.coset().shift(), this->vp_degree_);
}

template<typename FieldT>
FieldT vanishing_polynomial<FieldT>::evaluation_at_point(const FieldT &evalpoint) const {
    if (this->type_ == affine_subspace_type) {
        return this->linearized_polynomial_.evaluation_at_point(evalpoint);
    } else if (this->type_ == multiplicative_coset_type) {
        return libiop::power(evalpoint, this->vp_degree_) - this->vp_offset_;
    }
    return FieldT::zero();
}

template<typename FieldT>
std::vector<FieldT> vanishing_polynomial<FieldT>::evaluations_over_field_subset(const field_subset<FieldT> &S) const {
    if (S.type() == affine_subspace_type) {
        return this->evaluations_over_subspace(S.subspace());
    } else if (S.type() == multiplicative_coset_type) {
        return this->evaluations_over_coset(S.coset());
    }
    throw std::invalid_argument("field_subset type unsupported");
}

template<typename FieldT>
std::vector<FieldT> vanishing_polynomial<FieldT>::evaluations_over_subspace(const affine_subspace<FieldT> &S) const {
    if (this->type_ != affine_subspace_type) {
        throw std::invalid_argument("evaluations_over_subspace can only be used on subspace vanishing polynomials.");
    }
    return this->linearized_polynomial_.evaluations_over_subspace(S);
}

template<typename FieldT>
std::vector<FieldT> vanishing_polynomial<FieldT>::evaluations_over_coset(const multiplicative_coset<FieldT> &S) const {
    if (this->type_ != multiplicative_coset_type) {
        throw std::invalid_argument("evaluations_over_coset can only be used on multiplicative_coset vanishing polynomials.");
    }
    // P is of the form X^|G| - vp_offset_, for |G| some power of 2, and |S| is also a power of 2
    const std::size_t order_s = S.num_elements();
    const std::size_t order_g = this->vp_degree_;
    // Evaluations are the same as vanishing polynomial over subgroup case, but multiplied by shift^|G|
    const FieldT shift_to_order_g = libiop::power(S.shift(), order_g);
    std::vector<FieldT> evals;
    if (order_s <= order_g) {
        evals.resize(order_s, shift_to_order_g - this->vp_offset_);
        return evals;
    }
    // |G| divides |S|, therefore X^|G| is a homomorphism from S to a subgroup of order |S| / |G|.
    // Since P = x^|G| - 1, and since 1 is independent of X, it follows that there are
    // only |S| / |G| distinct evaluations, and these repeat.
    evals.reserve(order_s);
    const std::size_t order_s_over_g = order_s / order_g;
    const FieldT generator_to_order_g = libiop::power(S.generator(), order_g);
    FieldT cur = FieldT::one();
    for (std::size_t i = 0; i < order_s_over_g; i++) {
        evals.emplace_back(cur * shift_to_order_g - this->vp_offset_);
        cur = cur * generator_to_order_g;
    }
    // Place these |S| / |G| distinct evaluations in the remaining locations.
    for (std::size_t i = 1; i < order_g; i++) {
        for (std::size_t j = 0; j < order_s_over_g; j++) {
            evals.emplace_back(evals[j]);
        }
    }
    return evals;
}

template<typename FieldT>
std::size_t vanishing_polynomial<FieldT>::degree() const {
    return this->vp_degree_;
}

template<typename FieldT>
FieldT vanishing_polynomial<FieldT>::constant_coefficient() const {
    if (this->type() == affine_subspace_type) {
        return this->linearized_polynomial_.constant_coefficient();
    }
    // subgroup / coset type
    return -this->vp_offset_;
}

template<typename FieldT>
polynomial<FieldT> vanishing_polynomial<FieldT>::operator*(const polynomial<FieldT> &p) const
{
    if (this->type_ == affine_subspace_type) {
        return this->linearized_polynomial_ * p;
    }
    // in the multiplicative case just shift p, and subtract by p * this->vp_offset_
    std::vector<FieldT> result(p.degree() + this->vp_degree_ + 1, FieldT(0));
    const std::vector<FieldT> p_coeff = p.coefficients();
    add_scalar_multiple_at_offset(result, p_coeff, FieldT(1), this->vp_degree_);
    add_scalar_multiple_at_offset(result, p_coeff, FieldT(0) - this->vp_offset_, 0);
    return polynomial<FieldT>(std::move(result));
}

template<typename FieldT>
linearized_polynomial<FieldT> vanishing_polynomial<FieldT>::get_linearized_polynomial() const {
    if (this->type_ == multiplicative_coset_type) {
        std::vector<FieldT> coefficients(libiop::log2(this->vp_degree_) + 1, FieldT::zero());
        coefficients.emplace_back(FieldT::one());
        coefficients[0] = -this->vp_offset_;
        linearized_polynomial<FieldT> result(std::move(coefficients));
        return result;
    }
    return this->linearized_polynomial_;
}

template<typename FieldT>
polynomial<FieldT> vanishing_polynomial<FieldT>::get_polynomial() const {
    if (this->type_ == affine_subspace_type) {
        return this->linearized_polynomial_.expand_as_polynomial();
    }
    polynomial<FieldT> poly;
    poly.set_degree(this->vp_degree_);
    poly[0] = FieldT(-1);
    poly[this->vp_degree_] = FieldT(1);
    return poly;
}

template<typename FieldT>
field_subset_type vanishing_polynomial<FieldT>::type() const {
    return this->type_;
}

/* Returns P / Z as a pair.
   The first element in the pair is the quotient, the latter is the remainder. */
template<typename FieldT>
std::pair<polynomial<FieldT>,
          polynomial<FieldT> >
polynomial_over_multiplicative_vanishing_polynomial(const polynomial<FieldT> &P,
                                              const FieldT vp_offset,
                                              const size_t vp_degree)
{
    /* inverse of the leading term */
    const FieldT linv = FieldT::one().inverse();

    if (P.degree() < vp_degree)
    {
        polynomial<FieldT> quotient;
        std::vector<FieldT> coeff_copy = P.coefficients();
        polynomial<FieldT> remainder(std::move(coeff_copy));
        remainder.set_degree(vp_degree-1);

        return std::make_pair(std::move(quotient), std::move(remainder));
    }

    std::vector<FieldT> quotient(P.coefficients().begin() + vp_degree,
                          P.coefficients().end());
    std::vector<FieldT> remainder(P.coefficients().begin(),
                          P.coefficients().begin() + vp_degree);

    FieldT Z_0 = vp_offset;
    for (std::size_t i = quotient.size(); i--; )
    {
        // Z only has 2 terms, the leading term and the constant term.
        const FieldT twist = quotient[i] * linv;
        // correct the quotient for the leading term
        quotient[i] = twist;

        /* Handle the remainder by subtracting twist*Z[0] * y^i,
         * thus clearing the i-th term of P */
        if (i < vp_degree)
        {
            remainder[i] -= twist * Z_0;
        }
        else
        {
            quotient[i - vp_degree] -= twist * Z_0;
        }
    }

    return std::make_pair(std::move(polynomial<FieldT>(std::move(quotient))),
                          std::move(polynomial<FieldT>(std::move(remainder))));
}

template<typename FieldT>
std::pair<polynomial<FieldT>,
    polynomial<FieldT> >
polynomial_over_vanishing_polynomial(const polynomial<FieldT> &f,
                                     const vanishing_polynomial<FieldT> &Z)
{
    if (Z.type() == affine_subspace_type) {
        return polynomial_over_linearized_polynomial(f, Z.get_linearized_polynomial());
    } else {
        return polynomial_over_multiplicative_vanishing_polynomial(f, Z.constant_coefficient(), Z.degree());
    }
};

template<typename FieldT>
linearized_polynomial<FieldT> vanishing_polynomial_from_subspace(const affine_subspace<FieldT> &S)
{
    /* Vanishing polynomial for empty subspace is just y */
    linearized_polynomial<FieldT> poly({ FieldT(0), FieldT(1) });

    for (const FieldT &c : S.basis())
    {
        /* Note that since we are in a binary field,
          Z_{<S_1,...,S_k>}(y)
            = Z_{<S_1,...,S_{k-1}>}(y) * Z_{<S_1,...,S_{k-1}>}(y+S_k)
          By linearity, this is equivalent to:
            = Z_{<S_1,...,S_{k-1}>}(y)^2 +
                (Z_{<S_1,...,S_{k-1}>}(y)*Z_{<S_1,...,S_{k-1}>}(S_k)) */
        const FieldT poly_c = poly.evaluation_at_point(c);
        poly = poly.squared() + (poly * poly_c);
    }

    const FieldT poly_shift = poly.evaluation_at_point(S.offset());
    poly[0] += poly_shift;

    return poly;
}

} // namespace libiop
