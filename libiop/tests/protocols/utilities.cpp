#include <cstdint>
#include <memory>
#include <stdexcept>

#include <gtest/gtest.h>

#include "libiop/algebra/fields/gf64.hpp"
#include "libiop/algebra/field_subset/field_subset.hpp"
#include "libiop/algebra/polynomials/polynomial.hpp"
#include "libiop/iop/iop.hpp"

#include <libff/algebra/curves/edwards/edwards_pp.hpp>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>

namespace libiop {

/** This method is deprecated in favor of the latter method,
 * It is kept for debugging purposes. */
template<typename FieldT>
FieldT sum_over_field_subset(const polynomial<FieldT> &P,
                             const field_subset<FieldT> &S) {
    FieldT sum = FieldT::zero();
    for (auto &el : S.all_elements())
    {
        sum += P.evaluation_at_point(el);
    }
    return sum;
}

template<typename FieldT>
FieldT sum_over_default_field_subset(const polynomial<FieldT> &P,
                                     const field_subset<FieldT> &S) {
    // This assumes that S was created over the default basis.
    // It then creates an extended domain, and does an FFT to convert P to evaluations
    // in that extended domain.
    const std::size_t dim = std::max(log2(P.num_terms()), S.dimension());
    const field_subset<FieldT> extended_subset(1ull << dim);
    const std::vector<FieldT> evals = FFT_over_field_subset(P.coefficients(), extended_subset);
    FieldT sum = FieldT::zero();
    for (std::size_t i = 0; i < S.num_elements(); i++)
    {
        sum += evals[extended_subset.reindex_by_subset(S.dimension(),i)];
    }
    return sum;
}

template<typename FieldT>
std::size_t degree_bound_from_evals(const std::vector<FieldT> &evals,
                                    const field_subset<FieldT> &domain)
{
    const polynomial<FieldT> poly_from_evals(IFFT_over_field_subset<FieldT>(evals, domain));
    const std::size_t minimal_num_terms = poly_from_evals.minimal_num_terms();

    return minimal_num_terms;
}

template<typename FieldT>
void test_oracle_consistency(iop_protocol<FieldT> &IOP,
                             const oracle_handle_ptr &handle,
                             const std::vector<FieldT> &oracle_evals,
                             const field_subset<FieldT> &codeword_domain) {
    for (std::size_t i = 0; i < 10; i++) {
        std::size_t evaluation_index = std::rand() % codeword_domain.num_elements();
        const FieldT point_eval = IOP.get_oracle_evaluation_at_point(handle, evaluation_index, false);
        EXPECT_TRUE(point_eval == oracle_evals[evaluation_index]) << "evaluation at point was inconsistent at index " << 
            evaluation_index;
    }
}

template<typename FieldT>
void test_oracles_degree_and_consistency(iop_protocol<FieldT> &IOP,
                                        const std::vector<oracle_handle_ptr> &handles,
                                        const field_subset<FieldT> &codeword_domain,
                                        const bool exp_pass)
{
    bool passed = true;
    for (std::size_t i = 0; i < handles.size(); i++)
    {
        const oracle_handle_ptr handle = handles[i];
        const std::vector<FieldT> evals = IOP.get_oracle_evaluations(handle);
        const std::size_t expected_degree = IOP.get_oracle_degree(handle);
        const std::size_t oracle_degree = degree_bound_from_evals(evals, codeword_domain);
        if (exp_pass) {
            // TODO: Make handles have a concept of human readable name.
            ASSERT_EQ(expected_degree, oracle_degree) << "handle " << i << " did not have correct degree.";
        } else {
            passed = passed & (oracle_degree <= expected_degree);
        }
        test_oracle_consistency(IOP, handle, evals, codeword_domain);
    }
    if (!exp_pass) {
        ASSERT_FALSE(passed) << "All of the oracle handles had passing degrees.";
    }
}

template<typename FieldT>
std::vector<FieldT> make_codeword_zk(std::vector<FieldT> codeword,
                                     std::size_t query_bound,
                                     field_subset<FieldT> systematic_domain,
                                     field_subset<FieldT> codeword_domain) {
    const vanishing_polynomial<FieldT> constraint_vp(systematic_domain);
    const std::vector<FieldT> constraint_vp_over_codeword_domain =
        constraint_vp.evaluations_over_field_subset(codeword_domain);
    const polynomial<FieldT> R = polynomial<FieldT>::random_polynomial(query_bound);
    const std::vector<FieldT> R_over_codeword_domain = FFT_over_field_subset(
        R.coefficients(), codeword_domain);
    std::vector<FieldT> zk_codeword;
    for (std::size_t i = 0; i < codeword_domain.num_elements(); ++i)
    {
        zk_codeword.emplace_back(codeword[i] + 
            R_over_codeword_domain[i] * constraint_vp_over_codeword_domain[i]);
    }
    return zk_codeword;
}
}