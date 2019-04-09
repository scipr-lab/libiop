#include <cstddef>

#include <libfqfft/evaluation_domain/domains/basic_radix2_domain.hpp>

#include "libiop/algebra/fields/utils.hpp"
#include "libiop/algebra/utils.hpp"
#include "libiop/common/profiling.hpp"

namespace libiop {

/* Performs naive computation of the polynomial evaluation
   problem. Mostly useful for testing. */
template<typename FieldT>
std::vector<FieldT> naive_FFT(const std::vector<FieldT> &poly_coeffs,
                              const field_subset<FieldT> &domain)
{
    const std::size_t n = poly_coeffs.size();
    const std::size_t m = domain.dimension();

    const std::vector<FieldT> evalpoints = domain.all_elements();
    std::vector<FieldT> result;
    result.reserve(n);

    for (const FieldT &p : evalpoints)
    {
        FieldT v(0);
        for (size_t i = n; i--; )
        {
            v *= p;
            v += poly_coeffs[i];
        }

        result.emplace_back(v);
    }

    return result;
}

template<typename FieldT>
std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                 const affine_subspace<FieldT> &domain)
{
    std::vector<FieldT> S(poly_coeffs);
    S.resize(domain.num_elements(), FieldT::zero());

    const size_t n = S.size();
    const size_t m = domain.dimension();
    assert(n == (1ull<<m));

    std::vector<FieldT> recursed_betas((m+1)*m/2, FieldT(0));
    std::vector<FieldT> recursed_shifts(m, FieldT(0));
    size_t recursed_betas_ptr = 0;

    std::vector<FieldT> betas2(domain.basis());
    FieldT shift2 = domain.offset();
    for (size_t j = 0; j < m; ++j)
    {
        FieldT beta = betas2[m-1-j];
        FieldT betai(1);

        /* twist by beta. TODO: this can often be elided by a careful choice of betas */
        for (size_t ofs = 0; ofs < n; ofs += (1ull<<j))
        {
            for (size_t p = 0; p < (1ull<<j); ++p)
            {
                S[ofs + p] *= betai;
            }

            betai *= beta;
        }

        /* perform radix conversion */
        for (size_t stride = n/4; stride >= (1ul << j); stride >>= 1)
        {
            for (size_t ofs = 0; ofs < n; ofs += stride*4)
            {
                for (size_t i = 0; i < stride; ++i)
                {
                    S[ofs+2*stride+i] += S[ofs+3*stride+i];
                    S[ofs+1*stride+i] += S[ofs+2*stride+i];
                }
            }
        }

        /* compute deltas used in the reverse process */
        FieldT betainv = beta.inverse();
        for (size_t i = 0; i < m-1-j; ++i)
        {
            FieldT newbeta = betas2[i] * betainv;
            recursed_betas[recursed_betas_ptr++] = newbeta;
            betas2[i] = newbeta.squared() - newbeta;
        }

        FieldT newshift = shift2 * betainv;
        recursed_shifts[j] = newshift;
        shift2 = newshift.squared() - newshift;
    }

    bitreverse_vector<FieldT>(S);

    /* unwind the recursion */
    for (size_t j = 0; j < m; ++j)
    {
        recursed_betas_ptr -= j;
        /* note that this devolves to empty range for the first loop iteration */
        std::vector<FieldT> popped_betas = std::vector<FieldT>(recursed_betas.begin()+recursed_betas_ptr,
                                                               recursed_betas.begin()+recursed_betas_ptr+j);
        const FieldT popped_shift = recursed_shifts[m-1-j];
        std::vector<FieldT> sums = all_subset_sums<FieldT>(popped_betas, popped_shift);

        size_t stride = 1ull<<j;
        for (size_t ofs = 0; ofs < n; ofs += 2*stride)
        {
            for (size_t i = 0; i < stride; ++i)
            {
                S[ofs+i] += S[ofs+stride+i] * sums[i];
                S[ofs+stride+i] += S[ofs+i];
            }
        }
    }
    assert(recursed_betas_ptr == 0);

    return S;
}

template<typename FieldT>
std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &evals,
                                  const affine_subspace<FieldT> &domain)
{
    const size_t n = evals.size();
    const size_t m = domain.dimension();
    assert(n == (1ull<<m));

    std::vector<FieldT> S(evals);
    std::vector<FieldT> recursed_twists(m, FieldT(0));

    std::vector<FieldT> betas2(domain.basis());
    FieldT shift2 = domain.offset();
    for (size_t j = 0; j < m; ++j)
    {
        const FieldT beta = betas2[m-1-j];
        const FieldT betainv = beta.inverse();
        recursed_twists[j] = betainv;

        std::vector<FieldT> newbetas(m-1-j, FieldT(0));

        for (size_t i = 0; i < m-1-j; ++i)
        {
            FieldT newbeta = betas2[i] * betainv;
            newbetas[i] = newbeta;
            betas2[i] = newbeta.squared() - newbeta;
        }

        FieldT newshift = shift2 * betainv;
        shift2 = newshift.squared() - newshift;

        const std::vector<FieldT> sums = all_subset_sums<FieldT>(newbetas, newshift);

        const size_t half = 1ull<<(m-1-j);
        for (size_t ofs = 0; ofs < n; ofs += 2*half)
        {
            for (size_t p = 0; p < half; ++p)
            {
                S[ofs + half + p] += S[ofs + p];
                S[ofs + p] += S[ofs + half + p] * sums[p];
            }
        }
    }

    bitreverse_vector<FieldT>(S);

    for (size_t j = 0; j < m; ++j)
    {
        size_t N = 4ull<<(m-1-j);
        /* perform radix combinations */
        while (N <= n)
        {
            const size_t quarter = N/4;
            for (size_t ofs = 0; ofs < n; ofs += N)
            {
                for (size_t i = 0; i < quarter; ++i)
                {
                    S[ofs+1*quarter+i] += S[ofs+2*quarter+i];
                    S[ofs+2*quarter+i] += S[ofs+3*quarter+i];
                }
            }
            N *= 2;
        }

        /* twist by \beta^{-1} */
        const FieldT betainv = recursed_twists[m-1-j];
        FieldT betainvi(1);
        for (size_t ofs = 0; ofs < n; ofs += (1ull<<(m-1-j)))
        {
            for (size_t p = 0; p < (1ull<<(m-1-j)); ++p)
            {
                S[ofs + p] *= betainvi;
            }
            betainvi *= betainv;
        }
    }

    return S;
}

template<typename FieldT>
std::vector<FieldT> additive_FFT_wrapper(const std::vector<FieldT> &v,
                                         const affine_subspace<FieldT> &H)
{
    enter_block("Call to additive_FFT_wrapper");
    print_indent(); printf("* Vector size: %zu\n", v.size());
    print_indent(); printf("* Subspace size: %zu\n", H.num_elements());
    const std::vector<FieldT> result = additive_FFT(v, H);
    leave_block("Call to additive_FFT_wrapper");
    return result;
}

template<typename FieldT>
std::vector<FieldT> additive_IFFT_wrapper(const std::vector<FieldT> &v,
                                          const affine_subspace<FieldT> &H)
{
    enter_block("Call to additive_IFFT_wrapper");
    print_indent(); printf("* Vector size: %zu\n", v.size());
    print_indent(); printf("* Subspace size: %zu\n", H.num_elements());
    const std::vector<FieldT> result = additive_IFFT(v, H);
    leave_block("Call to additive_IFFT_wrapper");
    return result;
}

template<typename FieldT>
std::vector<FieldT> multiplicative_FFT_internal(
    const std::vector<typename libiop::enable_if<is_multiplicative<FieldT>::value, FieldT>::type> &poly_coeffs,
    const multiplicative_subgroup_base<FieldT> &domain, const FieldT shift)
{
    assert(poly_coeffs.size() <= domain.num_elements());

    libfqfft::basic_radix2_domain<FieldT> eval_domain = domain.FFT_eval_domain();

    std::vector<FieldT> vec = poly_coeffs;
    vec.resize(domain.num_elements(), FieldT::zero());
    // Handle separately, as cosetFFT requires more multiplications
    if (shift == FieldT::one()) {
        eval_domain.FFT(vec);
    } else {
        eval_domain.cosetFFT(vec, shift);
    }

    return vec;
}

template<typename FieldT>
std::vector<FieldT> multiplicative_FFT_internal(
    const std::vector<typename libiop::enable_if<is_additive<FieldT>::value, FieldT>::type> &poly_coeffs,
    const multiplicative_subgroup_base<FieldT> &domain, const FieldT shift)
{
    throw std::invalid_argument("attempting to perform multiplicative IFFT with non-multiplicative field type");
}

template<typename FieldT>
std::vector<FieldT> multiplicative_FFT(const std::vector<FieldT> &poly_coeffs,
                                       const multiplicative_coset<FieldT> &domain)
{
    return multiplicative_FFT_internal(poly_coeffs, domain, domain.shift());
}

template<typename FieldT>
std::vector<FieldT> multiplicative_IFFT_internal(
    const std::vector<typename libiop::enable_if<is_multiplicative<FieldT>::value, FieldT>::type> &evals,
    const multiplicative_subgroup_base<FieldT> &domain, const FieldT shift)
{
    assert(domain.num_elements() == evals.size());

    libfqfft::basic_radix2_domain<FieldT> eval_domain = domain.FFT_eval_domain();

    std::vector<FieldT> vec = evals;
    // Handle separately, as icosetFFT requires more multiplications
    if (shift == FieldT::one()) {
        eval_domain.iFFT(vec);
    } else {
        eval_domain.icosetFFT(vec, shift);
    }

    return vec;
}

template<typename FieldT>
std::vector<FieldT> multiplicative_IFFT_internal(
    const std::vector<typename libiop::enable_if<is_additive<FieldT>::value, FieldT>::type> &evals,
    const multiplicative_subgroup_base<FieldT> &domain, const FieldT shift)
{
    throw std::invalid_argument("attempting to perform multiplicative IFFT with non-multiplicative field type");
}

template<typename FieldT>
std::vector<FieldT> multiplicative_IFFT(const std::vector<FieldT> &evals,
                                        const multiplicative_coset<FieldT> &domain)
{
    return multiplicative_IFFT_internal(evals, domain, domain.shift());
}

template<typename FieldT>
std::vector<FieldT> multiplicative_FFT_wrapper(const std::vector<FieldT> &v,
                                               const multiplicative_coset<FieldT> &H)
{
    enter_block("Call to multiplicative_FFT_wrapper");
    print_indent(); printf("* Vector size: %zu\n", v.size());
    print_indent(); printf("* Subgroup size: %zu\n", H.num_elements());
    const std::vector<FieldT> result = multiplicative_FFT(v, H);
    leave_block("Call to multiplicative_FFT_wrapper");
    return result;
}

template<typename FieldT>
std::vector<FieldT> multiplicative_IFFT_wrapper(const std::vector<FieldT> &v,
                                                const multiplicative_coset<FieldT> &H)
{
    enter_block("Call to multiplicative_IFFT_wrapper");
    print_indent(); printf("* Vector size: %zu\n", v.size());
    print_indent(); printf("* Coset size: %zu\n", H.num_elements());
    const std::vector<FieldT> result = multiplicative_IFFT(v, H);
    leave_block("Call to multiplicative_IFFT_wrapper");
    return result;
}

template<typename FieldT>
std::vector<FieldT> FFT_over_field_subset(const std::vector<typename libiop::enable_if<is_multiplicative<FieldT>::value, FieldT>::type> coeffs,
                                          field_subset<FieldT> domain)
{
    return multiplicative_FFT_wrapper<FieldT>(coeffs, domain.coset());
}

template<typename FieldT>
std::vector<FieldT> FFT_over_field_subset(const std::vector<typename libiop::enable_if<is_additive<FieldT>::value, FieldT>::type> coeffs,
                                          field_subset<FieldT> domain)
{
    return additive_FFT_wrapper<FieldT>(coeffs, domain.subspace());
}

template<typename FieldT>
std::vector<FieldT> IFFT_over_field_subset(const std::vector<typename libiop::enable_if<is_multiplicative<FieldT>::value, FieldT>::type> evals,
                                           field_subset<FieldT> domain)
{
    return multiplicative_IFFT_wrapper<FieldT>(evals, domain.coset());
}

template<typename FieldT>
std::vector<FieldT> IFFT_over_field_subset(const std::vector<typename libiop::enable_if<is_additive<FieldT>::value, FieldT>::type> evals,
                                           field_subset<FieldT> domain)
{
    return additive_IFFT_wrapper<FieldT>(evals, domain.subspace());
}

} // namespace libiop
