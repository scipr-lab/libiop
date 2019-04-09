namespace libiop {

template<typename FieldT>
FieldT power(const FieldT &base, const std::size_t exponent)
{
    FieldT result = FieldT::one();

    bool found_one = false;

    for (long i = 8 * sizeof(exponent) - 1; i >= 0; --i)
    {
        if (found_one)
        {
            result = result.squared();
        }

        if (exponent & (1ull << i))
        {
            found_one = true;
            result *= base;
        }
    }

    return result;
}

template<typename FieldT>
std::vector<FieldT> subspace_element_powers(const affine_subspace<FieldT> &S,
                                            const std::size_t exponent)
{
    std::vector<FieldT> result(S.num_elements(), FieldT::one());

    for (std::size_t i = 0; i < 8 * sizeof(exponent); ++i)
    {
        if (!(exponent & (1ull << i)))
        {
            continue;
        }

        /* exponent has 2^i set, so we multiply in all subset sums of
           linearized polynomial 2^i. doing this for affine subspaces
           requires a trick; see comment at
           linearized_polynomial::evaluations_over_subspace */
        const std::size_t two_to_i = (1ull << i);
        std::vector<FieldT> basis_powers(S.basis());
        for (auto &el : basis_powers)
        {
            el = libiop::power(el, two_to_i);
        }

        const FieldT offset_power = libiop::power(S.offset(), two_to_i);
        const std::vector<FieldT> subspace_to_two_to_i =
            all_subset_sums<FieldT>(basis_powers, offset_power);

        for (std::size_t j = 0; j < S.num_elements(); ++j)
        {
            result[j] *= subspace_to_two_to_i[j];
        }
    }

    return result;
}

template<typename FieldT>
std::vector<FieldT> coset_element_powers(const multiplicative_coset<FieldT> &S,
                                         const std::size_t exponent)
{
    std::vector<FieldT> result;
    result.reserve(S.num_elements());

    const FieldT g_to_exp = libiop::power(S.generator(), exponent);
    FieldT current_term = libiop::power(S.shift(), exponent);

    for (std::size_t i = 0; i < S.num_elements(); i++) {
        result.emplace_back(current_term);
        current_term *= g_to_exp;
    }

    return result;
}

template<typename FieldT>
std::vector<FieldT> subset_element_powers(const field_subset<FieldT> &S,
                                                 const std::size_t exponent)
{
    if (S.type() == affine_subspace_type)
    {
        return subspace_element_powers(S.subspace(), exponent);
    }
    else
    {
        return coset_element_powers(S.coset(), exponent);
    }
}


} // libiop
