#include "libiop/algebra/lagrange.hpp"
#include "libiop/algebra/utils.hpp"

namespace libiop {

template<typename FieldT>
rowcheck_ABC_virtual_oracle<FieldT>::rowcheck_ABC_virtual_oracle(
    const field_subset<FieldT> &codeword_domain,
    const field_subset<FieldT> &constraint_domain) :
    codeword_domain_(codeword_domain),
    constraint_domain_(constraint_domain),
    Z_(constraint_domain)
{
}

/** Takes as input oracles Az, Bz, Cz. */
template<typename FieldT>
std::vector<FieldT> rowcheck_ABC_virtual_oracle<FieldT>::evaluated_contents(
    const std::vector<std::vector<FieldT> > &constituent_oracle_evaluations) const
{
    enter_block("rowcheck evaluated contents");
    if (constituent_oracle_evaluations.size() != 3)
    {
        throw std::invalid_argument("rowcheck_ABC has three constituent oracles.");
    }

    const std::vector<FieldT> &Az = constituent_oracle_evaluations[0];
    const std::vector<FieldT> &Bz = constituent_oracle_evaluations[1];
    const std::vector<FieldT> &Cz = constituent_oracle_evaluations[2];
    const std::vector<FieldT> &Z_inv =
        batch_inverse(this->Z_.evaluations_over_field_subset(this->codeword_domain_));

    std::vector<FieldT> result;
    result.reserve(this->codeword_domain_.num_elements());
    for (std::size_t i = 0; i < this->codeword_domain_.num_elements(); ++i)
    {
        result.emplace_back(Z_inv[i] * (Az[i] * Bz[i] - Cz[i]));
    }

    leave_block("rowcheck evaluated contents");
    return result;
}

template<typename FieldT>
FieldT rowcheck_ABC_virtual_oracle<FieldT>::evaluation_at_point(
    const std::size_t evaluation_position,
    const FieldT evaluation_point,
    const std::vector<FieldT> &constituent_oracle_evaluations) const
{
    UNUSED(evaluation_position);

    if (constituent_oracle_evaluations.size() != 3)
    {
        throw std::invalid_argument("rowcheck_ABC has three constituent oracles.");
    }

    const FieldT &A_X = constituent_oracle_evaluations[0];
    const FieldT &B_X = constituent_oracle_evaluations[1];
    const FieldT &C_X = constituent_oracle_evaluations[2];
    const FieldT &Z_X_inv = this->Z_.evaluation_at_point(evaluation_point).inverse();

    return (Z_X_inv*(A_X*B_X-C_X));
}

} // libiop
