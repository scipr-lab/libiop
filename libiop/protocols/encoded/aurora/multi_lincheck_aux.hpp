/**@file
*****************************************************************************
Aurora R1CS multi lincheck virtual oracle auxiliary classes.
*****************************************************************************
* @author     This file is part of libiop (see AUTHORS)
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/
#ifndef LIBIOP_PROTOCOLS_ENCODED_AURORA_MULTI_LINCHECK_AUX_HPP_
#define LIBIOP_PROTOCOLS_ENCODED_AURORA_MULTI_LINCHECK_AUX_HPP_

#include <cstring>
#include <cstddef>
#include <map>
#include <memory>
#include <vector>

#include "libiop/relations/sparse_matrix.hpp"
#include "libiop/algebra/lagrange.hpp"
#include "libiop/iop/iop.hpp"


namespace libiop {

template<typename FieldT>
class multi_lincheck_virtual_oracle : public virtual_oracle<FieldT> {
protected:
    const field_subset<FieldT> codeword_domain_;
    const field_subset<FieldT> constraint_domain_;
    const field_subset<FieldT> variable_domain_;
    const field_subset<FieldT> summation_domain_;
    const std::size_t input_variable_dim_;
    const std::vector<std::shared_ptr<sparse_matrix<FieldT> >> matrices_;

    std::vector<FieldT> alpha_powers_;
    std::vector<FieldT> r_Mz_;
    std::vector<FieldT> p_alpha_ABC_evals_;
    std::shared_ptr<lagrange_cache<FieldT> > lagrange_coefficients_cache_;

public:
    multi_lincheck_virtual_oracle(
        const field_subset<FieldT> &codeword_domain,
        const field_subset<FieldT> &constraint_domain,
        const field_subset<FieldT> &variable_domain,
        const field_subset<FieldT> &summation_domain,
        const std::size_t input_variable_dim,
        const std::vector<std::shared_ptr<sparse_matrix<FieldT> >> &matrices,
        std::shared_ptr<lagrange_cache<FieldT> > lagrange_coefficients_cache);

    void set_challenge(const FieldT &alpha, const std::vector<FieldT> r_Mz);

    virtual std::vector<FieldT> evaluated_contents(
        const std::vector<std::vector<FieldT> > &constituent_oracle_evaluations) const;
    virtual FieldT evaluation_at_point(
        const std::size_t evaluation_position,
        const FieldT evaluation_point,
        const std::vector<FieldT> &constituent_oracle_evaluations) const;
};

} // libiop

#include "libiop/protocols/encoded/aurora/multi_lincheck_aux.tcc"

#endif // LIBIOP_PROTOCOLS_ENCODED_AURORA_MULTI_LINCHECK_AUX_HPP_
