/**@file
*****************************************************************************
Aurora protocol for R1CS
*****************************************************************************
* @author     This file is part of libiop (see AUTHORS)
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/
#ifndef LIBIOP_PROTOCOLS_ENCODED_AURORA_AURORA_HPP_
#define LIBIOP_PROTOCOLS_ENCODED_AURORA_AURORA_HPP_

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>
#include "gtest/gtest_prod.h"

#include "libiop/iop/iop.hpp"
#include "libiop/protocols/encoded/aurora/rowcheck.hpp"
#include "libiop/protocols/encoded/aurora/sumcheck.hpp"
#include "libiop/protocols/encoded/aurora/multi_lincheck.hpp"
#include "libiop/relations/r1cs.hpp"

namespace libiop {

template<typename FieldT>
class fz_virtual_oracle;

template<typename FieldT>
class encoded_aurora_parameters {
protected:
    size_t interactive_security_parameter_;

    size_t codeword_domain_dim_;
    size_t constraint_domain_dim_;
    size_t summation_domain_dim_;
    size_t query_bound_;
    bool make_zk_;
    field_subset_type domain_type_;

public:
    encoded_aurora_parameters() {};
    encoded_aurora_parameters(const size_t interactive_security_parameter,
                              const size_t codeword_domain_dim,
                              const size_t constraint_domain_dim,
                              const size_t summation_domain_dim,
                              const size_t query_bound,
                              const bool make_zk,
                              const field_subset_type type);

    size_t max_tested_degree_bound() const;
    size_t max_constraint_degree_bound() const;
    size_t locality() const;

    bool make_zk() const;
    size_t query_bound() const;

    void print();

    multi_lincheck_parameters<FieldT> multi_lincheck_params_;
};

template<typename FieldT>
class encoded_aurora_protocol {
protected:
    iop_protocol<FieldT> &IOP_;
    domain_handle constraint_domain_handle_;
    domain_handle variable_domain_handle_;
    domain_handle codeword_domain_handle_;
    r1cs_constraint_system<FieldT> constraint_system_;
    encoded_aurora_parameters<FieldT> params_;

    domain_handle summation_domain_handle_; /* set either to constraint or variable domain depending which is larger */

    field_subset<FieldT> constraint_domain_,
        variable_domain_,
        codeword_domain_,
        summation_domain_,
        input_variable_domain_; /* mostly for convenience to avoid calling this->IOP_.get_domain(..._handle) every time */

    oracle_handle fw_handle_;
    oracle_handle fAz_handle_, fBz_handle_, fCz_handle_;

    std::shared_ptr<multi_lincheck<FieldT> > multi_lincheck_;

    std::size_t fw_mask_degree_;
    polynomial<FieldT> fw_mask_;
    std::vector<FieldT> f_1v_coefficients_;
    std::vector<FieldT> fw_over_codeword_domain_;

    std::vector<FieldT> Az_, Bz_, Cz_;
    polynomial<FieldT> R_Az_, R_Bz_, R_Cz_;
    std::vector<FieldT> fprime_Az_over_codeword_domain_,
        fprime_Bz_over_codeword_domain_,
        fprime_Cz_over_codeword_domain_;

    std::shared_ptr<r1cs_sparse_matrix<FieldT> > r1cs_A_, r1cs_B_, r1cs_C_;

    std::shared_ptr<fz_virtual_oracle<FieldT> > fz_oracle_;
    virtual_oracle_handle fz_oracle_handle_;

    std::shared_ptr<rowcheck_ABC_virtual_oracle<FieldT> > rowcheck_oracle_;
    virtual_oracle_handle rowcheck_oracle_handle_;

    void calculate_ABCz(const std::vector<FieldT> &fz_over_variable_domain);
    void compute_fprime_ABCz_over_codeword_domain();
    void register_witness_oracles();
    FRIEND_TEST(TestAdditiveR1CSComponents, R1CSTest);
    FRIEND_TEST(TestMultiplicativeR1CSComponents, R1CSTest);

public:
    /* Initialization and registration */
    encoded_aurora_protocol(iop_protocol<FieldT> &IOP,
                            const domain_handle &constraint_domain_handle,
                            const domain_handle &variable_domain_handle,
                            const domain_handle &codeword_domain_handle,
                            const r1cs_constraint_system<FieldT> &constraint_system,
                            const encoded_aurora_parameters<FieldT> &params);

    void register_challenge();
    void register_proof();

    /* Proving */
    void submit_witness_oracles(const r1cs_primary_input<FieldT> &primary_input,
                                const r1cs_auxiliary_input<FieldT> &auxiliary_input);
    void calculate_and_submit_proof();

    /* Verification */
    void construct_verifier_state(const r1cs_primary_input<FieldT> &primary_input);

    void print_soundness_error(const std::size_t LDT_query_repetitions,
                               const std::size_t RS_extra_dimensions);

    std::vector<oracle_handle_ptr> get_all_oracle_handles();
};

} // namespace libiop

#include "libiop/protocols/encoded/aurora/aurora.tcc"

#endif // LIBIOP_PROTOCOLS_ENCODED_AURORA_AURORA_HPP_
