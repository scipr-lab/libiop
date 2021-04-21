/**@file
*****************************************************************************
Full protocol for R1CS (encoded R1CS + FRI LDT)
*****************************************************************************
* @author     This file is part of libiop (see AUTHORS)
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/
#ifndef LIBIOP_PROTOCOLS_FRACTAL_HIOP_HPP_
#define LIBIOP_PROTOCOLS_FRACTAL_HIOP_HPP_

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

#include <libff/common/utils.hpp>
#include "libiop/iop/iop.hpp"
#include "libiop/protocols/encoded/r1cs_rs_iop/r1cs_rs_iop.hpp"
#include "libiop/protocols/encoded/r1cs_rs_iop/fractal_indexer.hpp"
#include "libiop/protocols/ldt/fri/fri_ldt.hpp"
#include "libiop/protocols/ldt/ldt_reducer.hpp"
#include "libiop/relations/r1cs.hpp"

namespace libiop {

/** Currently this file just implements the indexer.*/

/** TODO: Make the parameters actually output the field subsets. */
template<typename FieldT>
class fractal_iop_parameters {
    protected:
    size_t security_parameter_;
    size_t pow_bits_;
    size_t RS_extra_dimensions_;
    bool make_zk_;
    std::shared_ptr<r1cs_constraint_system<FieldT>> constraint_system_;

    field_subset<FieldT> index_domain_;
    field_subset<FieldT> matrix_domain_;
    field_subset<FieldT> codeword_domain_;
    size_t codeword_domain_dim_;

    size_t query_bound_;
    public:
    fractal_iop_parameters() {};
    fractal_iop_parameters(
        const size_t security_parameter,
        const size_t pow_bits,
        const size_t RS_extra_dimensions,
        const bool make_zk,
        const std::shared_ptr<r1cs_constraint_system<FieldT>> &constraint_system);

    void set_ldt_parameters(size_t localization_parameter,
                            FRI_soundness_type fri_soundness_type,
                            LDT_reducer_soundness_type ldt_reducer_soundness_type);
    void set_ldt_parameters(std::vector<size_t> localization_parameters,
                            FRI_soundness_type soundness_type,
                            LDT_reducer_soundness_type ldt_reducer_soundness_type);

    size_t RS_extra_dimensions() const;
    bool make_zk() const;
    std::vector<size_t> locality_vector() const;
    std::shared_ptr<r1cs_constraint_system<FieldT>> constraint_system() const;
    field_subset<FieldT> index_domain() const;
    field_subset<FieldT> matrix_domain() const;
    field_subset<FieldT> codeword_domain() const;

    long double achieved_soundness() const;
    void print() const;

    LDT_instance_reducer_params<FieldT> LDT_reducer_params_;
    FRI_protocol_parameters<FieldT> FRI_params_;
    encoded_aurora_parameters<FieldT> encoded_aurora_params_;
};

template<typename FieldT>
class fractal_iop {
protected:
    iop_protocol<FieldT> &IOP_;

    fractal_iop_parameters<FieldT> parameters_;

    domain_handle index_domain_handle_;
    domain_handle matrix_domain_handle_;
    domain_handle codeword_domain_handle_;

    std::vector<matrix_indexer<FieldT>> matrix_indexers_;
    std::vector<std::vector<oracle_handle_ptr>> indexed_handles_;

    std::shared_ptr<encoded_aurora_protocol<FieldT> > protocol_;
    std::shared_ptr<LDT_instance_reducer<FieldT, FRI_protocol<FieldT>> > LDT_reducer_;
    /* Called within initialization */
    void register_index_oracles();
public:
    /* Initialization and registration */
    fractal_iop(iop_protocol<FieldT> &IOP,
                           const fractal_iop_parameters<FieldT> &parameters);

    void register_interactions();
    void register_queries();

    /* Indexing */
    void produce_index();

    /* Proving */
    void produce_proof(const r1cs_primary_input<FieldT> &primary_input,
                       const r1cs_auxiliary_input<FieldT> &auxiliary_input,
                       iop_prover_index<FieldT> &index);

    /* Verification */
    bool verifier_predicate(const r1cs_primary_input<FieldT> &primary_input);
};


} // namespace libiop

#include "libiop/protocols/fractal_hiop.tcc"

#endif // LIBIOP_PROTOCOLS_FRACTAL_HIOP_HPP_
