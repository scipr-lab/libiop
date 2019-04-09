#include <stdexcept>

namespace libiop {

template<typename FieldT>
r1cs_sparse_matrix<FieldT>::r1cs_sparse_matrix(
    std::shared_ptr<r1cs_constraint_system<FieldT> > constraint_system,
    const r1cs_sparse_matrix_type matrix_type) :
    constraint_system_(constraint_system),
    matrix_type_(matrix_type)
{
}

template<typename FieldT>
linear_combination<FieldT> r1cs_sparse_matrix<FieldT>::get_row(const std::size_t row_index) const
{
    if (row_index >= this->num_rows())
    {
        throw std::invalid_argument("Requested row out of bounds.");
    }

    switch (this->matrix_type_)
    {
    case r1cs_sparse_matrix_A:
        return this->constraint_system_->constraints_[row_index].a_;
    case r1cs_sparse_matrix_B:
        return this->constraint_system_->constraints_[row_index].b_;
    case r1cs_sparse_matrix_C:
        return this->constraint_system_->constraints_[row_index].c_;
    default:
        throw std::logic_error("Invalid matrix type.");
    }
}

template<typename FieldT>
std::size_t r1cs_sparse_matrix<FieldT>::num_rows() const
{
    return this->constraint_system_->num_constraints();
}

template<typename FieldT>
std::size_t r1cs_sparse_matrix<FieldT>::num_columns() const
{
    return this->constraint_system_->num_variables() + 1;
}

} // libiop
