#include "libiop/relations/sparse_matrix.hpp"

namespace libiop {


std::vector<r1cs_sparse_matrix_type> all_r1cs_sparse_matrix_types(
    {r1cs_sparse_matrix_A, r1cs_sparse_matrix_B, r1cs_sparse_matrix_C});

} // libiop
