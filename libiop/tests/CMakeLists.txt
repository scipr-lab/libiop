#Tests 
include(CTest)

# algebra
add_executable(test_exponentiation algebra/test_exponentiation.cpp)
target_link_libraries(test_exponentiation iop gtest_main)

add_executable(test_fft algebra/test_fft.cpp)
target_link_libraries(test_fft iop gtest_main)

add_executable(test_lagrange algebra/test_lagrange.cpp)
target_link_libraries(test_lagrange iop gtest_main)

add_executable(test_merkle_tree bcs/test_merkle_tree.cpp)
target_link_libraries(test_merkle_tree iop gtest_main)

add_executable(test_linearized_polynomial algebra/test_linearized_polynomial.cpp)
target_link_libraries(test_linearized_polynomial iop gtest_main)

add_executable(test_vanishing_polynomial algebra/test_vanishing_polynomial.cpp)
target_link_libraries(test_vanishing_polynomial iop gtest_main)

add_executable(test_lagrange_polynomial algebra/test_lagrange_polynomial.cpp)
target_link_libraries(test_lagrange_polynomial iop gtest_main)

add_executable(test_bivariate_embedding algebra/test_bivariate_embedding.cpp)
target_link_libraries(test_bivariate_embedding iop gtest_main)

add_executable(test_successor_ordering algebra/test_successor_ordering.cpp)
target_link_libraries(test_successor_ordering iop gtest_main)

add_executable(test_algebra_utils algebra/test_algebra_utils.cpp)
target_link_libraries(test_algebra_utils iop gtest_main)

add_test(
  NAME test_exponentiation
  COMMAND test_exponentiation
)
add_test(
  NAME test_fft
  COMMAND test_fft
)
add_test(
  NAME test_lagrange
  COMMAND test_lagrange
)
add_test(
  NAME test_merkle_tree
  COMMAND test_merkle_tree
)
add_test(
  NAME test_linearized_polynomial
  COMMAND test_linearized_polynomial
)
add_test(
  NAME test_vanishing_polynomial
  COMMAND test_vanishing_polynomial
)
add_test(
  NAME test_lagrange_polynomial
  COMMAND test_lagrange_polynomial
)
add_test(
  NAME test_bivariate_embedding
  COMMAND test_bivariate_embedding
)
add_test(
  NAME test_successor_ordering
  COMMAND test_successor_ordering
)
add_test(
  NAME test_algebra_utils
  COMMAND test_algebra_utils
)

# iop
add_executable(test_iop iop/test_iop.cpp)
target_link_libraries(test_iop iop gtest_main)

add_executable(test_iop_query_position iop/test_iop_query_position.cpp)
target_link_libraries(test_iop_query_position iop gtest_main)

add_test(
  NAME test_iop
  COMMAND test_iop
)
add_test(
  NAME test_iop_query_position
  COMMAND test_iop_query_position
)

# protocols
add_executable(test_aurora_protocol protocols/test_aurora_protocol.cpp)
target_link_libraries(test_aurora_protocol iop gtest_main)

add_executable(test_r1cs_to_lincheck_reduction protocols/test_r1cs_to_lincheck_reduction.cpp)
target_link_libraries(test_r1cs_to_lincheck_reduction iop gtest_main)

add_executable(test_basic_lincheck protocols/test_basic_lincheck.cpp)
target_link_libraries(test_basic_lincheck iop gtest_main)

add_executable(test_boundary_constraint protocols/test_boundary_constraint.cpp)
target_link_libraries(test_boundary_constraint iop gtest_main)

add_executable(test_holographic_lincheck protocols/test_holographic_lincheck.cpp)
target_link_libraries(test_holographic_lincheck iop gtest_main)

add_executable(test_sumcheck protocols/test_sumcheck.cpp)
target_link_libraries(test_sumcheck iop gtest_main)

add_executable(test_direct_ldt protocols/test_direct_ldt.cpp)
target_link_libraries(test_direct_ldt iop gtest_main)

add_executable(test_fri protocols/test_fri.cpp)
target_link_libraries(test_fri iop gtest_main)

add_executable(test_fri_aux protocols/test_fri_aux.cpp)
target_link_libraries(test_fri_aux iop gtest_main)

add_executable(test_fri_optimizer protocols/test_fri_optimizer.cpp)
target_link_libraries(test_fri_optimizer iop gtest_main)

add_executable(test_fractal_protocol protocols/test_fractal_protocol.cpp)
target_link_libraries(test_fractal_protocol iop gtest_main)

add_executable(test_ligero_interleaved_lincheck_et protocols/test_ligero_interleaved_lincheck_et.cpp)
target_link_libraries(test_ligero_interleaved_lincheck_et iop gtest_main)

add_executable(test_ligero_interleaved_lincheck_ot protocols/test_ligero_interleaved_lincheck_ot.cpp)
target_link_libraries(test_ligero_interleaved_lincheck_ot iop gtest_main)

add_executable(test_ligero_protocol protocols/test_ligero_protocol.cpp)
target_link_libraries(test_ligero_protocol iop gtest_main)

add_executable(test_ligero_interleaved_rowcheck protocols/test_ligero_interleaved_rowcheck.cpp)
target_link_libraries(test_ligero_interleaved_rowcheck iop gtest_main)

add_executable(test_ldt_reducer protocols/test_ldt_reducer.cpp)
target_link_libraries(test_ldt_reducer iop gtest_main)

add_executable(test_fractal_indexer protocols/test_fractal_indexer.cpp)
target_link_libraries(test_fractal_indexer iop gtest_main)

add_executable(test_rational_linear_combination protocols/test_rational_linear_combination.cpp)
target_link_libraries(test_rational_linear_combination iop gtest_main)

add_executable(test_rational_sumcheck protocols/test_rational_sumcheck.cpp)
target_link_libraries(test_rational_sumcheck iop gtest_main)

add_executable(test_rowcheck protocols/test_rowcheck.cpp)
target_link_libraries(test_rowcheck iop gtest_main)

add_test(
  NAME test_aurora_protocol
  COMMAND test_aurora_protocol
)
add_test(
  NAME test_r1cs_to_lincheck_reduction
  COMMAND test_r1cs_to_lincheck_reduction
)
add_test(
  NAME test_basic_lincheck
  COMMAND test_basic_lincheck
)
add_test(
  NAME test_boundary_constraint
  COMMAND test_boundary_constraint
)
add_test(
  NAME test_holographic_lincheck
  COMMAND test_holographic_lincheck
)
add_test(
  NAME test_sumcheck
  COMMAND test_sumcheck
)
add_test(
  NAME test_direct_ldt
  COMMAND test_direct_ldt
)
add_test(
  NAME test_fri
  COMMAND test_fri
)
add_test(
  NAME test_fri_aux
  COMMAND test_fri_aux
)
add_test(
  NAME test_fri_optimizer
  COMMAND test_fri_optimizer
)
add_test(
  NAME test_fractal_protocol
  COMMAND test_fractal_protocol
)
add_test(
  NAME test_ligero_interleaved_lincheck_et
  COMMAND test_ligero_interleaved_lincheck_et
)
add_test(
  NAME test_ligero_interleaved_lincheck_ot
  COMMAND test_ligero_interleaved_lincheck_ot
)
add_test(
  NAME test_ligero_protocol
  COMMAND test_ligero_protocol
)
add_test(
  NAME test_ligero_interleaved_rowcheck
  COMMAND test_ligero_interleaved_rowcheck
)
add_test(
  NAME test_ldt_reducer
  COMMAND test_ldt_reducer
)
add_test(
  NAME test_fractal_indexer
  COMMAND test_fractal_indexer
)
add_test(
  NAME test_rational_linear_combination
  COMMAND test_rational_linear_combination
)
add_test(
  NAME test_rational_sumcheck
  COMMAND test_rational_sumcheck
)
add_test(
  NAME test_rowcheck
  COMMAND test_rowcheck
)

# relations

add_executable(test_r1cs relations/test_r1cs.cpp)
target_link_libraries(test_r1cs iop gtest_main)

add_executable(test_identity_matrices relations/test_identity_matrices.cpp)
target_link_libraries(test_identity_matrices iop gtest_main)

add_test(
  NAME test_r1cs
  COMMAND test_r1cs
)
add_test(
  NAME test_identity_matrices
  COMMAND test_identity_matrices
)

# snark

add_executable(test_poseidon snark/test_poseidon.cpp)
target_link_libraries(test_poseidon iop gtest_main)

add_executable(test_pow snark/test_pow.cpp)
target_link_libraries(test_pow iop gtest_main)

add_executable(test_bcs_transformation bcs/test_bcs_transformation.cpp)
target_link_libraries(test_bcs_transformation iop gtest_main)

add_executable(test_serialization snark/test_serialization.cpp)
target_link_libraries(test_serialization iop gtest_main)

add_executable(test_aurora_snark snark/test_aurora_snark.cpp)
target_link_libraries(test_aurora_snark iop gtest_main)

add_executable(test_fractal_snark snark/test_fractal_snark.cpp )
target_link_libraries(test_fractal_snark iop gtest_main)

add_executable(test_ligero_snark snark/test_ligero_snark.cpp)
target_link_libraries(test_ligero_snark iop gtest_main)

add_executable(test_linking snark/test_linking.cpp)
target_link_libraries(test_linking iop gtest_main)

add_test(
  NAME test_poseidon
  COMMAND test_poseidon
)
add_test(
  NAME test_pow
  COMMAND test_pow
)
add_test(
  NAME test_bcs_transformation
  COMMAND test_bcs_transformation
)
add_test(
  NAME test_serialization
  COMMAND test_serialization
)
add_test(
  NAME test_aurora_snark
  COMMAND test_aurora_snark
)
add_test(
  NAME test_fractal_snark
  COMMAND test_fractal_snark
)
add_test(
  NAME test_ligero_snark
  COMMAND test_ligero_snark
)
add_test(
  NAME test_linking
  COMMAND test_linking
)
add_custom_target(tests DEPENDS test_exponentiation test_fft test_lagrange test_merkle_tree test_linearized_polynomial test_vanishing_polynomial test_lagrange_polynomial test_bivariate_embedding test_successor_ordering test_algebra_utils test_iop test_iop_query_position test_aurora_protocol test_r1cs_to_lincheck_reduction test_basic_lincheck test_boundary_constraint test_holographic_lincheck test_sumcheck test_direct_ldt test_fri test_fri_aux test_fri_optimizer test_fractal_protocol test_ligero_interleaved_lincheck_et test_ligero_interleaved_lincheck_ot test_ligero_protocol test_ligero_interleaved_rowcheck test_ldt_reducer test_fractal_indexer test_rational_linear_combination test_rational_sumcheck test_rowcheck test_r1cs test_identity_matrices test_poseidon test_pow test_bcs_transformation test_serialization test_aurora_snark test_fractal_snark test_ligero_snark test_linking)
