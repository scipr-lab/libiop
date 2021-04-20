add_library(
  iop

  common/common.cpp
  
  bcs/hashing/blake2b.cpp
  protocols/ldt/ldt_reducer.cpp
  protocols/ldt/fri/fri_ldt.cpp
  protocols/ldt/fri/fri_aux.cpp
  relations/sparse_matrix.cpp
  iop/utilities/batching.cpp
  algebra/utils.cpp
)

# Cmake find modules
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(sodium)

target_link_libraries(
  iop

  sodium

  ff
)

target_include_directories(
  iop

  PUBLIC
  ../depends/libfqfft
  ../depends/libff
)

# BENCHMARKING

add_executable(benchmark_gf64 benchmarks/benchmark_gf64.cpp)
target_link_libraries(benchmark_gf64 iop benchmark)

add_executable(benchmark_gf128 benchmarks/benchmark_gf128.cpp)
target_link_libraries(benchmark_gf128 iop benchmark)

add_executable(benchmark_gf192 benchmarks/benchmark_gf192.cpp)
target_link_libraries(benchmark_gf192 iop benchmark)

add_executable(benchmark_edwards benchmarks/benchmark_edwards.cpp)
target_link_libraries(benchmark_edwards iop benchmark)

add_executable(benchmark_alt_bn128 benchmarks/benchmark_alt_bn128.cpp)
target_link_libraries(benchmark_alt_bn128 iop benchmark)

add_executable(benchmark_gf256 benchmarks/benchmark_gf256.cpp)
target_link_libraries(benchmark_gf256 iop benchmark)

add_executable(benchmark_hashes benchmarks/benchmark_hashes.cpp)
target_link_libraries(benchmark_hashes iop benchmark)

add_executable(benchmark_fft benchmarks/benchmark_fft.cpp)
target_link_libraries(benchmark_fft iop benchmark)

add_executable(benchmark_lagrange benchmarks/benchmark_lagrange.cpp)
target_link_libraries(benchmark_lagrange iop benchmark)

add_executable(benchmark_polynomials benchmarks/benchmark_polynomials.cpp)
target_link_libraries(benchmark_polynomials iop benchmark)

add_executable(benchmark_sumcheck benchmarks/benchmark_sumcheck.cpp)
target_link_libraries(benchmark_sumcheck iop benchmark gtest_main)

add_executable(benchmark_vector_op benchmarks/benchmark_vector_op.cpp)
target_link_libraries(benchmark_vector_op iop benchmark)

# INSTRUMENTATION

add_executable(instrument_algebra profiling/instrument_algebra.cpp)
add_executable(instrument_aurora_snark profiling/instrument_aurora_snark.cpp)
add_executable(instrument_fractal_snark profiling/instrument_fractal_snark.cpp)
add_executable(instrument_fri_snark profiling/instrument_fri_snark.cpp)
add_executable(instrument_ligero_snark profiling/instrument_ligero_snark.cpp)
if("${CPPDEBUG}")
  target_link_libraries(instrument_aurora_snark iop gtest_main)
  target_link_libraries(instrument_fractal_snark iop gtest_main)
  target_link_libraries(instrument_fri_snark iop)
  target_link_libraries(instrument_ligero_snark iop)
else()
  target_link_libraries(instrument_aurora_snark iop ${Boost_LIBRARIES} gtest_main)
  target_link_libraries(instrument_fractal_snark iop ${Boost_LIBRARIES} gtest_main)
  target_link_libraries(instrument_fri_snark iop ${Boost_LIBRARIES})
  target_link_libraries(instrument_ligero_snark iop ${Boost_LIBRARIES})
endif()
target_link_libraries(instrument_algebra iop ${Boost_LIBRARIES})

# TESTS

include(CTest)

# algebra
add_executable(test_exponentiation tests/algebra/test_exponentiation.cpp)
target_link_libraries(test_exponentiation iop gtest_main)

add_executable(test_fft tests/algebra/test_fft.cpp)
target_link_libraries(test_fft iop gtest_main)

add_executable(test_lagrange tests/algebra/test_lagrange.cpp)
target_link_libraries(test_lagrange iop gtest_main)

add_executable(test_merkle_tree tests/bcs/test_merkle_tree.cpp)
target_link_libraries(test_merkle_tree iop gtest_main)

add_executable(test_linearized_polynomial tests/algebra/test_linearized_polynomial.cpp)
target_link_libraries(test_linearized_polynomial iop gtest_main)

add_executable(test_vanishing_polynomial tests/algebra/test_vanishing_polynomial.cpp)
target_link_libraries(test_vanishing_polynomial iop gtest_main)

add_executable(test_lagrange_polynomial tests/algebra/test_lagrange_polynomial.cpp)
target_link_libraries(test_lagrange_polynomial iop gtest_main)

add_executable(test_bivariate_embedding tests/algebra/test_bivariate_embedding.cpp)
target_link_libraries(test_bivariate_embedding iop gtest_main)

add_executable(test_successor_ordering tests/algebra/test_successor_ordering.cpp)
target_link_libraries(test_successor_ordering iop gtest_main)

add_executable(test_algebra_utils tests/algebra/test_algebra_utils.cpp)
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
add_executable(test_iop tests/iop/test_iop.cpp)
target_link_libraries(test_iop iop gtest_main)

add_executable(test_iop_query_position tests/iop/test_iop_query_position.cpp)
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
add_executable(test_aurora_protocol tests/protocols/test_aurora_protocol.cpp)
target_link_libraries(test_aurora_protocol iop gtest_main)

add_executable(test_r1cs_to_lincheck_reduction tests/protocols/test_r1cs_to_lincheck_reduction.cpp)
target_link_libraries(test_r1cs_to_lincheck_reduction iop gtest_main)

add_executable(test_basic_lincheck tests/protocols/test_basic_lincheck.cpp)
target_link_libraries(test_basic_lincheck iop gtest_main)

add_executable(test_boundary_constraint tests/protocols/test_boundary_constraint.cpp)
target_link_libraries(test_boundary_constraint iop gtest_main)

add_executable(test_holographic_lincheck tests/protocols/test_holographic_lincheck.cpp)
target_link_libraries(test_holographic_lincheck iop gtest_main)

add_executable(test_sumcheck tests/protocols/test_sumcheck.cpp)
target_link_libraries(test_sumcheck iop gtest_main)

add_executable(test_direct_ldt tests/protocols/test_direct_ldt.cpp)
target_link_libraries(test_direct_ldt iop gtest_main)

add_executable(test_fri tests/protocols/test_fri.cpp)
target_link_libraries(test_fri iop gtest_main)

add_executable(test_fri_aux tests/protocols/test_fri_aux.cpp)
target_link_libraries(test_fri_aux iop gtest_main)

add_executable(test_fri_optimizer tests/protocols/test_fri_optimizer.cpp)
target_link_libraries(test_fri_optimizer iop gtest_main)

add_executable(test_fractal_protocol tests/protocols/test_fractal_protocol.cpp)
target_link_libraries(test_fractal_protocol iop gtest_main)

add_executable(test_ligero_interleaved_lincheck_et tests/protocols/test_ligero_interleaved_lincheck_et.cpp)
target_link_libraries(test_ligero_interleaved_lincheck_et iop gtest_main)

add_executable(test_ligero_interleaved_lincheck_ot tests/protocols/test_ligero_interleaved_lincheck_ot.cpp)
target_link_libraries(test_ligero_interleaved_lincheck_ot iop gtest_main)

add_executable(test_ligero_protocol tests/protocols/test_ligero_protocol.cpp)
target_link_libraries(test_ligero_protocol iop gtest_main)

add_executable(test_ligero_interleaved_rowcheck tests/protocols/test_ligero_interleaved_rowcheck.cpp)
target_link_libraries(test_ligero_interleaved_rowcheck iop gtest_main)

add_executable(test_ldt_reducer tests/protocols/test_ldt_reducer.cpp)
target_link_libraries(test_ldt_reducer iop gtest_main)

add_executable(test_fractal_indexer tests/protocols/test_fractal_indexer.cpp)
target_link_libraries(test_fractal_indexer iop gtest_main)

add_executable(test_rational_linear_combination tests/protocols/test_rational_linear_combination.cpp)
target_link_libraries(test_rational_linear_combination iop gtest_main)

add_executable(test_rational_sumcheck tests/protocols/test_rational_sumcheck.cpp)
target_link_libraries(test_rational_sumcheck iop gtest_main)

add_executable(test_rowcheck tests/protocols/test_rowcheck.cpp)
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

add_executable(test_r1cs tests/relations/test_r1cs.cpp)
target_link_libraries(test_r1cs iop gtest_main)

add_executable(test_identity_matrices tests/relations/test_identity_matrices.cpp)
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

add_executable(test_poseidon tests/snark/test_poseidon.cpp)
target_link_libraries(test_poseidon iop gtest_main)

add_executable(test_pow tests/snark/test_pow.cpp)
target_link_libraries(test_pow iop gtest_main)

add_executable(test_bcs_transformation tests/bcs/test_bcs_transformation.cpp)
target_link_libraries(test_bcs_transformation iop gtest_main)

add_executable(test_serialization tests/snark/test_serialization.cpp)
target_link_libraries(test_serialization iop gtest_main)

add_executable(test_aurora_snark tests/snark/test_aurora_snark.cpp)
target_link_libraries(test_aurora_snark iop gtest_main)

add_executable(test_fractal_snark tests/snark/test_fractal_snark.cpp )
target_link_libraries(test_fractal_snark iop gtest_main)

add_executable(test_ligero_snark tests/snark/test_ligero_snark.cpp)
target_link_libraries(test_ligero_snark iop gtest_main)

add_executable(test_linking tests/snark/test_linking.cpp)
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
