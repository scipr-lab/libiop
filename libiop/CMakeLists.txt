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

install(TARGETS iop)
# add_subdirectory(benchmarks EXCLUDE_FROM_ALL)
# add_subdirectory(profiling EXCLUDE_FROM_ALL)
# add_subdirectory(tests EXCLUDE_FROM_ALL)
