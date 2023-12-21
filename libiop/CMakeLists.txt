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
option(BUILD_BENCHMARKS "Option to build benchmarks")
set(benchmarks)
if(BUILD_BENCHMARKS)
   set(benchmarks EXCLUDE_FROM_ALL)
endif()

add_subdirectory(benchmarks ${benchmarks})


option(BUILD_PROFILING "Option to build profiling")
set(profiling)
if(BUILD_PROFILING)
  set(profiling EXCLUDE_FROM_ALL)
endif()

add_subdirectory(profiling ${profiling})

option(BUILD_TESTS "Option to build tests")
set(tests)
if(BUILD_TESTS)
  set(tests EXCLUDE_FROM_ALL)
endif()

add_subdirectory(tests ${tests})