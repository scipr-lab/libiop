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
set(IOP_HEADERS
  common/common.hpp
  
  bcs/hashing/blake2b.hpp
  protocols/ldt/ldt_reducer.hpp
  protocols/ldt/fri/fri_ldt.hpp
  protocols/ldt/fri/fri_aux.hpp
  relations/sparse_matrix.hpp
  iop/utilities/batching.hpp
  algebra/utils.hpp
)
message(STATUS "IOP_HEADERS: ${IOP_HEADERS}")
install(
    DIRECTORY ""
    DESTINATION "include/libiop"
    FILES "${IOP_HEADERS}"
)
install(TARGETS iop DESTINATION lib)
# add_subdirectory(benchmarks EXCLUDE_FROM_ALL)
# add_subdirectory(profiling EXCLUDE_FROM_ALL)
# add_subdirectory(tests EXCLUDE_FROM_ALL)
