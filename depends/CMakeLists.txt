add_subdirectory(gtest EXCLUDE_FROM_ALL)

option(BENCHMARK_ENABLE_GTEST_TESTS "" OFF)
set(BENCHMARK_ENABLE_TESTING OFF)
add_subdirectory(benchmark)

include_directories(ate-pairing/include)
include_directories(xbyak)
add_library(
        zm
        STATIC

        ate-pairing/src/zm.cpp
        ate-pairing/src/zm2.cpp
)

OPTION(IS_LIBFF_PARENT OFF)
add_subdirectory(libff)
if(NOT TARGET ff)
    add_subdirectory(libff)
endif()

OPTION(IS_LIBFQFFT_PARENT OFF)
if(NOT TARGET fqfft)
    add_subdirectory(libfqfft)
endif()
