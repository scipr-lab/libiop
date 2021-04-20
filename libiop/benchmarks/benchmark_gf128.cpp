#include <cmath>
#include <cstddef>
#include <vector>
#include <benchmark/benchmark.h>

#include <libff/algebra/fields/binary/gf128.hpp>
#include "libiop/algebra/utils.hpp"

using namespace libff;

namespace libiop {

static void BM_gf128_mul_vec(benchmark::State &state)
{
    const size_t sz = state.range(0);
    const std::vector<libff::gf128> avec = random_vector<libff::gf128>(sz);
    const std::vector<libff::gf128> bvec = random_vector<libff::gf128>(sz);

    std::vector<libff::gf128> cvec(sz);

    for (auto _ : state)
    {
        for (size_t i = 0; i < sz; ++i)
        {
            cvec[i] = avec[i] * bvec[i];
        }
    }

    state.SetItemsProcessed(state.iterations() * sz);
}

BENCHMARK(BM_gf128_mul_vec)->Range(1<<10, 1<<20)->Unit(benchmark::kMicrosecond);

static void BM_gf128_inverse_vec(benchmark::State& state)
{
    const size_t sz = state.range(0);
    const std::vector<libff::gf128> vec = random_vector<libff::gf128>(sz);

    std::vector<libff::gf128> result(sz);

    for (auto _ : state)
    {
        for (size_t i = 0; i < sz; ++i)
        {
            result[i] = vec[i].inverse();
        }
    }

    state.SetItemsProcessed(state.iterations() * sz);
}

BENCHMARK(BM_gf128_inverse_vec)->Range(1<<10, 1<<16)->Unit(benchmark::kMicrosecond);

}

BENCHMARK_MAIN();
