#include <vector>
#include <benchmark/benchmark.h>

#include <libff/algebra/field_utils/field_utils.hpp>
#include <libff/common/utils.hpp>
#include "libiop/algebra/utils.hpp"
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>

namespace libiop {

static void BM_alt_bn128_mul_vec(benchmark::State &state)
{
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;
    const size_t sz = state.range(0);
    const std::vector<FieldT> avec = random_vector<FieldT>(sz);
    const std::vector<FieldT> bvec = random_vector<FieldT>(sz);

    std::vector<FieldT> cvec(sz);

    for (auto _ : state)
    {
        for (size_t i = 0; i < sz; ++i)
        {
            cvec[i] = avec[i] * bvec[i];
        }
    }

    state.SetItemsProcessed(state.iterations() * sz);
}

BENCHMARK(BM_alt_bn128_mul_vec)->Range(1<<10, 1<<20)->Unit(benchmark::kMicrosecond);

static void BM_alt_bn128_mul_vec_data_dependency(benchmark::State &state)
{
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;
    const size_t sz = state.range(0);
    const std::vector<FieldT> avec = random_vector<FieldT>(sz);
    const std::vector<FieldT> bvec = random_vector<FieldT>(sz);

    FieldT sum = FieldT::zero();

    for (auto _ : state)
    {
        for (size_t i = 0; i < sz; ++i)
        {
            sum += avec[i] * bvec[i];
        }
    }

    state.SetItemsProcessed(state.iterations() * sz);
}

BENCHMARK(BM_alt_bn128_mul_vec_data_dependency)->Range(1<<10, 1<<20)->Unit(benchmark::kMicrosecond);

static void BM_alt_bn128_inverse_vec(benchmark::State& state)
{
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;
    const size_t sz = state.range(0);
    const std::vector<FieldT> vec = random_vector<FieldT>(sz);

    std::vector<FieldT> result(sz);

    for (auto _ : state)
    {
        for (size_t i = 0; i < sz; ++i)
        {
            result[i] = vec[i].inverse();
        }
    }

    state.SetItemsProcessed(state.iterations() * sz);
}

BENCHMARK(BM_alt_bn128_inverse_vec)->Range(1<<10, 1<<16)->Unit(benchmark::kMicrosecond);

}

BENCHMARK_MAIN();
