#include <vector>
#include <benchmark/benchmark.h>

#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>

#include "libiop/algebra/utils.hpp"
#include "libiop/bcs/hashing/blake2b.hpp"
#include "libiop/bcs/hashing/algebraic_sponge.hpp"
#include "libiop/bcs/hashing/poseidon.hpp"

namespace libiop {

static void BM_blake2b(benchmark::State &state)
{
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;
    
    const size_t sz = state.range(0);

    blake2b_leafhash<FieldT> leafhasher(128);
    const std::vector<FieldT> avec = random_vector<FieldT>(sz);

    std::vector<FieldT> cvec(sz);

    for (auto _ : state)
    {
        leafhasher.hash(avec);
    }

    state.SetItemsProcessed(state.iterations() * sz);
}

BENCHMARK(BM_blake2b)->RangeMultiplier(2)->Range(1, 16)->Unit(benchmark::kNanosecond);

static void BM_Starkware_poseidon(benchmark::State &state)
{
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;
    
    const size_t sz = state.range(0);

    poseidon_params<FieldT> params = default_128_bit_altbn_poseidon_params<FieldT>();
    poseidon<FieldT> hash(params);
    algebraic_leafhash<FieldT> leafhasher(std::make_shared<poseidon<FieldT>>(hash), 127);
    const std::vector<FieldT> avec = random_vector<FieldT>(sz);

    std::vector<FieldT> cvec(sz);

    for (auto _ : state)
    {
        leafhasher.hash(avec);
    }

    state.SetItemsProcessed(state.iterations() * sz);
}

BENCHMARK(BM_Starkware_poseidon)->RangeMultiplier(2)->Range(1, 16)->Unit(benchmark::kNanosecond);

static void BM_high_alpha_poseidon(benchmark::State &state)
{
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;

    const size_t sz = state.range(0);

    poseidon_params<FieldT> params = high_alpha_128_bit_altbn_poseidon_params<FieldT>();
    poseidon<FieldT> hash(params);
    algebraic_leafhash<FieldT> leafhasher(std::make_shared<poseidon<FieldT>>(hash), 127);
    const std::vector<FieldT> avec = random_vector<FieldT>(sz);

    std::vector<FieldT> cvec(sz);

    for (auto _ : state)
    {
        leafhasher.hash(avec);
    }

    state.SetItemsProcessed(state.iterations() * sz);
}

BENCHMARK(BM_high_alpha_poseidon)->RangeMultiplier(2)->Range(1, 32)->Unit(benchmark::kNanosecond);

static void BM_high_alpha_poseidon_state_size_4(benchmark::State &state)
{
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;

    const size_t sz = state.range(0);

    poseidon_params<FieldT> params = high_alpha_128_bit_altbn_poseidon_params<FieldT>(4);
    poseidon<FieldT> hash(params);
    algebraic_leafhash<FieldT> leafhasher(std::make_shared<poseidon<FieldT>>(hash), 127);
    const std::vector<FieldT> avec = random_vector<FieldT>(sz);

    std::vector<FieldT> cvec(sz);

    for (auto _ : state)
    {
        leafhasher.hash(avec);
    }

    state.SetItemsProcessed(state.iterations() * sz);
}

BENCHMARK(BM_high_alpha_poseidon_state_size_4)->RangeMultiplier(2)->Range(1, 32)->Unit(benchmark::kNanosecond);

}

BENCHMARK_MAIN();
