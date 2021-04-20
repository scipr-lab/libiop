#include <type_traits>
#include "libiop/bcs/bcs_common.hpp"
#include "libiop/bcs/hashing/blake2b.hpp"
#include "libiop/bcs/hashing/hash_enum.hpp"
#include "libiop/bcs/pow.hpp"

namespace libiop {

template<typename FieldT, typename MT_root_hash>
bcs_transformation_parameters<FieldT, MT_root_hash> default_bcs_params(
    const bcs_hash_type hash_type, const std::size_t security_parameter, const size_t dim_h)
{
    bcs_transformation_parameters<FieldT, MT_root_hash> params;
    params.security_parameter = security_parameter;
    params.hash_enum = hash_type;
    /* TODO: Push setting leaf hash into internal BCS code. Currently 2 is fine, as leaf size is internally unused. */
    const size_t leaf_size = 2;
    params.leafhasher_ = get_leafhash<FieldT, MT_root_hash>(hash_type, security_parameter, leaf_size);
    params.compression_hasher = get_two_to_one_hash<MT_root_hash, FieldT>(hash_type, security_parameter);
    params.hashchain_ =
        get_hashchain<FieldT, MT_root_hash>(hash_type, security_parameter);

    // Work per hash. Todo generalize this w/ proper explanations of work amounts
    const size_t work_per_hash = (hash_type == 1) ? 1 : 128;
    params.pow_params_ = pow_parameters(dim_h + 3 + libff::log2(work_per_hash), work_per_hash);
    return params;
}

} // namespace libiop
