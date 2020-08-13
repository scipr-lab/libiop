/**@file
 *****************************************************************************
 Specialized IOP that implements the BCS16 transformation.
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_SNARK_COMMON_COMMON_BCS_PARARMETERS_HPP_
#define LIBIOP_SNARK_COMMON_COMMON_BCS_PARARMETERS_HPP_

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <map>
#include <vector>

#include "libiop/bcs/bcs_common.hpp"
#include "libiop/bcs/hashing/hash_enum.hpp"
#include "libiop/bcs/hashing/hashing.hpp"
#include "libiop/bcs/hashing/blake2b.hpp"
#include "libiop/bcs/pow.hpp"

namespace libiop {

template<typename FieldT, typename MT_root_hash>
bcs_transformation_parameters<FieldT, MT_root_hash> default_bcs_params(
    const bcs_hash_type hash_type, const size_t security_parameter);

} // namespace libiop

#include "libiop/bcs/common_bcs_parameters.tcc"

#endif // LIBIOP_SNARK_COMMON_COMMON_BCS_PARARMETERS_HPP_
