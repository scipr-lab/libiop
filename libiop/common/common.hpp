/**@file
 *****************************************************************************
 Common routines not already present in libff.
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_COMMON_COMMON_HPP_
#define LIBIOP_COMMON_COMMON_HPP_

#include <cstddef>
#include <cstdio>
#include <cstdint>
#include <iostream>
#include <vector>
#include <math.h>

namespace libiop {

long double add_soundness_error_bits(const size_t bits1, const size_t bits2);
long double add_soundness_error_bits(const long double bits1, const long double bits2);

} // namespace libiop

#endif // LIBIOP_COMMON_COMMON_HPP_
