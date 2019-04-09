/**@file
 *****************************************************************************
  FRI localization parameter array optimizer for minimal proof size
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_PROTOCOLS_LDT_FRI_OPTIMIZER_HPP_
#define LIBIOP_PROTOCOLS_LDT_FRI_OPTIMIZER_HPP_

#include <vector>

#include "libiop/algebra/fields/utils.hpp"
#include "libiop/algebra/utils.hpp"

namespace libiop {

/** Returns the vector of FRI localization parameters that is predicted to produce the smallest
 *  argument size for these parameters. This is calculated by brute forcing all options.
 *  The first localization parameter is fixed as 1, to mitigate verifier time,
 *  and to lower argument size since there could be many oracles.
 *
 *  TODO: Make this optimize for multiple FRI interactive repetitions, multiple oracles passed into FRI
 *  TODO: Make this take in max_tested_degree, and include accounting for non-power-of-two rates.
 *        This will have a significant impact on zk Aurora.
 *  TODO: Make it take in a locality vector. This isn't necessary for now, since the first parameter is fixed at 1.
 *        In a similar vein, have it keep track of which input oracles are zk. */
std::vector<size_t> brute_force_optimal_localization_parameters(
    size_t codeword_dim,
    size_t query_repetitions,
    std::size_t locality,
    std::size_t field_size,
    std::size_t RS_extra_dimensions);

} // namespace libiop

#include "libiop/protocols/ldt/fri/optimizer.tcc"

#endif // LIBIOP_PROTOCOLS_LDT_FRI_OPTIMIZER_HPP_
