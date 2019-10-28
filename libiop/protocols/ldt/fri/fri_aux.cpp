#include <cstdint>
#include "fri_aux.hpp"

namespace libiop {


// -------------------------------------------------
// Optimizer utils

/** Generate all possible localization vectors that begin with starting,
 *  and reduce up to max_reducable_dimensions dimensions.
 */
std::vector<std::vector<size_t>> localization_vector_generator(
    size_t max_reducable_dimensions,
    std::vector<size_t> starting)
{
    std::vector<std::vector<size_t>> options;
    options.push_back(starting);
    if (max_reducable_dimensions == 0)
    {
        return options;
    }
    for (size_t i = 1; i <= max_reducable_dimensions; ++i)
    {
        std::vector<size_t> new_starting = starting;
        new_starting.push_back(i);
        std::vector<std::vector<size_t>> new_options =
            localization_vector_generator(max_reducable_dimensions - i, new_starting);
        options.insert(options.end(), new_options.begin(), new_options.end());
    }
    return options;
}

/* return all partitions of this number (that is, all possible FRI localization parameter vectors
   for this codeword domain dimension) */
std::vector<std::vector<size_t>> all_localization_vectors(size_t dimension_to_reduce)
{
    /* Fix the start as 1 */
    std::vector<size_t> starting({1});
    return localization_vector_generator(dimension_to_reduce - 1, starting);
}

} // namespace libiop
