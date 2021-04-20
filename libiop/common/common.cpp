#include <libiop/common/common.hpp>

namespace libiop {

using std::size_t;

long double add_soundness_error_bits(const size_t bits1, const size_t bits2)
{
    return add_soundness_error_bits((long double)bits1, (long double)bits2);
}

long double add_soundness_error_bits(const long double bits1, const long double bits2)
{
    long double result = exp2l(-1.0 * bits1) + exp2l(-1.0 * bits2);
    return -1 * log2l(result);
}

} // namespace libiop
