//#include <libff/common/utils.hpp>
#include "libiop/protocols/ldt/fri/fri_ldt.hpp"

namespace libiop {


const char* FRI_soundness_type_to_string(FRI_soundness_type soundness_type)
{
    if (soundness_type == FRI_soundness_type::heuristic)
    {
        return "heuristic";
    } else if (soundness_type == FRI_soundness_type::proven)
    {
        return "proven";
    }
    return "Invalid soundness type";
}

} // namespace libiop
