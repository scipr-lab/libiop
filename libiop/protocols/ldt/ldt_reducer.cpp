#include "libiop/protocols/ldt/ldt_reducer.hpp"
#include <cassert>
#include <stdexcept>

namespace libiop {


const char* LDT_reducer_soundness_type_to_string(LDT_reducer_soundness_type soundness_type)
{
    if (soundness_type == LDT_reducer_soundness_type::proven)
    {
        return "proven";
    } else if (soundness_type == LDT_reducer_soundness_type::optimistic_heuristic)
    {
        return "heuristic";
    }
    return "Invalid soundness type";
}


} // libiop
