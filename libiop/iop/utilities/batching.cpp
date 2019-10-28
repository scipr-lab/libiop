#include "batching.hpp"

namespace libiop {



std::vector<oracle_handle_ptr> virtual_oracle_handles_to_handle_ptrs(
    const std::vector<virtual_oracle_handle> handles)
{
    std::vector<oracle_handle_ptr> oracles;
    oracles.reserve(handles.size());
    for (size_t i = 0; i < handles.size(); i++)
    {
        oracles.emplace_back(std::make_shared<virtual_oracle_handle>(handles[i]));
    }
    return oracles;
}


} // namespace libiop
