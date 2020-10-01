#include "libiop/bcs/bcs_common.hpp"

#ifndef OPTIONS
#define OPTIONS

typedef struct options{
    std::size_t log_n_min;
    std::size_t log_n_max;
    std::size_t security_level;
    std::size_t field_size;
    std::size_t hash_enum_val;
    bool heuristic_ldt_reducer_soundness;
    bool is_multiplicative;
    bool make_zk;
    libiop::bcs_hash_type hash_enum;
} options;

#endif




