#include "libiop/bcs/bcs_common.hpp"

#ifndef OPTIONS
#define OPTIONS

typedef struct options{
    std::size_t log_n_min;
    std::size_t log_n_max;
    std::size_t security_level;
    std::size_t field_size;
    std::size_t RS_extra_dimensions;
    std::size_t num_localization_steps;
    std::size_t num_oracles;
    std::size_t num_interactive_repetitions;
    std::size_t num_query_repetitions;
    std::size_t hash_enum_val;
    std::size_t localization_parameter;
    float height_width_ratio;
    bool heuristic_ldt_reducer_soundness;
    bool heuristic_fri_soundness;
    bool is_multiplicative;
    bool make_zk;
    bool optimize_localization;
    libiop::bcs_hash_type hash_enum;
} options;

#endif




