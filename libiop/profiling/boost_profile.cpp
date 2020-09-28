#ifndef CPPDEBUG /* Ubuntu's Boost does not provide binaries compatible with libstdc++'s debug mode so we just reduce functionality here */
#include <boost/program_options.hpp>
#endif


#include "snark_types.hpp"
namespace po = boost::program_options;
using namespace libiop;

po::options_description gen_options(options &options)
{
	po::options_description base("Usage");
	

	base.add_options()
	    ("help", "print this help message")
	    ("log_n_min", po::value<std::size_t>(&options.log_n_min)->default_value(options.log_n_min))
	    ("log_n_max", po::value<std::size_t>(&options.log_n_max)->default_value(options.log_n_max))
	    ("security_level", po::value<std::size_t>(&options.security_level)->default_value(options.security_level))
	    ("field_size", po::value<std::size_t>(&options.field_size)->default_value(options.field_size))
	    ("heuristic_ldt_reducer_soundness", po::value<bool>(&options.heuristic_ldt_reducer_soundness)->default_value(options.heuristic_ldt_reducer_soundness))
	    ("heuristic_fri_soundness", po::value<bool>(&options.heuristic_fri_soundness)->default_value(options.heuristic_fri_soundness)) /* Find a better solution for this in the future */
	    ("is_multiplicative", po::value<bool>(&options.is_multiplicative)->default_value(options.is_multiplicative))
		("make_zk", po::value<bool>(&options.make_zk)->default_value(false))
		("optimize_localization", po::value<bool>(&options.optimize_localization)->default_value(false))
		("hash_enum", po::value<std::size_t>(&options.hash_enum_val)->default_value((size_t) blake2b_type))
  		("localization_parameter", po::value<std::size_t>(&options.localization_parameter)->default_value(2), "Only used when num_localization_steps is 0")
        ("num_localization_steps", po::value<std::size_t>(&options.num_localization_steps)->default_value(0))
        ("num_oracles", po::value<std::size_t>(&options.num_oracles)->default_value(1))
        ("num_interactive_repetitions", po::value<std::size_t>(&options.num_interactive_repetitions)->default_value(1))
        ("num_query_repetitions", po::value<std::size_t>(&options.num_query_repetitions)->default_value(64));

	return base;
}


