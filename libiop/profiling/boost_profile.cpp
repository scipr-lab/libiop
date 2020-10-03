#ifndef CPPDEBUG /* Ubuntu's Boost does not provide binaries compatible with libstdc++'s debug mode so we just reduce functionality here */
#include <boost/program_options.hpp>
#endif

#include "libiop/bcs/bcs_common.hpp"

namespace po = boost::program_options;
using namespace libiop;

// The hash_enum field must be initialized from the hash_enum_val field, as opposed to via the command line.
struct options{
    std::size_t log_n_min = 8;
    std::size_t log_n_max = 20;
    std::size_t security_level = 128;
    std::size_t field_size = 181;
    std::size_t hash_enum_val = (size_t) libiop::blake2b_type;
    bool heuristic_ldt_reducer_soundness = true;
    bool is_multiplicative = true;
    bool make_zk = false;
    libiop::bcs_hash_type hash_enum = blake2b_type;
};


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
		("is_multiplicative", po::value<bool>(&options.is_multiplicative)->default_value(options.is_multiplicative))
		("make_zk", po::value<bool>(&options.make_zk)->default_value(false))
		("hash_enum", po::value<std::size_t>(&options.hash_enum_val)->default_value((size_t) blake2b_type));


	return base;
}


