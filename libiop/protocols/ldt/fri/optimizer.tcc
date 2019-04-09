
namespace libiop {
/* helper functions for estimating argument size */

/* return the expected total number of unique tree elements in the authentication paths
   of a certain number of queries */
std::size_t num_unique_queries_at_depth(std::size_t num_queries, std::size_t depth)
{
    float width = 1ull << depth;

    /* expected number of elements at this level untouched by paths (balls and bins) */
    float proportion_empty = std::pow((width - 1) / width, num_queries);
    return width * (1 - proportion_empty);
}

/* return the expected total number of unique tree elements in the authentication paths
   of a certain number of queries. */
std::size_t num_unique_elements_in_query_paths(std::size_t num_queries, std::size_t depth)
{
    float sum = 0.0;
    for (size_t d = 1; d < depth; ++d)
    {
        /* for each level of the tree */
        sum += num_unique_queries_at_depth(num_queries, d);
    }

    // std::cout << "predicted num_unique_elements_in_query_paths("
    //           << num_queries << ", " << depth << ") = " << ((std::size_t) std::round(sum)) << "\n";

    return ((std::size_t) std::round(sum));
}

std::size_t num_membership_proof_hashes(std::vector<std::size_t> loc_vector,
                                        std::size_t num_queries,
                                        std::size_t codeword_dim,
                                        std::size_t RS_extra_dimensions)
{
    std::size_t depth = codeword_dim + RS_extra_dimensions + 1;

    std::size_t total_hashes = 0;
    total_hashes += 2 * num_unique_elements_in_query_paths(num_queries, depth) - num_queries; /* auxiliary hashes */
    total_hashes += 2 * num_unique_queries_at_depth(num_queries, depth); /* randomness hashes */

    for (size_t i = 1; i < loc_vector.size(); ++i)
    {
        depth -= loc_vector[i];
        total_hashes += num_unique_elements_in_query_paths(num_queries, depth) - num_queries; /* auxiliary hashes */
        total_hashes += num_unique_queries_at_depth(num_queries, depth); /* auxiliary hashes */
    }

    return total_hashes;
}

std::size_t num_query_responses(std::vector<std::size_t> loc_vector,
                                std::size_t num_queries,
                                std::size_t codeword_dim)
{
    std::size_t depth = codeword_dim + 3;

    std::size_t total = 0;
    total += 14 * num_unique_queries_at_depth(num_queries, depth);

    for (size_t i = 1; i < loc_vector.size(); ++i)
    {
        depth -= loc_vector[i];
        total += num_unique_queries_at_depth(num_queries, depth) * (1ull << loc_vector[i]);
    }

    return total;
}

/* return the estimated argument size for this vector of FRI localization parameters */
std::size_t predictor(std::vector<std::size_t> loc_vector,
                      std::size_t codeword_dim,
                      std::size_t num_queries,
                      std::size_t locality,
                      std::size_t field_size,
                      std::size_t RS_extra_dimensions)
{
    std::size_t end = codeword_dim;
    for (size_t i = 0; i < loc_vector.size(); ++i)
    {
        end -= loc_vector[i];
    }
    std::size_t prover_messages_size = (field_size / 8) * ((1ull << (end + 2)) + 1);

    std::size_t query_responses_size = (field_size / 8) * num_query_responses(loc_vector, num_queries, codeword_dim);//queries;

    std::size_t MT_roots_size = 32 * (loc_vector.size() + 1);

    std::size_t hash_size = 32;
    std::size_t total_hashes = num_membership_proof_hashes(loc_vector, num_queries, codeword_dim, RS_extra_dimensions);
    std::size_t MT_multi_membership_proofs_size = hash_size * total_hashes;

    return prover_messages_size + query_responses_size + MT_roots_size + MT_multi_membership_proofs_size;
}

std::size_t averaged_predictor(std::vector<std::size_t> loc_vector,
                               std::size_t codeword_dim,
                               std::size_t num_queries,
                               std::size_t locality,
                               std::size_t field_size,
                               std::size_t RS_extra_dimensions,
                               std::size_t num_repetitions)
{
    std::size_t sum = 0;
    for (size_t i = 0; i < num_repetitions; ++i)
    {
        sum += predictor(loc_vector, codeword_dim, num_queries, locality, field_size, RS_extra_dimensions);
    }
    return sum / num_repetitions;
}

std::vector<std::vector<std::size_t>> all_options_helper(std::size_t left, std::vector<std::size_t> starting)
{
    std::vector<std::vector<std::size_t>> options;
    if (left == 0)
    {
        options.push_back(starting);
        return options;
    }
    for (size_t i = 1; i <= left; ++i)
    {
        std::vector<std::size_t> new_starting = starting;
        new_starting.push_back(i);
        std::vector<std::vector<std::size_t>> new_options = all_options_helper(left - i, new_starting);
        options.insert(options.end(), new_options.begin(), new_options.end());
    }
    return options;
}

/* return all partitions of this number (that is, all possible FRI localization parameter vectors
   for this codeword domain dimension) */
std::vector<std::vector<std::size_t>> all_options(std::size_t left)
{
    std::vector<std::size_t> starting;
    return all_options_helper(left, starting);
}

/* return the vector of FRI localization parameters that is predicted to produce the smallest
   argument size for these parameters */
std::vector<std::size_t> brute_force_optimal_localization_parameters(
    std::size_t codeword_dim,
    std::size_t num_queries,
    std::size_t locality,
    std::size_t field_size,
    std::size_t RS_extra_dimensions)
{
    std::vector<std::vector<std::size_t>> options;
    for (size_t total_loc = 2; total_loc < codeword_dim - 1; ++total_loc)
    {
        std::vector<std::vector<std::size_t>> new_options = all_options(total_loc);
        options.insert(options.end(), new_options.begin(), new_options.end());
    }
    for (size_t i = 0; i < options.size(); ++i)
    {
        options[i].insert(options[i].begin(), 1);
    }

    std::size_t num_repetitions = 5; /* increase accuracy by averaging predictions */

    std::size_t min = averaged_predictor(options[0], codeword_dim, num_queries, locality,
                                         field_size, num_repetitions, RS_extra_dimensions);
    std::vector<std::size_t> best = options[0];
    for (size_t i = 1; i < options.size(); ++i)
    {
        std::size_t current = averaged_predictor(options[i], codeword_dim, num_queries, locality,
                                                 field_size, num_repetitions, RS_extra_dimensions);
        if (current < min)
        {
            min = current;
            best = options[i];
        }
    }

    return best;
}

}