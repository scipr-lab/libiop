#include "sodium/crypto_generichash_blake2b.h"
#include <cstring>
#include <stdexcept>

namespace libiop {

template<typename FieldT>
hash_digest blake2b_field_element_hash(const std::vector<FieldT> &data,
                                       const std::size_t digest_len_bytes)
{

    hash_digest result(digest_len_bytes, 'X');

    /* see https://download.libsodium.org/doc/hashing/generic_hashing.html */
    const int status = crypto_generichash_blake2b((unsigned char*)&result[0],
                                                  digest_len_bytes,
                                                  (result.empty() ? NULL : (unsigned char*)&data[0]),
                                                  sizeof(FieldT) * data.size(),
                                                  NULL, 0);
    if (status != 0)
    {
        throw std::runtime_error("Got non-zero status from crypto_generichash_blake2b. (Is digest_len_bytes correct?)");
    }


    return result;
}

template<typename FieldT>
std::vector<FieldT> blake2b_FieldT_randomness_extractor(const hash_digest &root,
                                                        const std::size_t index,
                                                        const std::size_t num_elements)
{
    const std::size_t root_plus_index_size = root.size() + sizeof(index);
    unsigned char* root_plus_index = (unsigned char*)(malloc(root_plus_index_size));
    memcpy(root_plus_index, &root[0], root.size());
    memcpy(root_plus_index + root.size(), &index, sizeof(index));

    std::vector<FieldT> result;
    result.reserve(num_elements);

    for (std::size_t i = 0; i < num_elements; ++i)
    {
        FieldT el;

        const int status = crypto_generichash_blake2b((unsigned char*)&el,
                                                      sizeof(el),
                                                      root_plus_index,
                                                      root_plus_index_size,
                                                      (unsigned char*)&i, sizeof(i));
        if (status != 0)
        {
            throw std::runtime_error("Got non-zero status from crypto_generichash_blake2b. (Is digest_len_bytes correct?)");
        }

        result.emplace_back(el);
    }

    free(root_plus_index);

    return result;
}

}
