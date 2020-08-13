/**@file
 *****************************************************************************
 Dummy algebraic hash function implementation
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_SNARK_COMMON_HASHING_DUMMY_ALGEBRAIC_SPONGE_HPP_
#define LIBIOP_SNARK_COMMON_HASHING_DUMMY_ALGEBRAIC_SPONGE_HPP_

#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include "libiop/bcs/hashing/algebraic_sponge.hpp"

namespace libiop {

/** This is an insecure algebraic sponge */
template<typename FieldT>
class dummy_algebraic_sponge : public algebraic_sponge<FieldTs>
{
    public:
        void apply_permutation()
        {
            for (size_t i = 0; i < this->state_.size(); i++)
            {
                this->state_[i] = (1 + this->state_[i]);
            }
        };
        void print_permutation_security_parameter();
        dummy_algebraic_sponge();

        /* Needed for C++ polymorphism*/
        virtual void reset();
        virtual std::shared_ptr<algebraic_sponge<FieldT>> new_sponge();
};

} // namespace libiop

#include "libiop/bcs/hashing/dummy_algebraic_hash.tcc"

#endif // LIBIOP_SNARK_COMMON_HASHING_DUMMY_ALGEBRAIC_HASH_HPP_
