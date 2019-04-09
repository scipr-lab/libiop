/**@file
*****************************************************************************
This file defines oracle interfaces
*****************************************************************************
* @author     This file is part of libiop (see AUTHORS)
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/
#ifndef LIBIOP_IOP_ORACLES_HPP_
#define LIBIOP_IOP_ORACLES_HPP_

#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace libiop {

/* Oracles */
template<typename FieldT>
class oracle {
protected:
    std::vector<FieldT> evaluated_contents_;

public:
    oracle() = default;
    oracle(const std::vector<FieldT> &evaluated_contents) : evaluated_contents_(evaluated_contents) {}
    oracle(std::vector<FieldT> &&evaluated_contents) : evaluated_contents_(std::move(evaluated_contents)) {}

    const std::vector<FieldT>& evaluated_contents() const {
        return this->evaluated_contents_;
    }
};

/* Virtual oracle */
template<typename FieldT>
class virtual_oracle : public oracle<FieldT> {
public:
    virtual std::vector<FieldT> evaluated_contents(
        const std::vector<std::vector<FieldT> > &constituent_oracle_evaluations) const = 0;

    virtual FieldT evaluation_at_point(
        const std::size_t evaluation_position,
        const FieldT evaluation_point,
        const std::vector<FieldT> &constituent_oracle_evaluations) const = 0;
    /* Below the IOP interface defines

       virtual_oracle_handle register_virtual_oracle(
       const domain_handle &domain,
       const std::size_t degree,
       const std::vector<oracle_handle> &constituent_oracles);

       At each evaluation_at_point call the virtual oracle also gets
       access to the evaluations of its constituent oracles.

       This assumes that the virtual oracle is composed of constituent
       oracles evaluated *at the same point*. Such assumption holds for
       the protocols we have in mind.

       If the need arises, we can accommodate this by letting
       virtual_oracle get access to entire IOP (maybe via proxy class)
       and let virtual_oracle query constituent oracles at arbitrary
       evaluation points.
    */

    /* Users will likely want to subclass virtual_oracle and add
       methods for sharing data. For example, many different virtual
       oracles might want to share the same Lagrange polynomial
       evaluations. */

    ~virtual_oracle() = default;
};

} // namespace libiop

#endif // LIBIOP_IOP_ORACLES_HPP_
