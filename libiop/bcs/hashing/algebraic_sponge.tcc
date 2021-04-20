#include <libff/algebra/field_utils/field_utils.hpp>
#include "libiop/bcs/hashing/hashing.hpp"
#include <libff/algebra/field_utils/bigint.hpp>
#include <cstring>
#include <sstream>
#include <stdexcept>

namespace libiop {

template<typename FieldT>
algebraic_sponge<FieldT>::algebraic_sponge(const size_t rate, const size_t capacity) :
    state_(std::vector<FieldT>(rate + capacity, FieldT::zero())),
    rate_(rate),
    capacity_(capacity)
{
}

template<typename FieldT>
void algebraic_sponge<FieldT>::absorb(const std::vector<FieldT> &new_input)
{
    /** If we have already absorbed, 
     * we need to permute state in order to securely absorb more. 
     * We could optimize this to continue where the prior absorb left off. */
    if (this->currently_absorbing)
    {
        this->apply_permutation();
    }
    this->absorb_internal(new_input, 0);
    this->currently_absorbing = true;
}

template<typename FieldT>
void algebraic_sponge<FieldT>::absorb_internal(
    const std::vector<FieldT> &new_input,
    // index of the first element in new_input we read
    const size_t begin_index)
{
    /** Add new_input to the elements corresponding to the rate in state.
     *  We split this into two cases, one is when the remainder of rate is
     *  sufficient to absorb all of the new_input. 
     *  If so absorb the remainder of the input.
     * 
     *  The other case is to abosrb |rate| elements of the new input,
     *  apply the permutation, and then recurse.
    */
    if (new_input.size() - begin_index <= this->rate_)
    {
        for (size_t i = 0; i < new_input.size() - begin_index; i++)
        {
            this->state_[i] += new_input[i + begin_index];
        }
        return;
    }
    /** Absorb the next |rate| elements, permute, and keep absorbing */
    for (size_t i = 0; i < this->rate_; i++)
    {
        this->state_[i] += new_input[i + begin_index];
    }
    this->apply_permutation();
    // Tail recurse
    this->absorb_internal(new_input, begin_index + this->rate_);
}

template<typename FieldT>
std::vector<FieldT> algebraic_sponge<FieldT>::squeeze_vector(size_t num_elements)
{
    std::vector<FieldT> output(num_elements, FieldT::zero());
    if (this->currently_absorbing) {
        this->next_unsqueezed_elem_ = 0;
        this->currently_absorbing = false;
    }
    size_t index_in_output = 0;
    this->squeeze_internal(output, index_in_output);
    return output;
}

template<typename FieldT>
void algebraic_sponge<FieldT>::squeeze_internal(std::vector<FieldT> &output, size_t begin_index)
{
    // We need to squeeze. TODO: Have a better signal for this.
    if (this->next_unsqueezed_elem_ == 0)
    {
        this->apply_permutation();
    }
    size_t index_in_output = begin_index;
    while (this->next_unsqueezed_elem_ < this->rate_ && index_in_output < output.size())
    {
        output[index_in_output] = this->state_[this->next_unsqueezed_elem_];
        index_in_output++;
        this->next_unsqueezed_elem_++;
    }
    /* We've squeezed as much as we need. */
    if (index_in_output == output.size())
    {
        return;
    }
    /* Otherwise we need to squeeze more. Apply a permutation and keep squeezing. */
    this->next_unsqueezed_elem_ = 0;
    this->squeeze_internal(output, index_in_output);
}

template<typename FieldT>
void algebraic_sponge<FieldT>::initialize_element_of_state(
    const FieldT elem, const size_t index)
{
    this->state_[index] = elem;
}

/* multiplicative case */
template<typename FieldT>
FieldT string_to_field_elem(
    typename libff::enable_if<libff::is_multiplicative<FieldT>::value, FieldT>::type dummy_field_elem,
    const zk_salt_type &zk_salt)
{
    libff::bigint<FieldT::num_limbs> num(0ul);
    assert(zk_salt.length() == 8 * FieldT::num_limbs);
    // copy string into big endian, word by word
    // This is because the bigint data is stored little endian
    for (size_t i = 0; i < FieldT::num_limbs; i++)
    {
        std::memcpy(&num.data[FieldT::num_limbs - i - 1], &zk_salt[i*sizeof(size_t)], sizeof(size_t));
    }
    return FieldT(num);
}

/* additive case */
template<typename FieldT>
FieldT string_to_field_elem(
    typename libff::enable_if<libff::is_additive<FieldT>::value, FieldT>::type dummy_field_elem,
    const zk_salt_type &zk_salt)
{
    throw std::invalid_argument("zk hash for binary fields is not yet implemented");
}

template<typename FieldT, typename MT_root_type>
algebraic_hashchain<FieldT, MT_root_type>::algebraic_hashchain(
    std::shared_ptr<algebraic_sponge<FieldT>> sponge,
    size_t security_parameter) :
    sponge_(sponge),
    security_parameter_(security_parameter)
{
    this->sponge_->reset();
    // We want to return a single field element from this for efficiency.
    // small field hashes are not yet supported.
    assert(this->sponge_->capacity_ == 1);
    // assert(this->sponge_->achieved_security_parameter() >= (double) security_parameter);
}

template<typename FieldT, typename MT_root_type>
void algebraic_hashchain<FieldT, MT_root_type>::absorb(
    const MT_root_type new_input)
{
    this->absorb_internal(new_input);
}

template<typename FieldT, typename MT_root_type>
void algebraic_hashchain<FieldT, MT_root_type>::absorb_internal(
    const typename libff::enable_if<std::is_same<MT_root_type, binary_hash_digest>::value, MT_root_type>::type new_input)
{
    FieldT new_input_as_FieldT = string_to_field_elem<FieldT>(new_input);
    this->sponge_->absorb(new_input_as_FieldT);
}

template<typename FieldT, typename MT_root_type>
void algebraic_hashchain<FieldT, MT_root_type>::absorb_internal(
    const typename libff::enable_if<std::is_same<MT_root_type, FieldT>::value, MT_root_type>::type new_input)
{
    this->sponge_->absorb(std::vector<FieldT>({new_input}));
}

template<typename FieldT, typename MT_root_type>
void algebraic_hashchain<FieldT, MT_root_type>::absorb(
    const std::vector<FieldT> &new_input)
{
    this->sponge_->absorb(new_input);
}

template<typename FieldT, typename MT_root_type>
std::vector<FieldT> algebraic_hashchain<FieldT, MT_root_type>::squeeze(
    size_t num_elements)
{
    return this->sponge_->squeeze_vector(num_elements);
}

template<typename FieldT, typename MT_root_type>
std::vector<size_t> algebraic_hashchain<FieldT, MT_root_type>::squeeze_query_positions(
    size_t num_elements, size_t range_of_positions)
{
    /** TODO: Squeeze more efficiently, we can get more than one query point per elem */
    std::vector<FieldT> squeezed_elems = this->sponge_->squeeze_vector(num_elements);
    std::vector<size_t> positions;
    /* Get least significant word of field element */
    size_t word_index = 0;
    for (size_t i = 0; i < squeezed_elems.size(); i++)
    {
        positions.emplace_back(libff::get_word_of_field_elem<FieldT>(squeezed_elems[i], word_index) % range_of_positions);
    }
    return positions;
}


template<typename FieldT, typename MT_root_type>
MT_root_type algebraic_hashchain<FieldT, MT_root_type>::squeeze_root_type()
{
    return this->sponge_->squeeze_vector(1)[0];
}

template<typename FieldT>
algebraic_leafhash<FieldT>::algebraic_leafhash(
    std::shared_ptr<algebraic_sponge<FieldT>> sponge,
    size_t security_parameter) :
    sponge_(sponge->new_sponge())
{
    this->sponge_->reset();
    // We want to return a single field element from this for efficiency.
    assert(this->sponge_->capacity_ == 1);
    // assert(this->sponge_->achieved_security_parameter() >= (double) security_parameter);
}

template<typename FieldT>
FieldT algebraic_leafhash<FieldT>::hash(
    const std::vector<FieldT> &leaf)
{
    this->sponge_->absorb(leaf);
    FieldT result = this->sponge_->squeeze_vector(1)[0];
    this->sponge_->reset();
    return result;
}

template<typename FieldT>
FieldT algebraic_leafhash<FieldT>::zk_hash(
    const std::vector<FieldT> &leaf,
    const zk_salt_type &zk_salt)
{
    std::vector<FieldT> leaf_copy(leaf);
    FieldT salt = string_to_field_elem<FieldT>(FieldT::zero(), zk_salt);
    leaf_copy.emplace_back(salt);
    this->sponge_->absorb(leaf_copy);
    FieldT result = this->sponge_->squeeze_vector(1)[0];
    this->sponge_->reset();
    return result;
}

template<typename FieldT>
algebraic_two_to_one_hash<FieldT>::algebraic_two_to_one_hash(
    std::shared_ptr<algebraic_sponge<FieldT>> sponge,
    size_t security_parameter) :
    sponge_(sponge)
{
    this->sponge_->reset();
    // We want to return a single field element from this for efficiency.
    assert(this->sponge_->capacity_ == 1);
    // assert(this->sponge_->achieved_security_parameter() >= (double) security_parameter);
}

template<typename FieldT>
FieldT algebraic_two_to_one_hash<FieldT>::hash(
    const FieldT &left, const FieldT &right)
{
    this->sponge_->initialize_element_of_state(left, 0);
    this->sponge_->initialize_element_of_state(right, 1);
    FieldT result = this->sponge_->squeeze_vector(1)[0];
    this->sponge_->reset();
    return result;
}

/** The following is unused code for generating arks and MDS.
 *  Currently we hardcode these, as obtained from sage.  */

// def generate_round_constant(fn_name, field, idx):
//     """
//     Returns a field element based on the result of sha256.
//     The input to sha256 is the concatenation of the name of the hash function
//     and an index.
//     For example, the first element for MiMC will be computed using the value
//     of sha256('MiMC0').
//     """
//     from hashlib import sha256
//     val = int(sha256('%s%d' % (fn_name, idx)).hexdigest(), 16)
//     if field.is_prime_field():
//         return field(val)
//     else:
//         return int2field(field, val % field.order())
// template<typename FieldT>
// std::vector<std::vector<FieldT>> generate_round_constant(std::string name, size_t id)
// {
//     // TODO: constant should be deterministic off of name, id
//     return FieldT::random_element();
// }

// template<typename FieldT>
// std::vector<std::vector<FieldT>> generate_mds_matrix(std::string name, size_t state_size)
// {
//     // TODO: Make randomness deterministic from name
//     const size_t num_attempts = 100;
//     std::string name_x = name + "x";
//     std::string name_y = name + "y";
//     for (size_t attempt = 0; attempt < num_attempts; attempt++)
//     {
//         std::vector<FieldT> x_vec_;
//         std::vector<FieldT> y_vec_;
//         for (size_t i = 0; i < state_size; i++)
//         {
//             x_vec_.emplace_back(generate_round_constant(name_x, attempt * state_size + i));
//             y_vec_.emplace_back(generate_round_constant(name_x, attempt * state_size + i));
//             // TODO: Add check for all x,y values being distinct.
//             // (Negligible probability of occuring with fields of interest)
//         }
//         // Sanity check: check the determinant of the matrix.
//         FieldT x_prod = FieldT::one();
//         FieldT y_prod = FieldT::one();
//         FieldT xy_prod = FieldT::one();
//         for (size_t i = 0; i < state_size; i++)
//         {
//             for (size_t j = 0; j < i; j++)
//             {
//                 x_prod *= x_vec_[i] - x_vec_[j];
//                 y_prod *= y_vec_[i] - y_vec_[j];
//             }
//             for (size_t j = 0; j < state_size; j++)
//             {
//                 xy_prod *= x_vec_[i] - y_vec_[j];
//             }
//         }
//         FieldT det = ((state_size % 4 < 2) ? 1 :-1) * x_prod * y_prod * (xy_prod.inverse());
//         assert(det != FieldT::zero());

//         // if len(mds.characteristic_polynomial().roots()) == 0:
//         //     // There are no eigenvalues in the field.
//         //     return mds
//     }
// }
}
