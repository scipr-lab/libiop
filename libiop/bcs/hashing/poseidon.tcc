#include <libff/algebra/field_utils/field_utils.hpp>
#include "libiop/bcs/hashing/algebraic_sponge.hpp"
#include <libff/common/profiling.hpp>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <cmath>

namespace libiop {

using libff::bigint;

template<typename FieldT>
poseidon_params<FieldT>::poseidon_params(
    size_t full_rounds, size_t partial_rounds, size_t alpha,
    size_t rate,
    std::vector<std::vector<FieldT>> &ark_matrix,
    bool supported_near_mds,
    std::vector<std::vector<FieldT>> &mds_matrix) :
    full_rounds_(full_rounds),
    partial_rounds_(partial_rounds),
    alpha_(alpha),
    state_size_(mds_matrix.size()),
    rate_(rate),
    capacity_(state_size_ - rate_),
    ark_matrix_(ark_matrix),
    supported_near_mds_(supported_near_mds),
    mds_matrix_(mds_matrix)
{
    if(ark_matrix.size() != full_rounds + partial_rounds)
    {
        throw std::invalid_argument("ark_matrix is of wrong dimension");
    }

    if(mds_matrix.size() != this->state_size_)
    {
        throw std::invalid_argument("mds_matrix is of wrong dimension");
    }
}

/** This follows Poseidon(https://eprint.iacr.org/2019/458.pdf) section 4.1
 *  for the soundness requirements.
 * 
 *  We implement the prime field case, with alpha = 3 or alpha = 5.
 */
template<typename FieldT>
double poseidon_params<FieldT>::achieved_soundness() const
{
    /** Rf >= 6 for statistical attacks
     * However the paper adds 25% number of rounds to full rounds for a given security level
     */
 
    if (this->full_rounds_ < 8)
    {
        return 0.0;
    }
    /* Switch to the variables of the paper */
    const size_t t = this->state_size_;
    /** Because the expressions are comparing if expr(num_rounds) > min(lambda, n),
     *  we set n to be ceil(log_2(|F|)), since that makes it larger.
     */
    const size_t n = libff::log_of_field_size_helper<FieldT>(FieldT::zero());
    /* build in the 25% extra full rounds, 7.5% partial rounds. */
    const double effective_full_rounds = ceil(.75 * this->full_rounds_);
    const double effective_partial_rounds = ceil(.925 * this->partial_rounds_);
    const double num_rounds = effective_full_rounds + effective_partial_rounds;

    /** Start at the maximum security level (half the bits in capacity), 
     *  and then lower it down with each inequality. 
     *  We use floor(log_2(|F|)) because this is the more conservative estimate. */
    double current_security_level = ((double) this->capacity_) * 
        ((double) libff::soundness_log_of_field_size_helper<FieldT>(FieldT::zero())) / 2.0;

    /** For the case of Poseidon we implement,
     *  our parameters must not satisfy equations (1) and (2) of section 3. */
    /** (1) yields: num_rounds > log_{alpha}(2) min(lambda, n) + log_2(t) 
     *  So 
     *      min(lambda, n) < (num_rounds - log_2(t)) / log_{alpha}(2)    
     *  So we solve for lambda in this equation. */
    double rhs_1 = num_rounds - (std::log(t)/std::log(2));
    rhs_1 = rhs_1 / (std::log(2) / std::log(this->alpha_));
    if (rhs_1 <= n)
    {
        current_security_level = std::min(rhs_1, current_security_level);
    }

    /** (2) yields two conditions:
     *  (2.1) num_rounds > [alpha_dependent_constant] min(lambda, n) 
     *  or equivalently
     *      min(lambda, n) < (num_rounds) / [alpha_dependent_constant]
     *  (2.2) (t - 1)*Rf + Rp > [alpha_dependent_constant] min(lambda, n)
     *  or equivalently
     *      min(lambda, n) < (t - 1)*Rf + Rp / [alpha_dependent_constant]
     *  So we solve for lambda in these equations. */

    // These are just pulled from the paper for 3, 5. It can't be larger than 1, hence that default.
    double alpha_dependent_constant_for_2_1 = 1;
    double alpha_dependent_constant_for_2_2 = 1;
    if (this->alpha_ == 3)
    {
        alpha_dependent_constant_for_2_1 = .32;
        alpha_dependent_constant_for_2_2 = .18;
    } 
    else if (this->alpha_ >= 5)
    {
        /* constants for alpha = 5, but these should decrease as alpha increases */
        alpha_dependent_constant_for_2_1 = .21;
        alpha_dependent_constant_for_2_2 = .14;
    }
    double rhs_2_1 = num_rounds / alpha_dependent_constant_for_2_1;
    double rhs_2_2 = ((t - 1) * effective_full_rounds + effective_partial_rounds) / alpha_dependent_constant_for_2_2;
    if (rhs_2_1 <= n)
    {
        current_security_level = std::min(rhs_2_1, current_security_level);
    }
    if (rhs_2_2 <= n)
    {
        current_security_level = std::min(rhs_2_2, current_security_level);
    }
    return current_security_level;
}

template<typename FieldT>
void poseidon_params<FieldT>::print() const
{
    libff::print_indent(); printf("\nPoseidon parameters\n");
    libff::print_indent(); printf("* State size = %lu\n", this->state_size_);
    libff::print_indent(); printf("* Rate = %lu\n", this->rate_);
    libff::print_indent(); printf("* Capacity = %lu\n", this->capacity_);
    libff::print_indent(); printf("* Alpha = %lu\n", this->alpha_);
    libff::print_indent(); printf("* Full rounds = %lu\n", this->full_rounds_);
    libff::print_indent(); printf("* Partial rounds = %lu\n", this->partial_rounds_);
    libff::print_indent(); printf("* Achieved security = %f\n", this->achieved_soundness());
}

template<typename FieldT>
poseidon<FieldT>::poseidon(
    const poseidon_params<FieldT> params):
    algebraic_sponge<FieldT>(params.rate_, params.capacity_),
    params_(params),
    zero_singleton_(FieldT::zero()),
    a_(params.mds_matrix_[1][1])
{   
    this->scratch_state_ = std::vector<FieldT>(params.state_size_, this->zero_singleton_);
}

template<typename FieldT>
double poseidon<FieldT>::achieved_security_parameter() const
{
    return this->params_.achieved_soundness();
}

template<typename FieldT>
void poseidon<FieldT>::print() const
{
    return this->params_.print();
}

template<typename FieldT>
FieldT poseidon<FieldT>::raise_to_alpha(const FieldT x) const
{

    if (this->params_.alpha_ == 17)
    {
        /* x^2 */
        FieldT intermediate = x * x;
        /* x^4 */
        intermediate *= intermediate;
        /* x^8 */
        intermediate *= intermediate;
        /* x^16 */
        intermediate *= intermediate;
        /* x^17 */
        intermediate *= x;
        return intermediate;
    }
    else if (this->params_.alpha_ == 5)
    {
        FieldT intermediate = x * x;
        intermediate *= intermediate;
        return x * intermediate;
    }
    else if (this->params_.alpha_ == 3)
    {
        FieldT intermediate = x * x;
        return x * intermediate;
    }
    else
    {
        return libff::power(x, this->params_.alpha_);
    }
    
}

template<typename FieldT>
void poseidon<FieldT>::apply_mix_layer()
{
    if (this->params_.supported_near_mds_ && this->params_.state_size_ == 3)
    {
        // assumes the near-MDS matrix is of form:
        // [[1, 0, 1],
        //  [1, 1, 0],
        //  [0, 1, 1]]
        // Let our state vector be called v = [x, y, z], and Mv = [x', y', z']
        const FieldT x_copy = this->state_[0];
        // Set x' = x + z
        this->state_[0] += this->state_[2];
        // Set z' = y + z
        this->state_[2] += this->state_[1];
        // Set y' = x + y
        this->state_[1] += x_copy;
    }
    else if (this->params_.supported_near_mds_ && this->params_.state_size_ == 4)
    {
        /**  We use a near MDS matrix of the form:
         *   [[0, 1, 1, 1],
         *    [1, 0, 1, 1],
         *    [1, 1, 0, 1],
         *    [1, 1, 1, 0]]
         */
        const FieldT complete_sum = (this->state_[0] + this->state_[1]) + (this->state_[2] + this->state_[3]);
        for (size_t i = 0; i < 4; i++)
        {
            this->state_[i] = complete_sum - this->state_[i];
        }
    }
    else
    {
        for (size_t row = 0; row < this->params_.state_size_; row++)
        {
            this->scratch_state_[row] = this->zero_singleton_;
            for (size_t col = 0; col < this->params_.state_size_; col++)
            {
                this->scratch_state_[row] += this->params_.mds_matrix_[row][col] * this->state_[col];
            }
        }
        std::swap(this->state_, this->scratch_state_);   
    }
}

template<typename FieldT>
void poseidon<FieldT>::apply_full_round(size_t round_id)
{
    const size_t state_size = this->params_.state_size_;
    // Add-Round Key & S-Box
    for (size_t i = 0; i < state_size; i++)
    {
        this->state_[i] += this->params_.ark_matrix_[round_id][i];
        this->state_[i] = this->raise_to_alpha(this->state_[i]);
    }

    // MixLayer
    // state = params.mds * state
    this->apply_mix_layer();
}

template<typename FieldT>
void poseidon<FieldT>::apply_partial_round(size_t round_id)
{
    const size_t state_size = this->params_.state_size_;
    // Add-Round Key & S-Box
    for (size_t i = 0; i < state_size; i++)
    {
        this->state_[i] += this->params_.ark_matrix_[round_id][i];
    }
    this->state_[state_size - 1] = this->raise_to_alpha(this->state_[state_size - 1]);

    // MixLayer
    // state = params.mds * state
    this->apply_mix_layer();
}

template<typename FieldT>
void poseidon<FieldT>::apply_permutation()
{
    size_t round = 0;
    bool full_round = true;
    for (size_t i = 0; i < this->params_.full_rounds_ / 2; i++)
    {
        this->apply_full_round(round);
        round++;
    }

    for (size_t i = 0; i < this->params_.partial_rounds_; i++)
    {
        this->apply_partial_round(round);
        round++;
    }

    for (size_t i = 0; i < this->params_.full_rounds_ / 2; i++)
    {
        this->apply_full_round(round);
        round++;
    }

    assert(round == this->params_.full_rounds_ + this->params_.partial_rounds_);
}

template<typename FieldT>
void poseidon<FieldT>::reset()
{
    // TODO: Replace with memset 0
    for (size_t i = 0; i < this->state_.size(); i++)
    {
        this->state_[i] = this->zero_singleton_;
    }
    this->next_unsqueezed_elem_ = 0;
    this->currently_absorbing = false;
}

template<typename FieldT>
poseidon_params<FieldT> default_128_bit_altbn_poseidon_params()
{
    const size_t alpha = 5;
    const size_t rate = 2;
    const size_t capacity = 1;
    /** Copied from the paper's table for 128 bit security */
    const size_t full_rounds = 8;
    const size_t partial_rounds = 56;
    const bool supported_near_mds = false;
    /** Magic constants obtained from script, but with an optimal near-MDS matrix */
    std::vector<std::vector<FieldT>> mds_matrix;
    mds_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("9478896780421655835758496955063136571251874317427585180076394551808670301829")),FieldT(bigint<FieldT::num_limbs>("1410220424381727336803825453763847584610565307685015130563813219659976870089")),FieldT(bigint<FieldT::num_limbs>("12324248147325396388933912754817224521085038231095815415485781874375379288849"))}));
    mds_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("5869197693688547188262203345939784760013629955870738354032535473827837048029")),FieldT(bigint<FieldT::num_limbs>("7027675418691353855077049716619550622043312043660992344940177187528247727783")),FieldT(bigint<FieldT::num_limbs>("12525656923125347519081182951439180216858859245949104467678704676398049957654"))}));
    mds_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("2393593257638453164081539737606611596909105394156134386135868506931280124380")),FieldT(bigint<FieldT::num_limbs>("21284282509779560826339329447865344953337633312148348516557075030360788076689")),FieldT(bigint<FieldT::num_limbs>("9426009547297688316907727916185688178981799243406990694957955230529774886223"))}));
    std::vector<std::vector<FieldT>> ark_matrix;
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("9478896780421655835758496955063136571251874317427585180076394551808670301829")),FieldT(bigint<FieldT::num_limbs>("1410220424381727336803825453763847584610565307685015130563813219659976870089")),FieldT(bigint<FieldT::num_limbs>("12324248147325396388933912754817224521085038231095815415485781874375379288849"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("5869197693688547188262203345939784760013629955870738354032535473827837048029")),FieldT(bigint<FieldT::num_limbs>("7027675418691353855077049716619550622043312043660992344940177187528247727783")),FieldT(bigint<FieldT::num_limbs>("12525656923125347519081182951439180216858859245949104467678704676398049957654"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("2393593257638453164081539737606611596909105394156134386135868506931280124380")),FieldT(bigint<FieldT::num_limbs>("21284282509779560826339329447865344953337633312148348516557075030360788076689")),FieldT(bigint<FieldT::num_limbs>("9426009547297688316907727916185688178981799243406990694957955230529774886223"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("5340930120720868177469579986808462005013697381998009281661327587975132166755")),FieldT(bigint<FieldT::num_limbs>("13224952063922250960936823741448973692264041750100990569445192064567307041002")),FieldT(bigint<FieldT::num_limbs>("5263772922985715307758718731861278699232625525745635678504665316187832057553"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("12905140589386545724352113723305099554526316070018892915579084990225436501424")),FieldT(bigint<FieldT::num_limbs>("3682692866591423277196501877256311982914914533289815428970072263880360882202")),FieldT(bigint<FieldT::num_limbs>("19681976272543335942352939522328215645129363120562038296975370569202780487598"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("5636115553781577891149626756201577064573186936824720926267940879716772984728")),FieldT(bigint<FieldT::num_limbs>("9501050736957980494328252533770324735114766672253853282051029963140075785396")),FieldT(bigint<FieldT::num_limbs>("2809392708032113981798687947163092027958611686144429680366467696224014505992"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("18433043696013996573551852847056868761017170818820490351056924728720017242180")),FieldT(bigint<FieldT::num_limbs>("1600424531609887868281118752288673305222025191763201214001133841689879221076")),FieldT(bigint<FieldT::num_limbs>("4077863666335789263839414578443702921867998881654209770993100224779179660280"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("10750183821931976144366649760909312224094512474274826686974526305203678408743")),FieldT(bigint<FieldT::num_limbs>("5876585841304782856135279046524906005004905983316552629403091395701737015709")),FieldT(bigint<FieldT::num_limbs>("13484299981373196201166722380389594773562113262309564134825386266765751213853"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("17382139184029132729706972098151128411278461930818849113274328379445169530719")),FieldT(bigint<FieldT::num_limbs>("20539300163698134245746932972121993866388520784246731402041866252259697791654")),FieldT(bigint<FieldT::num_limbs>("149101987103211771991327927827692640556911620408176100290586418839323044234"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3772300053282831651551351000101118094165364582047053942163129913249479587871")),FieldT(bigint<FieldT::num_limbs>("1859494671365748569037492975272316924127197175139843386363551067183747450207")),FieldT(bigint<FieldT::num_limbs>("6056775412522970299341516426839343188000696076848171109448990325789072743616"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("13535861576199801040709157556664030757939966797019046516538528720719863222691")),FieldT(bigint<FieldT::num_limbs>("3166287940256215995277151981354337176516077219230228956292184356796876826882")),FieldT(bigint<FieldT::num_limbs>("3878105211417696553129343540655091450996375987051865710523878345663272335218"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3234972450165117119793849127765475568944145932922109597427102281521349833458")),FieldT(bigint<FieldT::num_limbs>("4245107901241859301876588161430872878162557070919886440605528540426123750702")),FieldT(bigint<FieldT::num_limbs>("14797507122636944484020484450153618519329103538375805997650508264647579279513"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3893725073760673244819994221888005992135922325903832357013427303988853516024")),FieldT(bigint<FieldT::num_limbs>("21641836396029226240087625131527365621781742784615208902930655613239471409203")),FieldT(bigint<FieldT::num_limbs>("4622082908476410083286670201138165773322781640914243047922441301693321472984"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("14738633807199650048753490173004870343158648561341211428780666160270584694255")),FieldT(bigint<FieldT::num_limbs>("2635090520059500019661864086615522409798872905401305311748231832709078452746")),FieldT(bigint<FieldT::num_limbs>("19070766579582338321241892986615538320421651429118757507174186491084617237586"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("12622420533971517050761060317049369208980632120901481436392835424625664738526")),FieldT(bigint<FieldT::num_limbs>("4395637216713203985567958440367812800809784906642242330796693491855644277207")),FieldT(bigint<FieldT::num_limbs>("13856237567677889405904897420967317137820909836352033096836527506967315017500"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("2152570847472117965131784005129148028733701170858744625211808968788882229984")),FieldT(bigint<FieldT::num_limbs>("6585203416839617436007268534508514569040432229287367393560615429950244309612")),FieldT(bigint<FieldT::num_limbs>("2153122337593625580331500314713439203221416612327349850130027435376816262006"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("7340485916200743279276570085958556798507770452421357119145466906520506506342")),FieldT(bigint<FieldT::num_limbs>("12717879727828017519339312786933302720905962296193775803009326830415523871745")),FieldT(bigint<FieldT::num_limbs>("5392903649799167854181087360481925061021040403603926349022734894553054536405"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("7221669722700687417346373353960536661883467014204005276831020252277657076044")),FieldT(bigint<FieldT::num_limbs>("8259126917996748375739426565773281408349947402369855975457055235880500335093")),FieldT(bigint<FieldT::num_limbs>("9272385735015968356236075957906198733226196415690072035874639311675477515202"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("10999027991078055598627757097261950281899485771669414759870674222957875237568")),FieldT(bigint<FieldT::num_limbs>("15453393396765207016379045014101989306173462885430532298601655955681532648226")),FieldT(bigint<FieldT::num_limbs>("5478929644476681096437469958231489102974161353940993351588559414552523375472"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("6864274099016679903139678736335228538241825704814597078997020342617052506183")),FieldT(bigint<FieldT::num_limbs>("12133526413093116990739357861671284889661106676453313677855438696597541491864")),FieldT(bigint<FieldT::num_limbs>("4363234898901124667709814170397096827222883770682185860994495523839008586252"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("16799465577487943696587954846666404704275729737273450161871875150400464433797")),FieldT(bigint<FieldT::num_limbs>("3466902930973160737502426090330438125630820207992414876720169645462530526357")),FieldT(bigint<FieldT::num_limbs>("10062441698891350053170325824989022858836994651376301483266809451301259521913"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("5849282602749563270643968237860161465694876981255295041960826011116890638924")),FieldT(bigint<FieldT::num_limbs>("18460093993858702487671589299005229942046272739124591066186726570539410116617")),FieldT(bigint<FieldT::num_limbs>("9812100862165422922235757591915383485338044715409891361026651619010947646011"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3387849124775103843519196664933515074848119722071551419682472701704619249120")),FieldT(bigint<FieldT::num_limbs>("5283840871671971215904992681385681067319154145921438770232973796570506340281")),FieldT(bigint<FieldT::num_limbs>("14450974197863079729258614455552607708855872944526185987072755641686663205867"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("12613293459867195704822743599193025685229122593088639435739984309110321350551")),FieldT(bigint<FieldT::num_limbs>("6228273556621778927381918766322387348845347649737780310185999880647567569148")),FieldT(bigint<FieldT::num_limbs>("7482296435079443913598332362891173417094991594500715575107878549173583070413"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("18655449861670697203232484600163743308157596453845950955559776266093852537258")),FieldT(bigint<FieldT::num_limbs>("19948920146235041970991269588233091409704340607794045065548049409652881283328")),FieldT(bigint<FieldT::num_limbs>("13866078374565054775555309394949653928903776100036987352339975076159400168494"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("19398653685274645718325650121748668221118186023117741800737442235635318532994")),FieldT(bigint<FieldT::num_limbs>("4234154881267169381851681265196336178292466185695662916289548353755778788440")),FieldT(bigint<FieldT::num_limbs>("12763628380946395634691260884409562631856128057257959813602172954351304541746"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("7882453112990894293341171586279209575183467873317150236705310601775347127762")),FieldT(bigint<FieldT::num_limbs>("5669812778237054435250482766817044415794242063465169363632154286378940417646")),FieldT(bigint<FieldT::num_limbs>("16998738906020038479274018881471127087312245548341958049900081105113388112420"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3923902724726826782251513956816550869721438812970437824859252798290604500141")),FieldT(bigint<FieldT::num_limbs>("8649850619802776810849631749100283821801281306919958924112424995025830909252")),FieldT(bigint<FieldT::num_limbs>("11095642206650177249637693917287763476332497377393343056089442602164577098005"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("6935839211798937659784055008131602708847374430164859822530563797964932598700")),FieldT(bigint<FieldT::num_limbs>("7009671085960032501857416946339379996865118520008277046653124221544059312084")),FieldT(bigint<FieldT::num_limbs>("14361753917538892938870644779277430374939140280641641154553910654644462796654"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("6296738827713642491839335218022320853584196754765009910619998033694434027436")),FieldT(bigint<FieldT::num_limbs>("13849351053619304861036345979638534258290466678610892122310972291285921828452")),FieldT(bigint<FieldT::num_limbs>("434708832289952835651719825370636597763362139118091644948171210201038442144"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("16633750393567936099837698146248798150044883935695159627422586429892098538881")),FieldT(bigint<FieldT::num_limbs>("12944939557587269500508410478785174192748264930676627398550886896505925728421")),FieldT(bigint<FieldT::num_limbs>("13132297714437965464312509267711212830308064898189789451541658159340762509645"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3197382106307730326149017386920960267079843887376371149099833465681078850285")),FieldT(bigint<FieldT::num_limbs>("1219439673853113792340300173186247996249367102884530407862469123523013083971")),FieldT(bigint<FieldT::num_limbs>("3493891993991676033939225547105305872211028239751045376877382816726002847983"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("17474961424148900675164871904345354895260557993970869987490270849177572737815")),FieldT(bigint<FieldT::num_limbs>("14496326112831768456074139601688618143496262542471380389977686658437504436331")),FieldT(bigint<FieldT::num_limbs>("2924472580096769678506212811457662807142794313402961128576445038927398235897"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("4628296006426596599826873705217702584581936573072175641058168144816722698331")),FieldT(bigint<FieldT::num_limbs>("21191637522268746884323101636631937283436518241594045635071026927358145697662")),FieldT(bigint<FieldT::num_limbs>("16951212238971640283544926666565087199118390400059790490897089817025688673127"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("19613695336435411200907478310503966803576648245805018042761984388590288078910")),FieldT(bigint<FieldT::num_limbs>("19408817842355340096520725353160494939342325645253279486424056603334799168015")),FieldT(bigint<FieldT::num_limbs>("21454045045501902703155952158575095010854214688097850310899813261125869452799"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("7770328480231095569114093553841085793308707788942057894109603074902652929530")),FieldT(bigint<FieldT::num_limbs>("16464571997310094273270381226660568195148193554716113613093103468413654931642")),FieldT(bigint<FieldT::num_limbs>("17470702407108506528534764015553049093186219898758900659217736458688524875937"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("18550730212998825286534234924565339469725380540133305684933015562293032312245")),FieldT(bigint<FieldT::num_limbs>("2896017217286658654468296502214232988965841950467453595108246966331694256153")),FieldT(bigint<FieldT::num_limbs>("14675299739240143232464986549869467617250208852063994519435190317578889428919"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("18248678941898574458551739591299676471465479218953857934957966658325546357470")),FieldT(bigint<FieldT::num_limbs>("16024355190461954313137023461204350082170056518459846158694838747631538641098")),FieldT(bigint<FieldT::num_limbs>("7850865488269095796873132539790837027396511819869707882135427785654995088674"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3624709825280601739427579377772539527675482368902561427011111819459495093367")),FieldT(bigint<FieldT::num_limbs>("15151686747202278393509880888871669223363163423030702442491136636685867520291")),FieldT(bigint<FieldT::num_limbs>("13928718275408319129936199508884886190840388527360771842040395542259777472217"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("16438518367400433996754402147971541957362882260288877185296680838641299101277")),FieldT(bigint<FieldT::num_limbs>("14251644764165619982007152928697135211685596697488988945688571805257723263990")),FieldT(bigint<FieldT::num_limbs>("1558415498960552213241704009433360128041672577274390114589014204605400783336"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("10749063905018409753971835881025233122402531635641702583180695671822925094888")),FieldT(bigint<FieldT::num_limbs>("16354988262062662857299847129894034999364829278876526229109898088936579282653")),FieldT(bigint<FieldT::num_limbs>("18674506161104582501429570542497403127395303370792353981110261037140868411363"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("257179637705917619258734318345044380508175660880825307674152181190752217180")),FieldT(bigint<FieldT::num_limbs>("9940762366329880573269658025932923420224913533497554139322991519233208344887")),FieldT(bigint<FieldT::num_limbs>("1293239921425673430660897025143433077974838969258268884994339615096356996604"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3854755057793743998382749043935645076321648392889271506938052408708736403452")),FieldT(bigint<FieldT::num_limbs>("12412189818653171938404326534798862660266918801480759205769134043988251714987")),FieldT(bigint<FieldT::num_limbs>("14244337044908816699730279807127306408590045371642046203702435001248807053599"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("2664474558454245875259863292210158373587436806365583470716142566181853525483")),FieldT(bigint<FieldT::num_limbs>("13804808032959446920445446043982766933818740408606455713901597406432607762297")),FieldT(bigint<FieldT::num_limbs>("15382392877686271824754965441189716482120736348466561135774995773639134591702"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3837565647574361068761459036361038720977149057378981325662868867989693614839")),FieldT(bigint<FieldT::num_limbs>("15225462598913514623919729309237607260922732019503241269148802966711589239234")),FieldT(bigint<FieldT::num_limbs>("14538071498526966569560493408197490744750860055409168184834291065931235880474"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("5410684378713263186548337983252960352059120561867029137014627845977008289113")),FieldT(bigint<FieldT::num_limbs>("14667331447809547522064866822517742824275266629587157982725082691761010479709")),FieldT(bigint<FieldT::num_limbs>("2791759385703122291826173973326063535105256500818872868576708473427706638753"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3642448711195731276443459026400297797336444189377301111331173182906175572703")),FieldT(bigint<FieldT::num_limbs>("1762188042455633427137702520675816545396284185254002959309669405982213803405")),FieldT(bigint<FieldT::num_limbs>("2009136114070652405039722031642909652809813532808095287569356370144739527789"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("19791483249730874013718425513633316012832876840094804012613683750613457029295")),FieldT(bigint<FieldT::num_limbs>("11215847714185086315143282238988904989100073569422371901331104090526168332344")),FieldT(bigint<FieldT::num_limbs>("6136518117972887266318121761768078738116787014241169482775768363309364395343"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("18745129065505118158068407282927084003185002203481097964691803189568879488516")),FieldT(bigint<FieldT::num_limbs>("4425941001346546123680947386796057642506098630852783901377461933889171987192")),FieldT(bigint<FieldT::num_limbs>("18705881328407019529873437562261994664374735037477656748280619515075691198549"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("13047230314802854075744753527622973042529217055807151931674300909891623306887")),FieldT(bigint<FieldT::num_limbs>("6366778494392099924891219059630877678312150007619466052601251344235295985587")),FieldT(bigint<FieldT::num_limbs>("11984731463267744847120482801387059157065563300719282695907722257441888163935"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("6190354345246662069799625015516451435683819097684409839182088413597133300258")),FieldT(bigint<FieldT::num_limbs>("3453567349071653221618278943759136615746795332310201741158620901333366005083")),FieldT(bigint<FieldT::num_limbs>("3501451739665002417829002105335662012377395176791941899829101195688761736951"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("667342272997903647799602961987409419211738966452495719033033886694078041692")),FieldT(bigint<FieldT::num_limbs>("10475752001342638156955208359462181450676075289835626306234052334380574279894")),FieldT(bigint<FieldT::num_limbs>("1122105894868134874992852500882246664787309039693855267226551730451529279581"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("18579560865980555892476208730138769312184431461667680001251919743495048225398")),FieldT(bigint<FieldT::num_limbs>("11251487055721041066586341648434378410476221594972213638601515193325859147370")),FieldT(bigint<FieldT::num_limbs>("19196450369181962257928855170898976494942474416888497577021405056315702871427"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("18411083633425943022331023975608648608310348081260352933038968953840205960778")),FieldT(bigint<FieldT::num_limbs>("6327119822682179963667549434417325444351110952045620002522685857551651083472")),FieldT(bigint<FieldT::num_limbs>("17501079706299897104368639554646687210181409367601472897222937217901323545740"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("20969954387293127728287851205314979742002227914237776535291202019327184327937")),FieldT(bigint<FieldT::num_limbs>("8956086407974197848637272350322833148022602624785508354504602767996995213788")),FieldT(bigint<FieldT::num_limbs>("4081851985837245579578469133691002784866447850417410344029362852908793369007"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("4351604359086718117117597527743967124832774640980512524849354289224154141919")),FieldT(bigint<FieldT::num_limbs>("17615239559821402150637526649568971554876087064285249853519200056777204351029")),FieldT(bigint<FieldT::num_limbs>("4662985885051767286461661095506125353501801003851774063074483362173101455046"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("20179161646777445671702557506811137386924574568973382220913336875908534032741")),FieldT(bigint<FieldT::num_limbs>("16993030147580473888759301868545691152777330918124464719984004299784020276842")),FieldT(bigint<FieldT::num_limbs>("575664225980134190122457336356492934406612594003355445554139366898663257982"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("1900258475289806548211573591054356955406483174574792045064852691659255540421")),FieldT(bigint<FieldT::num_limbs>("19116096370228927568513878662761902411232415277576516248961941617882154589328")),FieldT(bigint<FieldT::num_limbs>("1951789632468349367267733649384193251827435198733250007703989540830646760455"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("411529621724849932999694270803131456243889635467661223241617477462914950626")),FieldT(bigint<FieldT::num_limbs>("4452093714558466154916571736623556422818297967689857034863500354807778364762")),FieldT(bigint<FieldT::num_limbs>("1144141104429089340350948910140160725026696285203534366803195498953182773627"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("17259840274319308134929561923186467767651283864853013159164761466977425646479")),FieldT(bigint<FieldT::num_limbs>("19864642394153108012585079380461314200434481284784864449963273042089975440262")),FieldT(bigint<FieldT::num_limbs>("16741663861181539595574639769219712417857243108953319368901488900165264615337"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("1661082462555315051331851196436613917688932192091465712488615203496731614940")),FieldT(bigint<FieldT::num_limbs>("3371894596850689293111972123461178477251217907955263260683014622910817304069")),FieldT(bigint<FieldT::num_limbs>("13889595442177754598892537082881169375497819074782697987852973860359070672134"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("291869714509212179917854069060736414045647593399252154880549638245744508325")),FieldT(bigint<FieldT::num_limbs>("10861058062577803187713313620936318927813115807431112216295356215798188415993")),FieldT(bigint<FieldT::num_limbs>("6171081416562712860652539550540300841165198894929904648173184407935220453004"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("2973172134096751681569118255635245936822679554194683352822238732713494744526")),FieldT(bigint<FieldT::num_limbs>("10707887432375252298829520257221060883356326036211033062627152466974818597932")),FieldT(bigint<FieldT::num_limbs>("3456931400961642010725970965069072634575668367039273500100334076757881622318"))}));
    ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("171504725307567237300892362739926597082177917037014797958996305423518168062")),FieldT(bigint<FieldT::num_limbs>("17390485639623180175939815091600029917135849399269454264368288075175620197711")),FieldT(bigint<FieldT::num_limbs>("5611514674429950228895331409253879072602192326841866325466758337494475684195"))}));
    return poseidon_params<FieldT>(full_rounds, partial_rounds, alpha, rate, 
        ark_matrix, supported_near_mds, mds_matrix);
}

template<typename FieldT>
poseidon_params<FieldT> high_alpha_128_bit_altbn_poseidon_params(const size_t state_size)
{
    if (state_size > 4 || state_size < 3)
    {
        throw std::invalid_argument("high_alpha_128_bit_altbn_poseidon_params only supports state size 3 or 4");
    }
    const size_t alpha = 17;
    const size_t capacity = 1;
    const size_t rate = state_size - capacity;
    /** Calculated using the same equations for interpolation, and grobner basis eq's for alpha=5
     *  Grobner basis attacks should only get harder as alpha increases, thus this is sound. */
    const size_t full_rounds = 8;
    const size_t partial_rounds = 29;
    const bool supported_near_mds = true;
    std::vector<std::vector<FieldT>> mds_matrix;
    std::vector<std::vector<FieldT>> ark_matrix;
    if (state_size == 3)
    {
        const size_t partial_rounds = 29;
        /** Magic constants obtained from script.
         *  We use a near MDS matrix of the form:
         * [[1, 0, 1],
         *  [1, 1, 0],
         *  [0, 1, 1]]
         */
        mds_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("1")),FieldT(bigint<FieldT::num_limbs>("0")),FieldT(bigint<FieldT::num_limbs>("1"))}));
        mds_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("1")),FieldT(bigint<FieldT::num_limbs>("1")),FieldT(bigint<FieldT::num_limbs>("0"))}));
        mds_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("0")),FieldT(bigint<FieldT::num_limbs>("1")),FieldT(bigint<FieldT::num_limbs>("1"))}));

        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("9478896780421655835758496955063136571251874317427585180076394551808670301829")),FieldT(bigint<FieldT::num_limbs>("1410220424381727336803825453763847584610565307685015130563813219659976870089")),FieldT(bigint<FieldT::num_limbs>("12324248147325396388933912754817224521085038231095815415485781874375379288849"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("5869197693688547188262203345939784760013629955870738354032535473827837048029")),FieldT(bigint<FieldT::num_limbs>("7027675418691353855077049716619550622043312043660992344940177187528247727783")),FieldT(bigint<FieldT::num_limbs>("12525656923125347519081182951439180216858859245949104467678704676398049957654"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("2393593257638453164081539737606611596909105394156134386135868506931280124380")),FieldT(bigint<FieldT::num_limbs>("21284282509779560826339329447865344953337633312148348516557075030360788076689")),FieldT(bigint<FieldT::num_limbs>("9426009547297688316907727916185688178981799243406990694957955230529774886223"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("5340930120720868177469579986808462005013697381998009281661327587975132166755")),FieldT(bigint<FieldT::num_limbs>("13224952063922250960936823741448973692264041750100990569445192064567307041002")),FieldT(bigint<FieldT::num_limbs>("5263772922985715307758718731861278699232625525745635678504665316187832057553"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("12905140589386545724352113723305099554526316070018892915579084990225436501424")),FieldT(bigint<FieldT::num_limbs>("3682692866591423277196501877256311982914914533289815428970072263880360882202")),FieldT(bigint<FieldT::num_limbs>("19681976272543335942352939522328215645129363120562038296975370569202780487598"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("5636115553781577891149626756201577064573186936824720926267940879716772984728")),FieldT(bigint<FieldT::num_limbs>("9501050736957980494328252533770324735114766672253853282051029963140075785396")),FieldT(bigint<FieldT::num_limbs>("2809392708032113981798687947163092027958611686144429680366467696224014505992"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("18433043696013996573551852847056868761017170818820490351056924728720017242180")),FieldT(bigint<FieldT::num_limbs>("1600424531609887868281118752288673305222025191763201214001133841689879221076")),FieldT(bigint<FieldT::num_limbs>("4077863666335789263839414578443702921867998881654209770993100224779179660280"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("10750183821931976144366649760909312224094512474274826686974526305203678408743")),FieldT(bigint<FieldT::num_limbs>("5876585841304782856135279046524906005004905983316552629403091395701737015709")),FieldT(bigint<FieldT::num_limbs>("13484299981373196201166722380389594773562113262309564134825386266765751213853"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("17382139184029132729706972098151128411278461930818849113274328379445169530719")),FieldT(bigint<FieldT::num_limbs>("20539300163698134245746932972121993866388520784246731402041866252259697791654")),FieldT(bigint<FieldT::num_limbs>("149101987103211771991327927827692640556911620408176100290586418839323044234"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3772300053282831651551351000101118094165364582047053942163129913249479587871")),FieldT(bigint<FieldT::num_limbs>("1859494671365748569037492975272316924127197175139843386363551067183747450207")),FieldT(bigint<FieldT::num_limbs>("6056775412522970299341516426839343188000696076848171109448990325789072743616"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("13535861576199801040709157556664030757939966797019046516538528720719863222691")),FieldT(bigint<FieldT::num_limbs>("3166287940256215995277151981354337176516077219230228956292184356796876826882")),FieldT(bigint<FieldT::num_limbs>("3878105211417696553129343540655091450996375987051865710523878345663272335218"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3234972450165117119793849127765475568944145932922109597427102281521349833458")),FieldT(bigint<FieldT::num_limbs>("4245107901241859301876588161430872878162557070919886440605528540426123750702")),FieldT(bigint<FieldT::num_limbs>("14797507122636944484020484450153618519329103538375805997650508264647579279513"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3893725073760673244819994221888005992135922325903832357013427303988853516024")),FieldT(bigint<FieldT::num_limbs>("21641836396029226240087625131527365621781742784615208902930655613239471409203")),FieldT(bigint<FieldT::num_limbs>("4622082908476410083286670201138165773322781640914243047922441301693321472984"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("14738633807199650048753490173004870343158648561341211428780666160270584694255")),FieldT(bigint<FieldT::num_limbs>("2635090520059500019661864086615522409798872905401305311748231832709078452746")),FieldT(bigint<FieldT::num_limbs>("19070766579582338321241892986615538320421651429118757507174186491084617237586"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("12622420533971517050761060317049369208980632120901481436392835424625664738526")),FieldT(bigint<FieldT::num_limbs>("4395637216713203985567958440367812800809784906642242330796693491855644277207")),FieldT(bigint<FieldT::num_limbs>("13856237567677889405904897420967317137820909836352033096836527506967315017500"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("2152570847472117965131784005129148028733701170858744625211808968788882229984")),FieldT(bigint<FieldT::num_limbs>("6585203416839617436007268534508514569040432229287367393560615429950244309612")),FieldT(bigint<FieldT::num_limbs>("2153122337593625580331500314713439203221416612327349850130027435376816262006"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("7340485916200743279276570085958556798507770452421357119145466906520506506342")),FieldT(bigint<FieldT::num_limbs>("12717879727828017519339312786933302720905962296193775803009326830415523871745")),FieldT(bigint<FieldT::num_limbs>("5392903649799167854181087360481925061021040403603926349022734894553054536405"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("7221669722700687417346373353960536661883467014204005276831020252277657076044")),FieldT(bigint<FieldT::num_limbs>("8259126917996748375739426565773281408349947402369855975457055235880500335093")),FieldT(bigint<FieldT::num_limbs>("9272385735015968356236075957906198733226196415690072035874639311675477515202"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("10999027991078055598627757097261950281899485771669414759870674222957875237568")),FieldT(bigint<FieldT::num_limbs>("15453393396765207016379045014101989306173462885430532298601655955681532648226")),FieldT(bigint<FieldT::num_limbs>("5478929644476681096437469958231489102974161353940993351588559414552523375472"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("6864274099016679903139678736335228538241825704814597078997020342617052506183")),FieldT(bigint<FieldT::num_limbs>("12133526413093116990739357861671284889661106676453313677855438696597541491864")),FieldT(bigint<FieldT::num_limbs>("4363234898901124667709814170397096827222883770682185860994495523839008586252"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("16799465577487943696587954846666404704275729737273450161871875150400464433797")),FieldT(bigint<FieldT::num_limbs>("3466902930973160737502426090330438125630820207992414876720169645462530526357")),FieldT(bigint<FieldT::num_limbs>("10062441698891350053170325824989022858836994651376301483266809451301259521913"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("5849282602749563270643968237860161465694876981255295041960826011116890638924")),FieldT(bigint<FieldT::num_limbs>("18460093993858702487671589299005229942046272739124591066186726570539410116617")),FieldT(bigint<FieldT::num_limbs>("9812100862165422922235757591915383485338044715409891361026651619010947646011"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3387849124775103843519196664933515074848119722071551419682472701704619249120")),FieldT(bigint<FieldT::num_limbs>("5283840871671971215904992681385681067319154145921438770232973796570506340281")),FieldT(bigint<FieldT::num_limbs>("14450974197863079729258614455552607708855872944526185987072755641686663205867"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("12613293459867195704822743599193025685229122593088639435739984309110321350551")),FieldT(bigint<FieldT::num_limbs>("6228273556621778927381918766322387348845347649737780310185999880647567569148")),FieldT(bigint<FieldT::num_limbs>("7482296435079443913598332362891173417094991594500715575107878549173583070413"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("18655449861670697203232484600163743308157596453845950955559776266093852537258")),FieldT(bigint<FieldT::num_limbs>("19948920146235041970991269588233091409704340607794045065548049409652881283328")),FieldT(bigint<FieldT::num_limbs>("13866078374565054775555309394949653928903776100036987352339975076159400168494"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("19398653685274645718325650121748668221118186023117741800737442235635318532994")),FieldT(bigint<FieldT::num_limbs>("4234154881267169381851681265196336178292466185695662916289548353755778788440")),FieldT(bigint<FieldT::num_limbs>("12763628380946395634691260884409562631856128057257959813602172954351304541746"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("7882453112990894293341171586279209575183467873317150236705310601775347127762")),FieldT(bigint<FieldT::num_limbs>("5669812778237054435250482766817044415794242063465169363632154286378940417646")),FieldT(bigint<FieldT::num_limbs>("16998738906020038479274018881471127087312245548341958049900081105113388112420"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3923902724726826782251513956816550869721438812970437824859252798290604500141")),FieldT(bigint<FieldT::num_limbs>("8649850619802776810849631749100283821801281306919958924112424995025830909252")),FieldT(bigint<FieldT::num_limbs>("11095642206650177249637693917287763476332497377393343056089442602164577098005"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("6935839211798937659784055008131602708847374430164859822530563797964932598700")),FieldT(bigint<FieldT::num_limbs>("7009671085960032501857416946339379996865118520008277046653124221544059312084")),FieldT(bigint<FieldT::num_limbs>("14361753917538892938870644779277430374939140280641641154553910654644462796654"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("6296738827713642491839335218022320853584196754765009910619998033694434027436")),FieldT(bigint<FieldT::num_limbs>("13849351053619304861036345979638534258290466678610892122310972291285921828452")),FieldT(bigint<FieldT::num_limbs>("434708832289952835651719825370636597763362139118091644948171210201038442144"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("16633750393567936099837698146248798150044883935695159627422586429892098538881")),FieldT(bigint<FieldT::num_limbs>("12944939557587269500508410478785174192748264930676627398550886896505925728421")),FieldT(bigint<FieldT::num_limbs>("13132297714437965464312509267711212830308064898189789451541658159340762509645"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3197382106307730326149017386920960267079843887376371149099833465681078850285")),FieldT(bigint<FieldT::num_limbs>("1219439673853113792340300173186247996249367102884530407862469123523013083971")),FieldT(bigint<FieldT::num_limbs>("3493891993991676033939225547105305872211028239751045376877382816726002847983"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("17474961424148900675164871904345354895260557993970869987490270849177572737815")),FieldT(bigint<FieldT::num_limbs>("14496326112831768456074139601688618143496262542471380389977686658437504436331")),FieldT(bigint<FieldT::num_limbs>("2924472580096769678506212811457662807142794313402961128576445038927398235897"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("4628296006426596599826873705217702584581936573072175641058168144816722698331")),FieldT(bigint<FieldT::num_limbs>("21191637522268746884323101636631937283436518241594045635071026927358145697662")),FieldT(bigint<FieldT::num_limbs>("16951212238971640283544926666565087199118390400059790490897089817025688673127"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("19613695336435411200907478310503966803576648245805018042761984388590288078910")),FieldT(bigint<FieldT::num_limbs>("19408817842355340096520725353160494939342325645253279486424056603334799168015")),FieldT(bigint<FieldT::num_limbs>("21454045045501902703155952158575095010854214688097850310899813261125869452799"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("7770328480231095569114093553841085793308707788942057894109603074902652929530")),FieldT(bigint<FieldT::num_limbs>("16464571997310094273270381226660568195148193554716113613093103468413654931642")),FieldT(bigint<FieldT::num_limbs>("17470702407108506528534764015553049093186219898758900659217736458688524875937"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("18550730212998825286534234924565339469725380540133305684933015562293032312245")),FieldT(bigint<FieldT::num_limbs>("2896017217286658654468296502214232988965841950467453595108246966331694256153")),FieldT(bigint<FieldT::num_limbs>("14675299739240143232464986549869467617250208852063994519435190317578889428919"))}));
        return poseidon_params<FieldT>(full_rounds, partial_rounds, alpha, rate, 
            ark_matrix, supported_near_mds, mds_matrix);
    }
    else if (state_size == 4)
    {
        const size_t partial_rounds = 30;
        /** Magic constants obtained from script.
         *  We use a near MDS matrix of the form:
         * [[0, 1, 1, 1],
         *  [1, 0, 1, 1],
         *  [1, 1, 0, 1],
         *  [1, 1, 1, 0]]
         */
        mds_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("0")),FieldT(bigint<FieldT::num_limbs>("1")),FieldT(bigint<FieldT::num_limbs>("1")),FieldT(bigint<FieldT::num_limbs>("1"))}));
        mds_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("1")),FieldT(bigint<FieldT::num_limbs>("0")),FieldT(bigint<FieldT::num_limbs>("1")),FieldT(bigint<FieldT::num_limbs>("1"))}));
        mds_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("1")),FieldT(bigint<FieldT::num_limbs>("1")),FieldT(bigint<FieldT::num_limbs>("0")),FieldT(bigint<FieldT::num_limbs>("1"))}));
        mds_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("1")),FieldT(bigint<FieldT::num_limbs>("1")),FieldT(bigint<FieldT::num_limbs>("1")),FieldT(bigint<FieldT::num_limbs>("0"))}));

        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("9478896780421655835758496955063136571251874317427585180076394551808670301829")),FieldT(bigint<FieldT::num_limbs>("1410220424381727336803825453763847584610565307685015130563813219659976870089")),FieldT(bigint<FieldT::num_limbs>("12324248147325396388933912754817224521085038231095815415485781874375379288849")),FieldT(bigint<FieldT::num_limbs>("5869197693688547188262203345939784760013629955870738354032535473827837048029"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("7027675418691353855077049716619550622043312043660992344940177187528247727783")),FieldT(bigint<FieldT::num_limbs>("12525656923125347519081182951439180216858859245949104467678704676398049957654")),FieldT(bigint<FieldT::num_limbs>("2393593257638453164081539737606611596909105394156134386135868506931280124380")),FieldT(bigint<FieldT::num_limbs>("21284282509779560826339329447865344953337633312148348516557075030360788076689"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("9426009547297688316907727916185688178981799243406990694957955230529774886223")),FieldT(bigint<FieldT::num_limbs>("5340930120720868177469579986808462005013697381998009281661327587975132166755")),FieldT(bigint<FieldT::num_limbs>("13224952063922250960936823741448973692264041750100990569445192064567307041002")),FieldT(bigint<FieldT::num_limbs>("5263772922985715307758718731861278699232625525745635678504665316187832057553"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("12905140589386545724352113723305099554526316070018892915579084990225436501424")),FieldT(bigint<FieldT::num_limbs>("3682692866591423277196501877256311982914914533289815428970072263880360882202")),FieldT(bigint<FieldT::num_limbs>("19681976272543335942352939522328215645129363120562038296975370569202780487598")),FieldT(bigint<FieldT::num_limbs>("5636115553781577891149626756201577064573186936824720926267940879716772984728"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("9501050736957980494328252533770324735114766672253853282051029963140075785396")),FieldT(bigint<FieldT::num_limbs>("2809392708032113981798687947163092027958611686144429680366467696224014505992")),FieldT(bigint<FieldT::num_limbs>("18433043696013996573551852847056868761017170818820490351056924728720017242180")),FieldT(bigint<FieldT::num_limbs>("1600424531609887868281118752288673305222025191763201214001133841689879221076"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("4077863666335789263839414578443702921867998881654209770993100224779179660280")),FieldT(bigint<FieldT::num_limbs>("10750183821931976144366649760909312224094512474274826686974526305203678408743")),FieldT(bigint<FieldT::num_limbs>("5876585841304782856135279046524906005004905983316552629403091395701737015709")),FieldT(bigint<FieldT::num_limbs>("13484299981373196201166722380389594773562113262309564134825386266765751213853"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("17382139184029132729706972098151128411278461930818849113274328379445169530719")),FieldT(bigint<FieldT::num_limbs>("20539300163698134245746932972121993866388520784246731402041866252259697791654")),FieldT(bigint<FieldT::num_limbs>("149101987103211771991327927827692640556911620408176100290586418839323044234")),FieldT(bigint<FieldT::num_limbs>("3772300053282831651551351000101118094165364582047053942163129913249479587871"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("1859494671365748569037492975272316924127197175139843386363551067183747450207")),FieldT(bigint<FieldT::num_limbs>("6056775412522970299341516426839343188000696076848171109448990325789072743616")),FieldT(bigint<FieldT::num_limbs>("13535861576199801040709157556664030757939966797019046516538528720719863222691")),FieldT(bigint<FieldT::num_limbs>("3166287940256215995277151981354337176516077219230228956292184356796876826882"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3878105211417696553129343540655091450996375987051865710523878345663272335218")),FieldT(bigint<FieldT::num_limbs>("3234972450165117119793849127765475568944145932922109597427102281521349833458")),FieldT(bigint<FieldT::num_limbs>("4245107901241859301876588161430872878162557070919886440605528540426123750702")),FieldT(bigint<FieldT::num_limbs>("14797507122636944484020484450153618519329103538375805997650508264647579279513"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3893725073760673244819994221888005992135922325903832357013427303988853516024")),FieldT(bigint<FieldT::num_limbs>("21641836396029226240087625131527365621781742784615208902930655613239471409203")),FieldT(bigint<FieldT::num_limbs>("4622082908476410083286670201138165773322781640914243047922441301693321472984")),FieldT(bigint<FieldT::num_limbs>("14738633807199650048753490173004870343158648561341211428780666160270584694255"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("2635090520059500019661864086615522409798872905401305311748231832709078452746")),FieldT(bigint<FieldT::num_limbs>("19070766579582338321241892986615538320421651429118757507174186491084617237586")),FieldT(bigint<FieldT::num_limbs>("12622420533971517050761060317049369208980632120901481436392835424625664738526")),FieldT(bigint<FieldT::num_limbs>("4395637216713203985567958440367812800809784906642242330796693491855644277207"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("13856237567677889405904897420967317137820909836352033096836527506967315017500")),FieldT(bigint<FieldT::num_limbs>("2152570847472117965131784005129148028733701170858744625211808968788882229984")),FieldT(bigint<FieldT::num_limbs>("6585203416839617436007268534508514569040432229287367393560615429950244309612")),FieldT(bigint<FieldT::num_limbs>("2153122337593625580331500314713439203221416612327349850130027435376816262006"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("7340485916200743279276570085958556798507770452421357119145466906520506506342")),FieldT(bigint<FieldT::num_limbs>("12717879727828017519339312786933302720905962296193775803009326830415523871745")),FieldT(bigint<FieldT::num_limbs>("5392903649799167854181087360481925061021040403603926349022734894553054536405")),FieldT(bigint<FieldT::num_limbs>("7221669722700687417346373353960536661883467014204005276831020252277657076044"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("8259126917996748375739426565773281408349947402369855975457055235880500335093")),FieldT(bigint<FieldT::num_limbs>("9272385735015968356236075957906198733226196415690072035874639311675477515202")),FieldT(bigint<FieldT::num_limbs>("10999027991078055598627757097261950281899485771669414759870674222957875237568")),FieldT(bigint<FieldT::num_limbs>("15453393396765207016379045014101989306173462885430532298601655955681532648226"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("5478929644476681096437469958231489102974161353940993351588559414552523375472")),FieldT(bigint<FieldT::num_limbs>("6864274099016679903139678736335228538241825704814597078997020342617052506183")),FieldT(bigint<FieldT::num_limbs>("12133526413093116990739357861671284889661106676453313677855438696597541491864")),FieldT(bigint<FieldT::num_limbs>("4363234898901124667709814170397096827222883770682185860994495523839008586252"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("16799465577487943696587954846666404704275729737273450161871875150400464433797")),FieldT(bigint<FieldT::num_limbs>("3466902930973160737502426090330438125630820207992414876720169645462530526357")),FieldT(bigint<FieldT::num_limbs>("10062441698891350053170325824989022858836994651376301483266809451301259521913")),FieldT(bigint<FieldT::num_limbs>("5849282602749563270643968237860161465694876981255295041960826011116890638924"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("18460093993858702487671589299005229942046272739124591066186726570539410116617")),FieldT(bigint<FieldT::num_limbs>("9812100862165422922235757591915383485338044715409891361026651619010947646011")),FieldT(bigint<FieldT::num_limbs>("3387849124775103843519196664933515074848119722071551419682472701704619249120")),FieldT(bigint<FieldT::num_limbs>("5283840871671971215904992681385681067319154145921438770232973796570506340281"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("14450974197863079729258614455552607708855872944526185987072755641686663205867")),FieldT(bigint<FieldT::num_limbs>("12613293459867195704822743599193025685229122593088639435739984309110321350551")),FieldT(bigint<FieldT::num_limbs>("6228273556621778927381918766322387348845347649737780310185999880647567569148")),FieldT(bigint<FieldT::num_limbs>("7482296435079443913598332362891173417094991594500715575107878549173583070413"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("18655449861670697203232484600163743308157596453845950955559776266093852537258")),FieldT(bigint<FieldT::num_limbs>("19948920146235041970991269588233091409704340607794045065548049409652881283328")),FieldT(bigint<FieldT::num_limbs>("13866078374565054775555309394949653928903776100036987352339975076159400168494")),FieldT(bigint<FieldT::num_limbs>("19398653685274645718325650121748668221118186023117741800737442235635318532994"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("4234154881267169381851681265196336178292466185695662916289548353755778788440")),FieldT(bigint<FieldT::num_limbs>("12763628380946395634691260884409562631856128057257959813602172954351304541746")),FieldT(bigint<FieldT::num_limbs>("7882453112990894293341171586279209575183467873317150236705310601775347127762")),FieldT(bigint<FieldT::num_limbs>("5669812778237054435250482766817044415794242063465169363632154286378940417646"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("16998738906020038479274018881471127087312245548341958049900081105113388112420")),FieldT(bigint<FieldT::num_limbs>("3923902724726826782251513956816550869721438812970437824859252798290604500141")),FieldT(bigint<FieldT::num_limbs>("8649850619802776810849631749100283821801281306919958924112424995025830909252")),FieldT(bigint<FieldT::num_limbs>("11095642206650177249637693917287763476332497377393343056089442602164577098005"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("6935839211798937659784055008131602708847374430164859822530563797964932598700")),FieldT(bigint<FieldT::num_limbs>("7009671085960032501857416946339379996865118520008277046653124221544059312084")),FieldT(bigint<FieldT::num_limbs>("14361753917538892938870644779277430374939140280641641154553910654644462796654")),FieldT(bigint<FieldT::num_limbs>("6296738827713642491839335218022320853584196754765009910619998033694434027436"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("13849351053619304861036345979638534258290466678610892122310972291285921828452")),FieldT(bigint<FieldT::num_limbs>("434708832289952835651719825370636597763362139118091644948171210201038442144")),FieldT(bigint<FieldT::num_limbs>("16633750393567936099837698146248798150044883935695159627422586429892098538881")),FieldT(bigint<FieldT::num_limbs>("12944939557587269500508410478785174192748264930676627398550886896505925728421"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("13132297714437965464312509267711212830308064898189789451541658159340762509645")),FieldT(bigint<FieldT::num_limbs>("3197382106307730326149017386920960267079843887376371149099833465681078850285")),FieldT(bigint<FieldT::num_limbs>("1219439673853113792340300173186247996249367102884530407862469123523013083971")),FieldT(bigint<FieldT::num_limbs>("3493891993991676033939225547105305872211028239751045376877382816726002847983"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("17474961424148900675164871904345354895260557993970869987490270849177572737815")),FieldT(bigint<FieldT::num_limbs>("14496326112831768456074139601688618143496262542471380389977686658437504436331")),FieldT(bigint<FieldT::num_limbs>("2924472580096769678506212811457662807142794313402961128576445038927398235897")),FieldT(bigint<FieldT::num_limbs>("4628296006426596599826873705217702584581936573072175641058168144816722698331"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("21191637522268746884323101636631937283436518241594045635071026927358145697662")),FieldT(bigint<FieldT::num_limbs>("16951212238971640283544926666565087199118390400059790490897089817025688673127")),FieldT(bigint<FieldT::num_limbs>("19613695336435411200907478310503966803576648245805018042761984388590288078910")),FieldT(bigint<FieldT::num_limbs>("19408817842355340096520725353160494939342325645253279486424056603334799168015"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("21454045045501902703155952158575095010854214688097850310899813261125869452799")),FieldT(bigint<FieldT::num_limbs>("7770328480231095569114093553841085793308707788942057894109603074902652929530")),FieldT(bigint<FieldT::num_limbs>("16464571997310094273270381226660568195148193554716113613093103468413654931642")),FieldT(bigint<FieldT::num_limbs>("17470702407108506528534764015553049093186219898758900659217736458688524875937"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("18550730212998825286534234924565339469725380540133305684933015562293032312245")),FieldT(bigint<FieldT::num_limbs>("2896017217286658654468296502214232988965841950467453595108246966331694256153")),FieldT(bigint<FieldT::num_limbs>("14675299739240143232464986549869467617250208852063994519435190317578889428919")),FieldT(bigint<FieldT::num_limbs>("18248678941898574458551739591299676471465479218953857934957966658325546357470"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("16024355190461954313137023461204350082170056518459846158694838747631538641098")),FieldT(bigint<FieldT::num_limbs>("7850865488269095796873132539790837027396511819869707882135427785654995088674")),FieldT(bigint<FieldT::num_limbs>("3624709825280601739427579377772539527675482368902561427011111819459495093367")),FieldT(bigint<FieldT::num_limbs>("15151686747202278393509880888871669223363163423030702442491136636685867520291"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("13928718275408319129936199508884886190840388527360771842040395542259777472217")),FieldT(bigint<FieldT::num_limbs>("16438518367400433996754402147971541957362882260288877185296680838641299101277")),FieldT(bigint<FieldT::num_limbs>("14251644764165619982007152928697135211685596697488988945688571805257723263990")),FieldT(bigint<FieldT::num_limbs>("1558415498960552213241704009433360128041672577274390114589014204605400783336"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("10749063905018409753971835881025233122402531635641702583180695671822925094888")),FieldT(bigint<FieldT::num_limbs>("16354988262062662857299847129894034999364829278876526229109898088936579282653")),FieldT(bigint<FieldT::num_limbs>("18674506161104582501429570542497403127395303370792353981110261037140868411363")),FieldT(bigint<FieldT::num_limbs>("257179637705917619258734318345044380508175660880825307674152181190752217180"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("9940762366329880573269658025932923420224913533497554139322991519233208344887")),FieldT(bigint<FieldT::num_limbs>("1293239921425673430660897025143433077974838969258268884994339615096356996604")),FieldT(bigint<FieldT::num_limbs>("3854755057793743998382749043935645076321648392889271506938052408708736403452")),FieldT(bigint<FieldT::num_limbs>("12412189818653171938404326534798862660266918801480759205769134043988251714987"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("14244337044908816699730279807127306408590045371642046203702435001248807053599")),FieldT(bigint<FieldT::num_limbs>("2664474558454245875259863292210158373587436806365583470716142566181853525483")),FieldT(bigint<FieldT::num_limbs>("13804808032959446920445446043982766933818740408606455713901597406432607762297")),FieldT(bigint<FieldT::num_limbs>("15382392877686271824754965441189716482120736348466561135774995773639134591702"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("3837565647574361068761459036361038720977149057378981325662868867989693614839")),FieldT(bigint<FieldT::num_limbs>("15225462598913514623919729309237607260922732019503241269148802966711589239234")),FieldT(bigint<FieldT::num_limbs>("14538071498526966569560493408197490744750860055409168184834291065931235880474")),FieldT(bigint<FieldT::num_limbs>("5410684378713263186548337983252960352059120561867029137014627845977008289113"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("14667331447809547522064866822517742824275266629587157982725082691761010479709")),FieldT(bigint<FieldT::num_limbs>("2791759385703122291826173973326063535105256500818872868576708473427706638753")),FieldT(bigint<FieldT::num_limbs>("3642448711195731276443459026400297797336444189377301111331173182906175572703")),FieldT(bigint<FieldT::num_limbs>("1762188042455633427137702520675816545396284185254002959309669405982213803405"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("2009136114070652405039722031642909652809813532808095287569356370144739527789")),FieldT(bigint<FieldT::num_limbs>("19791483249730874013718425513633316012832876840094804012613683750613457029295")),FieldT(bigint<FieldT::num_limbs>("11215847714185086315143282238988904989100073569422371901331104090526168332344")),FieldT(bigint<FieldT::num_limbs>("6136518117972887266318121761768078738116787014241169482775768363309364395343"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("18745129065505118158068407282927084003185002203481097964691803189568879488516")),FieldT(bigint<FieldT::num_limbs>("4425941001346546123680947386796057642506098630852783901377461933889171987192")),FieldT(bigint<FieldT::num_limbs>("18705881328407019529873437562261994664374735037477656748280619515075691198549")),FieldT(bigint<FieldT::num_limbs>("13047230314802854075744753527622973042529217055807151931674300909891623306887"))}));
        ark_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("6366778494392099924891219059630877678312150007619466052601251344235295985587")),FieldT(bigint<FieldT::num_limbs>("11984731463267744847120482801387059157065563300719282695907722257441888163935")),FieldT(bigint<FieldT::num_limbs>("6190354345246662069799625015516451435683819097684409839182088413597133300258")),FieldT(bigint<FieldT::num_limbs>("3453567349071653221618278943759136615746795332310201741158620901333366005083"))}));
        return poseidon_params<FieldT>(full_rounds, partial_rounds, alpha, rate, 
            ark_matrix, supported_near_mds, mds_matrix);
    }
    else
    {
        throw std::invalid_argument("case already handled earlier");
    }
    
}

}
