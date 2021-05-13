#include <cstdint>

#include <gtest/gtest.h>

#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>

#include "libiop/bcs/hashing/hashing.hpp"
#include "libiop/bcs/hashing/poseidon.hpp"
#include "libiop/bcs/hashing/algebraic_sponge.hpp"
#include "libiop/bcs/hashing/hash_enum.hpp"

namespace libiop {

template<typename FieldT>
poseidon_params<FieldT> default_params()
{
    const size_t full_rounds = 6;
    const size_t partial_rounds = 6;
    const size_t alpha = 5;
    const size_t state_size = 3;
    const size_t rate = 2;
    const size_t capacity = state_size - rate;
    const bool supported_near_mds = false;
    std::vector<std::vector<FieldT>> mds_matrix;
    mds_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("8825587153275961389267911596064203545440508865794516775909719303807406035696")),FieldT(bigint<FieldT::num_limbs>("10687073437533017977898012181774660682972097507111480941190856895646559858118")),FieldT(bigint<FieldT::num_limbs>("5875859313364614113662523387083632673129873319244972245111852429828028656578"))}));
    mds_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("13019783996669809080274352004727116845804189744832555486912893723817400034038")),FieldT(bigint<FieldT::num_limbs>("8721110299256301289286672573897056895986105760217886314397823390102897986931")),FieldT(bigint<FieldT::num_limbs>("16503839855514153794855421052971932410310146883135771153026909577369893787129"))}));
    mds_matrix.push_back(std::vector<FieldT>({FieldT(bigint<FieldT::num_limbs>("8677639611242538259717954365467396258187948401773964682625641787971040155167")),FieldT(bigint<FieldT::num_limbs>("2088218044364962781629079623906443285931855688684056637520251653768172553894")),FieldT(bigint<FieldT::num_limbs>("18238362791116776129630677184049089192261109185329641710406304984813506470419"))}));
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
    poseidon_params<FieldT> params(full_rounds, partial_rounds, alpha, rate, ark_matrix, supported_near_mds, mds_matrix);
    return params;
}

TEST(PermutationTest, PoseidonTest) {
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;

    poseidon_params<FieldT> params = default_params<FieldT>();
    poseidon<FieldT> poseidon_sponge(params);

    /* Empty state */
    FieldT hash_input = FieldT::zero();
    FieldT result = poseidon_sponge.squeeze_vector(params.capacity_)[0];
    FieldT expected = FieldT(bigint<FieldT::num_limbs>("11513774210489128719203754000419293109474869123660673521809718785157314013443"));
    ASSERT_TRUE(result == expected);

    poseidon_sponge.reset();
    result = poseidon_sponge.squeeze_vector(params.capacity_)[0];
    ASSERT_TRUE(result == expected);

    poseidon_params<FieldT> params2 = high_alpha_128_bit_altbn_poseidon_params<FieldT>();
    poseidon<FieldT> poseidon_sponge2(params2);
    result = poseidon_sponge2.squeeze_vector(params.capacity_)[0];
    expected = FieldT(bigint<FieldT::num_limbs>("19745903574422741006139475519330790957027605504238596103618584028691101830733"));
    ASSERT_TRUE(result == expected);
}

TEST(LeafTest, PoseidonTest) {
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;

    poseidon_params<FieldT> params = default_params<FieldT>();
    poseidon<FieldT> poseidon_sponge(params);
    algebraic_leafhash<FieldT> leafhasher(std::make_shared<poseidon<FieldT>>(poseidon_sponge), 128);

    /* Empty hash - same as permutation on empty state */
    std::vector<FieldT> hash_input(1, FieldT::zero());
    FieldT result = leafhasher.hash(hash_input);
    FieldT expected = FieldT(bigint<FieldT::num_limbs>("11513774210489128719203754000419293109474869123660673521809718785157314013443"));
    ASSERT_TRUE(result == expected);

    /* Check that it resets state */
    result = leafhasher.hash(hash_input);
    ASSERT_TRUE(result == expected);

    /* Check that zk salt being 0 works */
    using namespace std::string_literals;
    zk_salt_type salt("\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00"s);
    result = leafhasher.zk_hash(hash_input, salt);
    ASSERT_TRUE(result == expected);

    /** Check that zk salt being 'A'*8 + 'B'*8 + 'C'*8 + 'D'*8 works as expected
     *  This is ensures we parse the salt correctly. */
    salt = std::string("AAAAAAAABBBBBBBBCCCCCCCCDDDDDDDD"s);
    // Actual number
    FieldT salt_as_field_elem = FieldT(bigint<FieldT::num_limbs>("29515630589904128245248592656078826240104804215908883401742459362001266426948"));
    expected = leafhasher.hash(std::vector<FieldT>({hash_input[0], salt_as_field_elem}));    
    result = leafhasher.zk_hash(hash_input, salt);
    ASSERT_TRUE(result == expected);
}

TEST(TwoToOneTest, PoseidonTest) {
    libff::alt_bn128_pp::init_public_params();
    typedef libff::alt_bn128_Fr FieldT;

    std::shared_ptr<leafhash<FieldT, FieldT>> leafhasher = get_leafhash<FieldT, FieldT>(starkware_poseidon_type, 128, 2);
    two_to_one_hash_function<FieldT> twoToOneHash = get_two_to_one_hash<FieldT, FieldT>(starkware_poseidon_type, 128);

    std::vector<FieldT> hash_input(2, FieldT::zero());
    FieldT expected = leafhasher->hash(hash_input);
    FieldT result = twoToOneHash(hash_input[0], hash_input[1], 32);
    ASSERT_TRUE(result == expected);
    /** Check determinism */
    result = twoToOneHash(hash_input[0], hash_input[1], 32);
    ASSERT_TRUE(result == expected);
    result = leafhasher->hash(hash_input);
    ASSERT_TRUE(result == expected);
}

}