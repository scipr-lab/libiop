use algebra::{
    biginteger::BigInteger256 as BigInteger,
    fields::{Field, PrimeField, Fp256, Fp256Parameters, FpParameters},
};

use num_traits::identities::One;

pub type Fr = Fp256<FrParameters>;

pub struct FrParameters;

impl Fp256Parameters for FrParameters {}
impl FpParameters for FrParameters {
    type BigInt = BigInteger;

    // MODULUS = 21888242871839275222246405745257275088548364400416034343698204186575808495617
    const MODULUS: BigInteger = BigInteger([
        4891460686036598785u64,
        2896914383306846353u64,
        13281191951274694749u64,
        3486998266802970665u64,
    ]);

    const MODULUS_BITS: u32 = 254;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = 2;

    // 6350874878119819312338956282401532410528162663560392320966563075034087161851
    const R: BigInteger = BigInteger([
        12436184717236109307u64,
        3962172157175319849u64,
        7381016538464732718u64,
        1011752739694698287u64,
    ]);

    const R2: BigInteger = BigInteger([
        1997599621687373223u64,
        6052339484930628067u64,
        10108755138030829701u64,
        150537098327114917u64,
    ]);

    const INV: u64 = 14042775128853446655u64;

    // 5
    const GENERATOR: BigInteger = BigInteger([
        0u64,
        0u64,
        0u64,
        5u64,
    ]);

    const TWO_ADICITY: u32 = 28;

    // 9103219067921713944291392827692070036145651957329286315305642004821462161904
    const ROOT_OF_UNITY: BigInteger = BigInteger([
        11229192882073836016u64,
        14463146568527388372u64,
        17977229792282476866u64,
        1450226466237278415u64,
    ]);

    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        11669102379873075200u64,
        10671829228508198984u64,
        15863968012492123182u64,
        1743499133401485332u64,
    ]);

    const T: BigInteger = BigInteger([
        11211439779908376895u64,
        1735440370612733063u64,
        1376415503089949544u64,
        12990080814u64,
    ]);

    const T_MINUS_ONE_DIV_TWO: BigInteger = BigInteger([
        14829091926808964255u64,
        867720185306366531u64,
        688207751544974772u64,
        6495040407u64,
    ]);
}

#[test]
fn test_alt_bn128_fr() {
    // use algebra::fields::tests::{field_test, primefield_test};
    use crate::alt_bn128::fr::Fr;
    use std::str::FromStr;

    let a: Fr = Fr::one();
    let b: Fr = Fr::one();
    let c: Fr = a + b;
    let exp: Fr = Fr::from_str("2").map_err(|_| ()).unwrap();
    assert_eq!(c, exp);
}
