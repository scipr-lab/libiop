// use algebra::fields::Field;
use crate::alt_bn128::Fr;
use digest::Digest;

use super::PRF;
use crate::CryptoError;

#[cfg(feature = "r1cs")]
pub mod constraints;

#[derive(Clone)]
pub struct Poseidon;

impl PRF for Poseidon {
    type Input = [Fr; 3];
    type Output = [Fr; 3];
    

    fn evaluate(input: &Self::Input) -> Result<Self::Output, CryptoError> {
        let eval_time = start_timer!(|| "Poseidon::Eval");
        let mut h = b2s::new();
        h.input(seed.as_ref());
        h.input(input.as_ref());
        let mut result = [0u8; 32];
        result.copy_from_slice(&h.result());
        end_timer!(eval_time);
        Ok(result)
    }
}
