use algebra::{
    bytes::ToBytes, 
    };
use std::{fmt, rc::Rc};

#[cfg(feature = "r1cs")]
pub mod constraints;

pub struct MerkleTreePath<F: Field> {
    pub(crate) path: Vec<(
        <P::H as FixedLengthCRH>::Output,
        <P::H as FixedLengthCRH>::Output,
    )>,
}

/// Stores the hashes of a particular path (in order) from leaf to root.
/// Our path `is_left_child()` if the boolean in `path` is true.
#[derive(Derivative)]
#[derivative(
    Clone(bound = "P: MerkleTreeConfig"),
    Debug(bound = "P: MerkleTreeConfig, <P::H as FixedLengthCRH>::Output: fmt::Debug")
)]
pub struct MerkleTreePath<P: MerkleTreeConfig> {
    pub(crate) path: Vec<(
        <P::H as FixedLengthCRH>::Output,
        <P::H as FixedLengthCRH>::Output,
    )>,
}
