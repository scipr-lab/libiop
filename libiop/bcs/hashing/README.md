# Hashing

This folder contains the hashes currently supported by libiop. Currently it contains Blake2B, Poseidon, and Rescue. 

Three types of hashes are required for usage in the BCS transformation.
* Leaf hashes - These are used to compress MT leaves. Ideally they should support optimizations for arbitrary input sizes.
* Two to one Hashes - These are used to compress MT nodes.
* Hashchain 


## Adding new hashes

To add a new hash, you have to update several locations within hash enum. Namely, the enum itself, and the helper methods to obtain the various hashes of the particular enum.

There is generic infrastructure for converting a secure one way permutation over F^t into a secure hash, via an algebraic sponge construction.