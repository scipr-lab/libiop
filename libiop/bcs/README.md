# BCS Transformation

This directory contains code for a slightly modified variant of the BCS transformation from [BCS16] and [COS20], which is a post-quantum, hash-based compiler from IOP to SNARK. 

The way this code is structured is that it contains definitions for BCS Prover and BCS Verifier classes. These classes are used by the calling SNARK, who uses it throughout the IOP. All of the 'interactive' steps in the IOP (getting randomness, queries, submitting prover mesasges and oracles) are done as calls to these classes. The interface for these calls is as defined in [/libiop/iop/iop.hpp].

The BCS transformation has essentially two parts, a hashchain to absorb all prover messages and squeeze verifier randomness, and merkle trees to compress the oracles emitted by the prover.

The BCS transformation supports using distinct hashes for MT leaves, MT compression and the hashchain, the interfaces are as defined in [/libiop/bcs/hashing/hashing.hpp]. The library supports using both algebraic and binary hashes, and using different hashes for each type.

## Variations over the standard BCS Transformation

* We implement a Proof of Work for the final round.
This directly boosts the round-by-round soundness of the final round, which in parameterization lets you then lower the number of queries and thereby overall verifier time/proof size.
A more intuitive explanation for how this works is that if the adversary fails at breaking interactive soundness, then its best bet is to grind over queries until they probabilistically all fall into places that make the IOP happen to be true.
However the proof of work significantly increases the cost of each trying a new query set. 
If the adversary does W work on average to solve the proof of work, then in time T, it can only brute force over T/W different query sets. Therefore we can lower the number of queries per set, since the T-time adversary gets fewer tries.

* Coset-hashing: In FRI, the verifier will query at structured related points. E.g. x, g*x. We make these points be in the same leaf, since they will always be queried together. This is a notable savings on BCS proving time.

* If you want a zero knowledge SNARK, the papers prove zero knowledge by adding a salt to every leaf, in all rounds. We instead implement salts only for rounds that have oracles that must be kept zero knowledge, as specified by the IOP. This is sufficient for zero-knowledge.

## TODO For usability

A major feature that is missing from this library at the moment is a method of serializing the final transcript. The transcript for all protocols is as defined in bcs_common.hpp, and just needs a standardized method for encoding it.