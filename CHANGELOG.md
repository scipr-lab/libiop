# Changelog

## Pending

## v0.2.1

* Upgrade Cmake to version 3.1
* Bug fixes in Fractal/Aurora zero knowledge degree bound settings.

## v0.2.0

This update implements Poseidon, and adds significantly algebraic hash infrastructure.
The BCS transformation now has an option for proofs of work for shrinking proof size.
The Fractal RS-IOP has been updated, and is now significantly faster than before.

Many other changes throughout the library.

## v0.1.0

This update implements Fractal, a transparent preprocessing SNARK with a polylog verifier.
It also includes significant performance and memory improvements.

Aurora has gotten significantly faster. The prover is 2.2 times faster, and the verifier is 2.4 times faster.
Ligero has gotten 25% faster in both the prover and verifier times.
The Fractal prover has gotten 20% faster than the stated numbers in the paper.

Some highlights of the optimizations from the update:

Algebra Speedups:
* FFTs are now Nlog(d) instead of NlogN
* FFTs now utilize a cache, which halves the number of multiplications

IOP:
* Switch all codewords to being passed around as shared pointers,
  which significantly reduces memory consumption

Merkle Tree:
* It now handles serializing by coset, reducing time to copy memory

LDT Reducer:
* Improve efficiency of shifting sub-maximal degree oracles

FRI:
* The prover does a single batch inversion per round,
  instead of one batch inversion per coset of size 2^{localization param}
* Optimized the FRI lagrange interpolation equations to further reduce the linear factor
* The FRI argument size optimizer now accounts for MT pruning, and tends to yield 3% better argument sizes

Rowcheck:
* Take advantage of Z_H | L being an |H| to 1 map, to speedup computing (1 / Z_H) over L.

Aurora Lincheck:
* Switched the verifier lagrange interpolation to IFFTs.
  Despite the worse asymptotics, this concretely performs much better.

Sumcheck:
* Make the IFFT over a domain of minimal size, instead of the entire codeword
* Improve efficiency of sumcheck virtual oracles,
  by checking an equi-satisfiable rational function that is easier to compute

## v0.0.0

Initial release