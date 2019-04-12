<h1 align="center">libiop: a C++ library for IOP-based zkSNARKs</h1>
<p align="center">
   <a href="https://github.com/scipr-lab/libiop/blob/master/AUTHORS"><img src="https://img.shields.io/badge/authors-SCIPR%20Lab-orange.svg"></a>
   <a href="https://github.com/scipr-lab/libiop/blob/master/LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
</p>

This library provides zkSNARK constructions that are __transparent__ and __post-quantum__, and moreover rely only on __lightweight symmetric cryptography__ (any [cryptographic hash function](https://en.wikipedia.org/wiki/Cryptographic_hash_function)).

The library provides a tool chain for transforming certain types of probabilistic proofs (see below) into zkSNARKs that have the above properties. The library includes two zkSNARK constructions that follow this paradigm:

* The __Aurora__ protocol from [[BCRSVW]](https://eprint.iacr.org/2018/828), whose argument size is O(log<sup>2</sup> N);
* The __Ligero__ protocol from [[AHIV]](https://acmccs.github.io/papers/p2087-amesA.pdf), whose argument size is O(N<sup>0.5</sup>).

Both of these zkSNARKs support R1CS, an NP-complete relation that generalizes arithmetic circuit satisfiability. An important component of Aurora, which is of independent interest, is the [FRI proximity test](https://eccc.weizmann.ac.il/report/2017/134/). These protocols all support binary extension fields, and smooth prime fields.

<span style="color:red">**WARNING:**</span> This is an academic proof-of-concept prototype, and in particular has not received careful code review. <br> This implementation is NOT ready for production use.


## Non-interactive arguments from IOPs

[Interactive oracle proofs](https://eprint.iacr.org/2016/116) (IOPs) are a multi-round generalization of [probabilistically checkable proofs](https://en.wikipedia.org/wiki/Probabilistically_checkable_proof) (PCPs). Known IOP protocols have much better (asymptotic and concrete) performance compared to PCP protocols.

A method known as the [BCS transformation](https://eprint.iacr.org/2016/116) uses a cryptographic hash function (modeled as a [random oracle](https://en.wikipedia.org/wiki/Random_oracle)) to transform any public-coin IOP into a non-interactive argument. This transformation has the benefit that the resulting non-interactive argument is:

* __transparent__, because the only global parameter needed to produce/validate proof strings is the hash function; and
* __post-quantum__, because no attacks are known on the BCS transformation.

Moreover, since no cryptography beyond the hash function is used, non-interactive arguments obtained via the BCS transformation are relatively __lightweight__ (in terms of computational resources). Finally, the BCS transformation __preserves zero knowledge__, in the sense that if the underlying IOP is (honest verifier) zero knowledge then the resulting non-interactive argument is zero knowledge.

This library provides zkSNARKs obtained via the BCS transformation.

## IOP protocols

The folder [`libiop/iop`](libiop/iop) contains infrastructure for writing IOP protocols.

The folder [`libiop/protocols`](libiop/protocols) contains several protocols written using this infrastructure. These include:

<table align="center">
  <tr>
    <th></th>
    <th>language</th>
    <th>round<br>complexity</th>
    <th>oracle length<br>(field elts)</th>
    <th>query<br>complexity</th>
    <th>prover time<br>(field ops)</th>
    <th>verifier time<br>(field ops)</th>
  </tr>
  <tr>
    <td>Ligero-IOP</td>
    <td>R1CS</td>
    <td>2</td>
    <td>O(N)</td>
    <td>O(N<sup>0.5</sup>)</td>
    <td>O(N logN)</td>
    <td>O(N)</td>
  </tr>
  <tr>
    <td>Aurora-IOP</td>
    <td>R1CS</td>
    <td>O(log N)</td>
    <td>O(N)</td>
    <td>O(log N)</td>
    <td>O(N logN)</td>
    <td>O(N)</td>
  </tr>
</table>

The first is an IOP from the [Ligero paper](https://eprint.iacr.org/2018/828),<sup>[1]</sup> and the second is an IOP from the [Aurora paper](https://eprint.iacr.org/2018/828).

Efficient IOP protocols such as the above are obtained by combining two components: (1) RS-encoded IOP, and a (2) proximity test for the RS code. (See [this paper](https://eprint.iacr.org/2018/828) for more details.) The codebase in this library provides infrastructure that enables generically composing these components.

* The folder [`libiop/protocols/encoded`](libiop/protocols/encoded) contains RS-encoded IOPs. This includes the RS-encoded IOPs that form the core of the Aurora and Ligero protocols.
* The folder [`libiop/protocols/ldt`](libiop/protocols/ldt) contains proximity tests for the RS code. This includes a _direct test_ (used by Ligero) and the _[FRI protocol](https://eccc.weizmann.ac.il/report/2017/134/)_ (used by Aurora and other IOPs).

<sup>[1]</sup>: More precisely, the Ligero paper only describes a construction for _arithmetic circuits_. An appendix of the Aurora paper explains how to extend the construction to support R1CS. The latter is the implemented protocol.

## BCS transformation

The folder [`libiop/snark/common`](libiop/snark/common) contains the BCS transformation as a standalone component.

The folder [`libiop/snark`](libiop/snark) contains zkSNARKs obtained by applying the BCS transformation to the IOP protocols above.

<table align="center">
  <tr>
    <th></th>
    <th>language</th>
    <th>prover time</th>
    <th>argument size</th>
    <th>verifier time</th>
  </tr>
  <tr>
    <td>Ligero-SNARK</td>
    <td>R1CS</td>
    <td>O<sub>&kappa;</sub>(N logN)</td>
    <td>O<sub>&kappa;</sub>(N<sup>0.5</sup>)</td>
    <td>O<sub>&kappa;</sub>(N)</td>
  </tr>
  <tr>
    <td>Aurora-SNARK</td>
    <td>R1CS</td>
    <td>O<sub>&kappa;</sub>(N logN)</td>
    <td>O<sub>&kappa;</sub>(log<sup>2</sup> N)</td>
    <td>O<sub>&kappa;</sub>(N)</td>
  </tr>
</table>
&kappa; is used to denote the asymptotics dependence on the security parameter explicitly.

A flag `make_zk` can be set to indicate that the transformation should preserve zero knowledge, or not set to indicate that the IOP being transformed is not zero knowledge and so there is no need to preserve zero knowledge.

## Installation

Please follow the [installation instructions](INSTALL.md).

## Testing

Test files are in [`libiop/tests`](libiop/tests).

For example, to run all of the tests for the Aurora protocol, do the following:

```bash
	$ ./test_aurora_snark
	$ ./test_aurora_protocol # Tests the RS-encoded IOP
	$ ./test_aurora_protocol_components # More granular tests
	$ ./test_aurora_multi_lincheck
	$ ./test_aurora_rowcheck
	$ ./test_aurora_sumcheck
```

## Profiling

The folder [`libiop/profiling`](libiop/profiling) contains tooling to produce protocol execution traces with timing and argument size information.
For example, we can create traces for Aurora and Ligero over a 181 bit prime field, with RS extra dimensions=3 with the following commands:

```bash
  $ ./instrument_aurora_snark --make_zk 1 --is_multiplicative 1 --field_size=181 --optimize_localization=1
  $ ./instrument_ligero_snark --make_zk 1 --is_multiplicative 1 --field_size=181 --RS_extra_dimensions=3
```

We use traces generated from the above commands to create the following plots:
<p align="center"><img src="https://user-images.githubusercontent.com/6440154/55777698-83b66700-5a55-11e9-9146-f2237bc681bd.jpg" alt="argument size" width="40%"/></p>
<p align="center"><img src="https://user-images.githubusercontent.com/6440154/55777717-9a5cbe00-5a55-11e9-8d7a-7b3245c45012.jpg" alt="prover time" width="40%"/><img src="https://user-images.githubusercontent.com/6440154/55777724-a34d8f80-5a55-11e9-81fa-c695b875accf.jpg" alt="verifier time" width="40%"/></p>

## License

This library is licensed under the [MIT License](LICENSE).

## Acknowledgements

This work was supported by:
a Google Faculty Award;
the Israel Science Foundation;
the UC Berkeley Center for Long-Term Cybersecurity;
and donations from the Ethereum Foundation, the Interchain Foundation, and Qtum.
