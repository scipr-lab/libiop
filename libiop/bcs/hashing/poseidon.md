# Poseidon

This file details our configuration of Poseidon used. We use Poseidon as a choice of hash function for our Merkle Trees, and when it is deployed, it is the bottleneck in time for the entire prover. Consequently, we implement several speed optimizations / protocol changes in our deployment of Poseidon.

Following the notation of the paper, let `t` be the state size, `n` be `log_2(|F|)`, and `alpha` be the exponent used within the S-Box.

## Raising alpha

First we notice that the majority of the partial rounds are needed for security against interpolation attacks, not from security against Grobner attacks. This means that raising alpha potentially lowers the number of partial rounds. We consider the effect raising alpha has on each attack:

* Interpolation attacks: The soundness inequality for a generic alpha (with alpha^{-1} mod (p -1) > alpha) is:
    `R_f + R_p > log_alpha(2) * min(log_2(|F|), lambda) + log_2(t)`
* Grobner attacks: The coefficients of the grobner attack decrease as alpha increases (discussed in the appendices), however for these changes it suffices to assume the constants don't increase. (As the grobner basis attack is nowhere near the limiting factor)
* Statistical attacks: The maximum differential probability of each S-Box increases. The term that changes is the `C` term. In the paper it is given by hand for alpha=3, and alpha=5. In the general case, it is `C = log_2(alpha - 1)`. This follows from the reasoning presented in the paper, that `f(x) = (x + d)^alpha - x^alpha` is a degree `alpha - 1` polynomial, which then has a differential characteristic with probability `(alpha - 1)/n`. Thus we can simply use this definition of C in the prior statistical soundness equations.

We use alpha=17, with the suggested 25% padding for R_f, and a marginally reduced 7.4% padding for R_p. This means that at state size 3, Poseidon (with padding) needs 8 full rounds, and 29 partial rounds. This nearly halves the number of partial rounds required, and concretely improves computation time for our fields.

## Near-MDS matrices

Another optimization we use is Near-MDS matrices. A Near-MDS matrix is a matrix with branch number n, whereas MDS matrices have branch number n + 1. Near-MDS matrices allow for one zero entry per row and column. This notably lowers the in-circuit hash cost under the Fractal metric. 

In the statistical analysis, this changes the differential characteristic over two rounds from `[2^(-n+C)]^(t + 1)` to `[2^(-n+C)]^t`. This changes the statistical soundness equations to be

If `C * t <= nt - lambda` then `R_F = 6 `, else `R_f = 10`.

We thank Lorenzo Grassi for very helpful explanation of how Near-MDS matrices can be used, and how they affect Poseidon's security.

There are multiple works on finding efficient near-MDS matrices, for t > 4: https://tosc.iacr.org/index.php/ToSC/article/view/588, and for t <= 4: https://eprint.iacr.org/2008/395.pdf.

As we currently only use state size <= 4, we use the latter work. It details how at this state size there exist  optimal near-MDS matrix which consist of a maximal number of 0's, with the only non-zero entries being 1. We use circulant matrices in both cases. This means our MDS matrix multiplications no longer require field multiplications, and just require field additions. This results in a 2.2x speed up over normal MDS matrices and alpha=17 at state size 3, and a 4x speedup at alpha=5.

Alternatively, one could use a lightweight MDS matrix, with full branch number. These provide most of the speed benefit (removing multiplications in favor of small exponents that can be computed via additions)