# Exact mixed-cut obstruction near the eight-CNOT tail

Fix the six-CNOT prefix of `topology14_exact_matchgate_rational.json` and
the exact output eigenspace gauge in
`search13_fulltail_invariant_exact_witness.py`. Let $T$ be the resulting
eight-CNOT tail operator. The gauge preserves exact diagonalization of the
four-qubit shift. This note concerns **synthesis of this fixed $T$**,
not all possible seven-CNOT tails or output gauges.

The original chronological tail pair schedule is

`02,02,13,13,01,01,23,23`.

Form a seven-CNOT pair schedule by deleting one of its eight entries and
changing at most one remaining pair to any other physical pair. There
are **124 distinct ordered unoriented schedules**. Each pair allows both
CNOT directions and arbitrary intervening one-qubit gates.

## Target rank certificates

The exact target has operator-Schmidt ranks `(4,4,4,4,14,14,8)` across
single-wire cuts and the balanced cuts `01|23`, `02|13`, `03|12`.
For mixed input/output cuts, use boundary-axis order
`out0,out1,out2,out3,in0,in1,in2,in3`. Its realignment rank is **at
least 14** across each of these two cuts:

| Bit mask | Left boundary axes | Target rank lower bound |
|---|---|---:|
| 83 | `out0,out1,in0,in2` | 14 |
| 163 | `out0,out1,in1,in3` | 14 |

`search13_rank8_neighborhood_mixedcut.py` computes these lower bounds
by evaluating the exact ℚ(ζ₉₆) entries modulo 97 and 193, using
primitive 96th roots 5 and 25 respectively. Every denominator occurring
in the target matrix is invertible at these primes, and the roots satisfy
Φ₉₆. Thus any nonzero minor after reduction certifies a nonzero
characteristic-zero minor. The script records the modular rank for every
one of 127 complementary boundary-cut pairs. An independent review
reconstructed rank 14 at both displayed cuts for both primes.

## Schedule exclusion

A CNOT circuit is a tensor network with dimension-two wire bonds.
Assigning its seven CNOT vertices to either side of a boundary cut gives
a rank upper bound of $2^c$, where $c$ is the number of crossed wire
bonds. Minimizing over assignments gives the ordinary graph mincut.
For the 124 schedules, exact ordinary output/input cut capacities reject
94. The gauge-independent commutant condition rejects none of the
remaining 30. Mixed-cut wire-bond mincuts reject 22 more: each has a
cut of size two at mask 83 or 163, giving rank at most four against the
target's rank at least 14.

For the eight remaining schedules, split every CNOT tensor exactly as

\[
\operatorname{CX}_{c\to t}
 =\sum_{b=0}^{1} |b\rangle\!\langle b|_c\otimes X_t^b.
\]

This creates a further dimension-two bond between the control and
target factors. One-qubit gates may be absorbed into wire segments.
`search13_rank8_directional_mincut.py` checks all $2^7=128$ direction
assignments for each of the eight schedules. Every assignment has a
mincut of **three at both masks 83 and 163**, hence realignment rank at
most eight. Since the target rank is at least 14, all eight schedules
are excluded. The direction-aware graph is an exact tensor-network
factorization, so the bound holds for arbitrary one-qubit gates.

Together, the two exact screens exclude **all 124 schedules in this
one-deletion/one-rewire neighborhood** from synthesizing the fixed
seven-CNOT tail. This does **not** exclude other seven-CNOT schedules,
other output eigenspace gauges, a changed six-CNOT prefix, or a general
13-CNOT cycle diagonalizer. The global interval remains
$6\le C_{\min}\le14$.

Reproduce with:

```bash
python3 search13_fulltail_invariant_exact_witness.py
python3 search13_rank8_neighborhood_mixedcut.py
python3 search13_rank8_directional_mincut.py
python3 experiment_log.py show 284
python3 experiment_log.py show 286
```
