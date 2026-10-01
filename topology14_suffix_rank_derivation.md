# Fixed final Fourier suffix needs five CNOTs

Use four wires `0,1,2,3`, with `0` most significant. The final five-CNOT
suffix of `baseline_18.json` begins at gate 38. Its CNOTs are `2→1`,
`0→1` twice, and `2→3` twice.

The initial `H_1 CX(2,1) H_1` is `CZ(1,2)`. The remaining gates form
independent two-qubit blocks `F_01` and `F_23`. Their intervening `Rz`
rotations on wires 1 and 3 commute through gates on the other pair. Hence,
up to an input-side product of one-qubit gates, the suffix is

`S = (F_01 ⊗ F_23) CZ(1,2)`.

One-qubit gates on either side do not change operator-Schmidt rank across any
wire partition. Computing `F_01` and `F_23` exactly from the rational-π gates
in `baseline_18.json` over `Q(ζ_16)`, the 16-by-16 realignment determinants
of `S` are exactly **1** for each cut `02|13` and `03|12`. Thus both cuts have
operator-Schmidt rank 16. Across `01|23`, the two `F` blocks are local and
the sole crossing operator is `CZ(1,2)`, so the rank is at most 2 (and is 2).

Suppose four CNOTs synthesized this fixed suffix. Rank 16 on each of the
first two cuts requires all four CNOTs to cross each cut. The only wire pairs
crossing *both* cuts are `01` and `23`; therefore each of the four CNOTs
would connect either `0,1` or `2,3`. Such a circuit has rank 1 across
`01|23`, contradicting the suffix's rank 2. So the fixed suffix requires at
least **five CNOTs**. Its existing circuit meets that bound.

This is a bound on this fixed suffix operator, not on every full `V4`
diagonalizer. Numerical rank checks for the four different 15-CNOT catalog
candidates (`fivemask15_topology_0.json` through `_3.json`) found the same
obstruction at their final five-CNOT suffixes. For a replacement of the
*final six-CNOT suffix* by five CNOTs, the rank-cut filter admits four
undirected five-edge patterns for catalog circuits 0/1 and eight for 2/3.
Passing the filter is necessary, not sufficient. The first viable cut in
catalog circuits 2/3 is at gate 27, leaving nine CNOTs in the prefix.

Reproduce the exact certificate with
`python3 topology14_tail_cut_exact.py` (ledger run 102). Numeric checks and
topology enumeration are in `topology14_tail_cut_bound.py` (run 97),
`topology14_alt_suffix_rank.py` (run 105), and
`topology14_earlier_split_rank.py` (run 107). The actual suffix factorization
was independently checked as matrices to maximum entry error `2.48e-16`.
