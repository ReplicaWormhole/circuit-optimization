# Coordinator review of the router-deletion obstruction

Board131, by-hand structural review only. No matrix or symbolic search.
After the Bell decoder and the retained router, the unchanged suffix commutes
with Z on each control-register wire u=q1 and v=q3. Its paired X1 wrappers
implement a negative control and do not alter this fact. Its other gates either
act on x=q0/y=q2 or are controlled on u/v.

If the conjugated shift has a nonzero block between two u (or v) sectors,
conjugation by this block-diagonal unitary tail preserves the rank of that
block. It therefore cannot produce a computational diagonal. A simultaneous
output row/column permutation preserves whether off-diagonal entries are zero.
The collaborators' explicit label-exchange witnesses establish such blocks
when either individual XOR router is omitted from the literal construction.

This argument is restricted to omitting a router with the same controlled
suffix. It does not exclude a refit with arbitrary local rotations on u/v,
which need not commute with Z_u/Z_v; nor does it exclude an exact resynthesis
of the router and suffix together, a changed Bell prefix, degenerate-eigenspace
gauge freedom, the original full98 word, or a global eleven-CNOT diagonalizer.
No lower bound or incumbent change follows from this structural obstruction.

## Joint router conjugation and a restricted synthesis obstruction

The collaborators independently confirmed the genuinely modified target
R M R, R=CX(y,v). Its x=1 block is SWAP_yv CZ_yv Sdg_y; its x=0 block is I.
The final selector also changes to F'=R F R, so saving one interaction in M
alone would not establish an eleven-CNOT complete circuit.

For u=0 the transformed final-selector target T' has Y on the even-parity
subspace spanned by00/11 and H on the odd subspace01/10. Thus

T'=(XY+YX)/2+(XX+YY+ZI-IZ)/(2 sqrt(2)).

Its operator coefficient matrix in the Pauli product basis has an X/Y block
with determinant -1/8, plus independent ZI and IZ entries. Hence its operator
Schmidt rank across y|v is four. By contrast, conjugating a single-qubit Pauli
by one CNOT and local gates has rank at most two: a CNOT maps its three axes
into terms with at most two distinct spectator factors; local conjugation
preserves Schmidt rank. The simple three-CNOT template consisting of a one-CX
pair diagonalizer, one controlled local Pauli, and the inverse diagonalizer
cannot implement F'. The independent verifier confirmed signs and rank.

This rules out only that template. It excludes neither general three-CNOT
controlled-matchgate synthesis nor a joint resynthesis of the router, merged
block and final selector. No full circuit was simulated in this review.
