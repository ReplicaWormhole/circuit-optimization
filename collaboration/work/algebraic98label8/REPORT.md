# Run395 attenuated-root residual structure

Hypothesis91; read-only interpretation of saved run395 outputs and the two
run395 reports. I performed no new matrix or gradient evaluation, numerical or
symbolic search, optimization, or certificate attempt.

## Exact relation between the two reported losses

Let `B=U V4 U†`, which is unitary. The independent run395 audit reports eight
diagonal entries at the four exact roots (counts +1:4, -1:2, +i:1, -i:1) and
eight entries aligned with those roots but attenuated by

`a=(2+sqrt(2))/4 = cos²(pi/8) ~= 0.8535533905932737`.

The attenuated entries have two rows per root. This is consistent with the
full spectrum multiplicities: assigning the remaining two slots per root to
these eight rows gives the minimizing diagonal labels. The report values are
numerical observations from the saved output, not new symbolic certification
of exact equalities.

For any unitary `B`, row normalization gives

`sum_(s != r) |B_rs|² = 1-|B_rr|²`.

Thus the eight unit-modulus diagonal rows contribute no offdiagonal norm, and
the eight attenuated rows give

`L_off = [8(1-a²)]/16 = (1-a²)/2 ~= 0.135723304703363`.

With the optimal assigned root on each active row, the diagonal mismatch is
`1-a`, so

`L_spec-L_off = [8(1-a)²]/16 = (1-a)²/2 ~= 0.010723304703363`,

and consequently

`L_spec = [(1-a²)+(1-a)²]/2 = 1-a = (2-sqrt(2))/4 ~= 0.146446609406726`.

This explains why the labeled loss exceeds the original loss while both are
nonzero at the run395 endpoint. The optimizer returned after one iteration and
three evaluations, so this structure describes the saved stationary warm
start, not a successful refinement.

The actual raw per-row assignment cost in run395 is
`sum_s |(UV4)_rs-lambda U_rs|²`; summing selected raw row costs already equals
`||UV4-DU||_F²` and includes all offdiagonal contributions. It must not have
`L_off` added again. The reduced Hungarian cost `|B_rr-lambda|²` contains only
the assignment-dependent diagonal term; adding `L_off` recovers the full
normalized minimum. This corrects any terse handoff that could be read as
saying raw row costs omit `L_off`.

## Qualified residual-block hypothesis

If the eight reported unit-modulus diagonal entries are exact, unitarity forces
their rows' offdiagonal entries to vanish: a unit vector row with one component
of modulus one has all others zero. The corresponding columns vanish off the
diagonal as well. After a common row/column permutation, `B` therefore splits
as an eight-dimensional diagonal spectator block and an eight-dimensional
unitary active block. The saved checker values are within floating-point error
of unit modulus, so approximate spectator decoupling is supported numerically;
this report does not assert an exact invariant subspace.

In the active block, the observed diagonal is `a` times its assigned roots,
with multiplicity two for each of `+1,-1,+i,-i`. This suggests a structural
synthesis hypothesis: first realize an eight-dimensional spectator eigenspace
and then synthesize the balanced active eight-dimensional mixing that spreads
the remaining root sectors while maintaining these diagonal overlaps. The
native 4-CX/4-XXYY schedule could be inspected for a gauge or parameter
submanifold preserving the spectator sector, then its remaining freedom used
to implement the active block. This may be more actionable than treating all
98 coordinates as an undifferentiated fit, but no such parameter submanifold
or decomposition has been established.

Unproved parts: the displayed near-unit diagonals may only be approximate;
the current endpoint's full `B` is not certified block diagonal; the active
block's offdiagonal pattern and its decomposition into schedule gates are
unknown; and no evidence shows a nearby exact diagonalizer or that the fixed
schedule realizes a chosen active-block target. The observed `maxoff` near
`1/(2 sqrt(2))` is compatible with distributed coupling but does not identify
its graph or prove a block ansatz. The candidate remains invalid and no lower
bound or topology exclusion follows.

## Actionable next checks

1. In an independently audited endpoint, classify rows whose diagonal
   magnitude is numerically one, check their leakage and corresponding column
   leakage, and permute them into an explicit spectator/active partition.
2. If the partition persists across fresh starts, derive the active-block
   matrix pattern from saved full lists and ask whether the schedule's local
   layers can preserve the spectator sector exactly.
3. Treat an exact reduced synthesis as a new hypothesis requiring its own
   reservation and independent full gate-list and exact-certificate checks.

No run was performed or requested by this read-only hypothesis.
