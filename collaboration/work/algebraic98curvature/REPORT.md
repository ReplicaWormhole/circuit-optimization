# Curvature of the minimum spectrum-assignment objective

Hypothesis97. This is a read-only derivation for interpreting a fixed-label
Hessian diagnostic. No matrices, assignments, Hessians, derivatives, optimizer
steps, or certificates were evaluated here.

## When a fixed-label Hessian is the Hessian of the minimum

For circuit coordinates `x`, let `f_j(x)=||U(x)V4-D_jU(x)||_F²/16` for every
admissible root-label vector `j`, with exact root capacities 6,4,3,3. The
spectrum-assignment objective is the finite lower envelope

`f(x)=min_j f_j(x)`.

If one root-label vector `j0` is the unique minimizer at `x0`, define its
assignment gap

`Delta(x0)=min_(j != j0) [ f_j(x0)-f_j0(x0) ] > 0`.

There are finitely many distinct label vectors, each branch is continuous
(and smooth for a smooth circuit chart), so strict positive gap implies that
`j0` remains minimizing throughout some neighborhood of `x0`. In that
neighborhood `f=f_j0`, and its ordinary Hessian is exactly the Hessian of the
minimum-assignment objective. A small positive gap means the neighborhood may
be correspondingly small; assignment stability should not be extrapolated
past it.

The Hungarian implementation expands each root into repeated identical
columns. Permuting columns that carry the same root leaves the label vector
`(d_r)` and the matrix `D` unchanged. Such duplicate-slot ties are not distinct
branches of `f`; they do not invalidate the fixed-D Hessian. A genuine branch
tie requires a different row-to-root vector, hence a different `D`.

At a tie between distinct root-label vectors, `f` is the minimum of smooth
branches and need not be differentiable. Its directional first derivative is
the minimum of the directional derivatives of the active branches. If active
gradients differ, there is a kink, and a classical Hessian of `f` is not
defined there. A fixed-D Hessian then describes only that selected smooth
branch; it is not automatically the curvature of the lower envelope. If the
active gradients happen to agree, first differentiability may hold, but second
order behavior is still governed by the competing active branches and may
fail to have a single Hessian. One should establish a strict root-label gap to
claim ordinary local Hessian validity.

A useful distinct-root gap diagnostic must exclude an entire row-to-root
choice, not one expanded Hungarian slot. For each row `r`, forbid assigning
that row to its currently selected root (all duplicate columns carrying that
root), solve the capacity-constrained assignment again, and take the smallest
cost increase over rows. This is the minimum excess cost among all root-label
vectors differing in at least one row. Forbidding one duplicate slot instead
can return zero by swapping equal-root columns and measures only representational
tie, not a nonsmooth branch boundary. This explains the intended meaning of
the root-banned assignment diagnostic; no gap was computed here.

The label-free objective `L_off=||offdiag(UV4U†)||_F²/16` is smooth and
independent of the assignment branch. Since `f=L_off+g`, where `g` is the
minimum diagonal mismatch, a strict assignment gap gives
`H_f=H_Loff+H_g,j0`. At a distinct-root tie, `H_Loff` remains well-defined,
but the extra lower-envelope term `g` can be nonsmooth. Thus a fixed-label
Hessian should not be reported as the original offdiagonal Hessian or as a
branch-independent Hessian at a tie.

## Meaning of a negative-curvature trigger

Suppose the point is in a strict-gap neighborhood and the fixed-label gradient
is small. A reliably negative Hessian eigenvalue below the prespecified
numerical trigger (for example, less than `-1e-6`) identifies a candidate
negative-curvature direction for that local smooth branch. For a unit
eigenvector `v`, Taylor expansion gives

`f_j0(x0+t v)=f_j0(x0)+t grad(f_j0)·v + (t²/2)lambda + o(t²)`.

At a stationary point, `lambda<0` predicts local decrease for sufficiently
small nonzero `t`; with a small but nonzero gradient, compare both signs and
verify the branch value because the linear term can dominate at finite `t`. A
strict assignment gap is not required to transfer a verified branch decrease
to the minimum objective: since `j0` is active at `x0`, `f_j0(x0)=f(x0)`, and
for any trial point `x1`, `f(x1)=min_j f_j(x1)<=f_j0(x1)`. Therefore, if
`f_j0(x1)<f_j0(x0)`, then `f(x1)<f(x0)`, even if another assignment becomes
active at `x1`. A gap is required to interpret the selected Hessian as the
ordinary Hessian of the lower envelope in a neighborhood, not for this
pointwise descent implication. Reevaluate the branch and the reassigned
minimum at the actual trial point; finite-step decreases are hypotheses about
optimization only, not diagonalization or exact certification.

A negative trigger is not a certificate that the family contains a valid
circuit, and it says nothing global about all starts or the topology. Conversely,
a computed nonnegative floating-point Hessian does not prove a local minimum:
finite precision, chart/gauge directions, a nonzero gradient, conditioning,
assignment ties, and the chosen threshold all limit that inference. Even a
positive-definite Hessian would establish only a strict local minimum for a
smooth branch under the usual regularity assumptions, not global optimality;
semidefinite curvature is weaker still. A diagnostic that does not cross its
negative threshold is therefore inconclusive, not a no-go result.

## Implementation caveats for root-gap reporting

- Report the active root vector as well as raw expanded column indices.
- State whether the constrained solve forbids a selected root for one row or
  merely one duplicate slot; only the former measures a distinct-root gap.
- Report the gap in the same normalized or unnormalized units consistently
  used for branch objectives, and account for numerical tolerance separately.
- Pair Hessian interpretation with gradient norm, symmetry/eigensolver
  residual, assignment gap, and a direct finite-step reassignment check.

No numerical diagnostic or search was run under this hypothesis.
