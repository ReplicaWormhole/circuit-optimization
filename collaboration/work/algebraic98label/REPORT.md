# Spectrum-assignment residual for the full98 family

Hypothesis85; read-only derivation for the proposed full98 polish. No matrix
evaluation, numerical or symbolic search, optimization, derivative evaluation,
or certificate attempt was performed here. Let `U` be any circuit from the
unitary 98-coordinate family, let `V=V4`, and let `B=U V U†`.

## Exact spectrum and admissible diagonals

The computational strings split into two fixed cycles, one length-two cycle,
and three length-four cycles under right cyclic shift. The fixed strings are
`0000` and `1111`. The length-two cycle is `{0101,1010}`. Its two Fourier
vectors have eigenvalues `+1` and `-1`; each length-four cycle contributes one
copy of each fourth root of unity. Therefore the exact multiplicities are

`+1: 6,  -1: 4,  +i: 3,  -i: 3`.

Let `D` range over diagonal matrices whose sixteen entries are exactly this
multiset, in arbitrary row order. This allows every eigenvalue ordering and
every orthonormal basis choice inside a degenerate eigenspace; it does not
assign labels to particular computational rows in advance.

## Objective and zero-set equivalence

Define

`L_spec(U) = min_D ||U V - D U||_F^2 / 16`.

Because `U` is unitary, right multiplication by `U†` preserves the Frobenius
norm, so

`L_spec(U) = min_D ||B-D||_F^2 / 16`.

For each fixed `D`,

`||B-D||_F^2 = sum_(r != s) |B_rs|^2 + sum_r |B_rr-d_r|^2`.

Thus, if the existing label-free objective is

`L_off(U)=||offdiag(B)||_F^2/16`,

then

`L_spec(U) = L_off(U) + (1/16) min_D sum_r |B_rr-d_r|^2`.

The added term is nonnegative, hence `L_spec >= L_off`. Both objectives have
the same exact zero set for unitary `U`: `L_spec=0` iff some admissible `D`
satisfies `UV=DU`, equivalently `UVU†=D`, which is exactly a diagonalizer
with the correct spectrum. Conversely, if `U` diagonalizes `V`, similarity
preserves the spectrum, so its diagonal target is in the admissible set and
`L_spec=0`. Also, `L_off=0` alone already forces the same spectrum because
`B` is unitarily similar to `V`; the assignment objective changes the
off-solution landscape, not the mathematical target set.

## Hungarian costs as row overlaps

For a candidate diagonal entry `b_r=B_rr` and root `lambda`, an assignment
cost is

`C(r,lambda)=|b_r-lambda|^2`.

The minimum-cost capacity-constrained assignment with capacities `(6,4,3,3)`
is a standard 16-by-16 Hungarian problem with repeated root columns (or an
equivalent capacitated assignment). Since `|b_r|^2+1` is independent of the
assigned root, the same assignment minimizes `-2 Re(conj(lambda)*b_r)`.

There is also a row-overlap form. Both row vectors `U[r,:]` and `(UV)[r,:]`
have norm one, and

`<U[r,:],(UV)[r,:]> = (U V U†)[r,r] = b_r`.

Therefore

`sum_s |(UV)[r,s]-lambda U[r,s]|^2 = 2-2 Re(conj(lambda)*b_r)`.

Summing over rows and dividing by16 gives the fixed-`D` objective. This
explains why a raw row residual and the diagonal cost produce the same
assignment: the off-diagonal row contribution is independent of the label.

## Gradient and ties

For a fixed selected `D`, write `E=UV-DU`. The smooth branch gradient with
respect to a real circuit parameter `theta` is

`d L_D/d theta = (2/16) Re tr(E† ( (dU/dtheta)V - D(dU/dtheta) ))`.

When the minimizing label vector is unique, the selected branch remains
active in a neighborhood and its gradient is the derivative of `L_spec`.
The minimum of finitely many smooth branches is continuous and piecewise
smooth. At a tie between distinct label vectors, a classical gradient may
not exist; the gradient obtained while holding one selected `D` fixed is
only that branch's gradient. Permutations of duplicate columns carrying the
same root do not change `D` or the branch; relevant ties exchange different
roots between rows.

For reproducible optimization, keep the primary floating costs unmodified,
fix root-column order `(+1,-1,+i,-i)` with capacities `(6,4,3,3)`, pin the
assignment-library version, define a deterministic secondary rule for equal
cost optima, and record the selected label vector and assignment changes.
Avoid adding an arbitrary epsilon to the costs: it can alter the true
minimum at near-ties. Recompute the actual minimum-assignment objective at
every line-search point while differentiating the currently selected
fixed-`D` branch. Branch switches can make quasi-Newton smoothness assumptions
fail and may invalidate gradient checks close to ties. Final candidates must
still be judged by the original off-diagonal loss and independent gate-list
checks; no numerical objective is an exact certificate.
