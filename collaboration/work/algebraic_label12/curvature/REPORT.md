# Bounded label-free curvature audit: structural interpretation

Hypothesis33; structural math/workflow only. No optimizer, Hessian computation,
symbolic screening or certification was performed by this researcher. Source
is the invalid twelve-CX first joint target saved by run376. The numerical
researcher owns one reserved round:100label-free refinement iterations, one
full Hessian, then only if least eigenvalue<-1e-7 two signed Euclidean-norm
0.1 perturbations and150iterations each,120seconds total and one thread.
Full168axis coordinates remain free. Ancilla-free all-to-all directed CX and
arbitrary local gates define the model; q0 is most significant.

Let f(x)=||offdiag(U(x) V U(x)†)||_F²/16. This objective has no fixed output
eigenlabels; zero is the diagonalizer criterion for unitary U. Previous frozen
target fitting minimizes a different function, so its stopping criterion says
nothing by itself about grad f. Record gradient norm at the actual auditpoint.
The least eigenpair of the symmetric Hessian H should include its normalization,
symmetry error and residual ||Hv-lambda v||. The threshold is a numerical
policy, not an exact sign certificate.

If gradient vanishes and H has a strictly negative eigenvalue, that point is
not a local minimum: Taylor expansion along v has decreasing quadratic term.
Calling it a stationary saddle requires stationarity and some increasing
direction or other conventional saddle evidence; a maximum can also have
negative curvature. At a nonstationary point, a negative eigenvalue is simply
a negative second derivative in that coordinate direction. Both signed
perturbations are useful because the linear term can dominate at finite step.
No fit failure establishes that every circuit on the topology has positive
residual. A positive semidefinite Hessian is only a second-order necessary
condition for a minimum, not sufficient in the presence of zero modes; even
a strict local minimum cannot prove a global or topology-wide lower bound.

There are exact objective symmetries with continuous tangent directions.
Output product Rz gates conjugate B=UVU† by a diagonal product unitary and
preserve every offdiagonal entry magnitude. Collective input rotations
G tensor G tensor G tensor G commute with the cycle and leave B unchanged.
They can be absorbed into first/last local layers without adding CX. These
four output-axis and three collective-input directions can overlap or become
coordinate-degenerate; no generic rank is asserted. For a symmetry path
x(t), 0=d²f/dt²=v†Hv+grad f dot x''(0). Only at stationarity does the second
term disappear. Thus treating an arbitrary gauge tangent as a Hessian zero
mode at a nonstationary point would be incorrect. Axis-angle coordinate
singularities supply another reason small eigenvalues need careful interpretation.

Verifier's independent finite differences along the saved v check
[f(x+h v)-2f(x)+f(x-h v)]/h² at two declared h values. They test directional
curvature numerically without retraining; they do not prove exact negativity.
Compare scaling and precision sensitivity rather than claiming agreement from
a single step. All final gate lists still need independent chronological
matrix evaluation, count and unitarity checks; any passing circuit would need
new independent exact certification before promotion.

If this schedule repeatedly reaches the same bounded basin, an unscheduled
analytic next route is to inspect the dominant residual pair(10,14) reported
for the improved source, not infer a forced angle from its size. Current labels
there are0 and2. A target swap preserves root capacities; it is a distinct
frozen-label mechanism from label-free curvature fitting. However, the adjacent
pair(11,15) has equal labels0, so a paired-swap target would be identical to
the single(10,14) target. This must not be sold as an additional distinct label
target. The resulting q1 reduced target spectrum returns to{1-i,1+i}, equal
to the original pre-joint base up to order; the earlier single-wire invariant
alone therefore does not prove it is locally inequivalent to that base.
An eigenspace-structured initializer could instead exploit arbitrary orthonormal
bases within root multiplicities(6,3,4,3), but requires a concrete circuit
synthesis or newly reserved search before any cost or success claim.

Adversarial review: no confirmed defect; symmetry acceleration term retained,
stationarity separated from curvature, PSD condition kept necessary-only,
and proposed paired target correctly recognized as identical at label level.
No optimum or twelve-CX exclusion follows. The exact13 incumbent persists.
Report and actual chat record are intentional retained artifacts.
