# Exact13 acceptance and separate native-cost review — hypothesis18

Researcher: verifier_round2. Coordinator persisted returned review because
child shell calls failed. Researcher directly read HANDOFF.md and collaboration
README, then reviewed coordinator-transmitted exact source excerpts and audit
summaries. This was not independent access to the complete verifier file, and
no new exact execution or duplicate same-verifier reproduction was launched.
Findings and gate-model distinctions were exchanged directly with both peers.

Recommendation: ACCEPT exact13 in arbitrary single-qubit plus directed-CX model.
Source canonical_algebraic_ansatz.json SHA256:
d500ff2a79c4717ff38d935d71b67bd798b06a7c6899a275fec25570b86b8d13.
Source has73 chronological gates,60local rotations and13CX.

Reviewed arithmetic: Fraction coefficients in Gaussian Q(sqrt2), exact positive
square tests and signed trig parsing, dictionary multiplication with XOR masks
and generator-square products on AND masks. Empty dictionaries map to exact
zero even if generators are dependent. Do not assert degree32 or independence.
Chronological left row updates use MSB bit masks; residual column permutation
rightrotate(j) gives precisely UV4=DU with labels(1,i,-1,-i).

Local exact signed trig values satisfy s²+c²=1. Strict radicand bounds establish
realness; A=(1+c)I-i sP satisfies A*A=2(1+c)I>0 when c!=-1. PureP handles
c=-1. Nonzero projective factors and any negative/global scalar transfer the
identity to normalized unitary rotations. Current source hash agrees with
both existing370/371 exact result files, each256empty residual dictionaries,
60local normalization checks and13CX. Separate tower relation and sign/scalar
audits and100digit normalized simulation provide independent corroboration.
Runs370/371 reproduce the SAME full-matrix implementation; no second full exact
matrix implementation is claimed. No blocker identified in transmitted evidence.

Gate models:
- Arbitrary one-qubit+directedCX: accepted circuit13native entanglers,13compiledCX.
- Same gates plus arbitrary RZZ(theta)=exp(-i theta Z tensor Z/2): unchanged
  exact13 remains13native/13compiledCX; improved regrouping not established.
- Existing independently tunable XX/YY+CX hybrid: older exact14-derived baseline
  has10native entanglers, depth8, compiledCX upper14. Native JSON is numerically
  matched to exact compiled source with analytic provenance, not a direct exact
  native JSON certificate or ten-CX result.
- Arbitrary SU(4) elementary gates are a different model; exact13 trivially
  supplies upper13, no optimized count certified.

No lower bound transfers automatically to other entangler models. No improved
exact native regrouping was established. Preserve all old certificates.
