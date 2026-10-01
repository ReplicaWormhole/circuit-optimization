# Orbit phase ramps and an explicit twelve-CX family hypothesis

Hypothesis53; structural mathematics only, no numerical or symbolic search,
exact certificate execution or candidate promotion. Conventions: ancilla-free
four wires, q0 MSB, right cycle, chronological gates, arbitrary local gates
plus directed all-to-all CX. The accepted exact13 incumbent remains unchanged.

## Exact phase freedom for this Fourier basis

On a cycle orbit of length m, write S|j>=|j+1 mod m> and
Fm[k,j]=omega^(k j)/sqrt(m), omega=exp(2pi i/m). Let R be diagonal unitary
with entries rj. The weighted shift R S R† has successor coefficient
r(j+1) conjugate(rj). A matrix is diagonalized by Fm exactly when it is
circulant: inverse Fourier conjugation of diagonal matrices gives circulants,
and Fourier conjugation of circulants is diagonal. With the weighted shift's
single cyclic support, being circulant is equivalent to all successor
coefficients having the same value rho. Therefore

Fm R S R† Fm† diagonal iff rj=r0 rho^j with rho^m=1.

Necessity uses unit modulus and the wraparound relation; sufficiency gives
R S R†=rho S directly. Constants r0 are arbitrary unit phases. Thus singleton
phases are arbitrary, the pair ratio is ±1, and quartet ratios are1,i,-1,-i.
Constant phases are sufficient but are not necessary. For the full fixed
orbit Fourier basis, a computational-diagonal R preserves diagonalization
iff this condition holds separately on each orbit. This does not characterize
arbitrary changes of eigenbasis or arbitrary nondiagonal commuting unitaries.

For U0=F P in input coordinates, apply the criterion to phases ordered by
each original orbit. If a cheaper encoder implements Q=diag(r_address)P,
the same criterion uses each quartet address in cyclic order before Fourier.
Output bit reversal merely permutes roots. The actual orbit67 quartet labels
are[0,2,1,3]; rho=i^s changes these to[s,s+2,s+1,s+3] modulo4. The pair
labels[0,2] become[2s,2+2s] modulo4 for rho=(-1)^s. Singleton roots stay1.
There are128 discrete sector shifts(4³*2) and six independent orbit phase
constants; no enumeration or fitting was performed.

Independent verifier reviewed the weighted-shift sign, Fourier/circulant
equivalence, closure and special cases. Relative-phase Toffoli substitutions
must either obey these ramps or have explicit Fourier compensation. A truth
table alone would miss this issue. Even assuming three-CX relative-phase
replacements for all three CCX gates in the seven-CX encoder, its cost would
be16CX before mixing. Beating13 requires at least four plus the entire mixing
cost in joint savings; optimizing the high67 circuit alone is not a credible
sub13 claim.

## One competitive but unproved construction family

Use an exact parity prefix CX0→1,CX2→1,CX3→1 (threeCX), four XX/YY blocks
F(c,t;a,b)=exp(i a XX+i b YY), and one bridge CX1→2. Primitive chronological
order is:

1. The three parity CX gates;
2. F(0,2), F(1,3);
3. CX1→2;
4. F(0,1), F(2,3).

Permit arbitrary local layers between primitives and at both endpoints.
Each F compiles by the exact two-CX identity in
topology14_exact_matchgate_derivation.md. Thus every member has compiled
CX upper3+4*2+1=12. The optional free local layers around prefix gates
generally destroy a literal parity interpretation; preserving the exact
parity-prefix subfamily means inserting only the common input/output layer
around the full three-gate prefix. The family can begin with that restricted
parity version and relax deliberately under a later reserved budget.

In the separate native {CX,F(a,b),arbitrary locals} model it has four CX and
four F primitives, eight native entanglers. This is not an eight-CX cost.
The phase-ramp freedom above gives exact target variants only if an actual
orbit-address implementation is obtained; it does not establish that this
short parity prefix implements the nonlinear P or that the remaining blocks
finish diagonalization. Relative phases from any replacement must be tracked.

This is an explicit bounded candidate-family hypothesis with a genuine
twelve-CX resource ceiling, not a discovered circuit. It replaces the earlier
six-CX parity-correction prefix with a three-CX parity tag and one bridge,
retaining four exactly costed two-CX blocks. An independently reserved future
round could test one fixed primitive schedule with finite initializations
and complete gate lists, or derive its invariant subspace action analytically.
No success, topology exclusion, exact bound below13, optimality or novelty
claim follows here. Failed testing would leave other architectures and bases
open. Exact certification and promotion would require separate independent work.

Adversarial review: ramp constant condition includes wraparound; labels use
actual bit-reversed order; phase freedom separated from free circuit cost;
parity meaning scoped under extra local layers; each entangler explicitly
costed. No confirmed defect. Remaining gap is constructing any diagonalizing
member of the twelve-CX family. Report and actual chats are retained evidence.

## Separate unscheduled global-gap probe

The global interval remains6<=Cmin<=13. A concrete next mathematical probe
could freeze one lexicographically first degree4422 survivor from the existing
144class, then choose the two output Z observables with smallest backward
lightcones, using stable wire ties. For A=U†ZjU and B=U†ZkU, diagonalization
requires Hermiticity, support within those lightcones, [A,V]=[B,V]=0,
A²=B²=I, [A,B]=0 and TrA=TrB=TrAB=0. Expand in the invariant Pauli-orbit
basis to obtain exact joint polynomial constraints. This couples observables
instead of merely asking that each support admits some commuting operator.
Review existing screens first; only if demonstrably distinct reserve a single
120second, one-thread symbolic elimination/certificate attempt. Infeasibility
with a valid exact certificate would exclude only that frozen schedule, or
an explicitly proved equivalent class. An inconclusive attempt excludes
nothing. No schedule selection, symbolic computation or new screen was
executed here; all3024 six-CX survivors retain their current status.

Verifier independently confirms the actual output-order shifts and the
arithmetic12CX resource ceiling; realizability remains conjectural.

Correction to the unscheduled lower-bound proposal: coordinator identified
that commuting cycle observables require full support, so surviving schedules
typically have full four-wire lightcones. With full supports, the proposed
joint polynomial constraints omit circuit reachability and contain the known
exact13 diagonalizer's observables. That system is feasible independently of
the six-CX schedule and is therefore vacuous as an exclusion test. Do not
reserve elimination on it. A future probe must encode actual shared six-gate
template reachability or a proved commuting-through-axis constraint, beyond
support and commutant membership, before any budget is spent. No elimination
was executed and no survivor status changed.
