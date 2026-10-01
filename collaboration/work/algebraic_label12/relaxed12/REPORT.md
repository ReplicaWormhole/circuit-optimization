# Actual-template operator reachability and relaxed twelve-CX scope

Hypothesis58; structural derivation only, no numerical or symbolic search.
Four wires, q0 MSB, chronological directed all-to-all CX, arbitrary locals,
no ancillas. No incumbent certification duplicated or candidate promoted.

## Operator-Schmidt necessary bounds

For a2|2 cut S|Sc, expand a Hermitian operator A in tensor Pauli bases.
Its real coefficient matrix is16by16; its rank is operator Schmidt rank rS.
Local gates or any unitary acting entirely on one side preserve this rank.
A crossing CX has operator Schmidt rank2:
CX=|0><0| tensor I + |1><1| tensor X, with orientation adapted to the cut.
Expanding CX A CX† gives at most four products per original Schmidt term,
so rS_after<=min(16,4rS_before). Initially an output local Pauli has rank1.
With s crossing gates, the generic cap is min(16,4^s).

If the observable is still supported on one single wire immediately before
a crossing CX, the improved cap2 applies: a control axis aX+bY+cZ becomes
(aX+bY) tensor X+cZ tensor I, and a target axis becomes
I tensor aX+Z tensor(bY+cZ). The other wires are identity. Importantly,
earlier same-side entanglers may spread the observable within that side while
keeping cut rank1; then the single-wire cap2 need not hold, and generic4
must be used. A cut-only gate count is insufficient to invoke that improvement.

For A=U†ZjU in any diagonalizer, [A,V]=0. Cyclic permutation maps the cut
01|23 to12|30 (equivalently03|12), so their Schmidt ranks must agree.
The opposite cut02|13 maps to itself. Combine these equalities with actual
schedule-specific upper caps; no necessary positive lower rank beyond1 has
been proved. Indeed Z0Z1Z2Z3 is a traceless Hermitian involution commuting
with V, with full support and rank1 across every cut. It is reachable from
a local axis by a parity CX network. Thus full-support commutant constraints
or rank caps alone cannot establish a universal rank obstruction.

## Nonvacuous actual-template test

Freeze one six-CX skeleton and one output Zj. Build its backward-evolved
Pauli coefficients using the literal reversed gate sequence. Each local
conjugation acts as a real SO3 matrix on that wire's X,Y,Z coefficients and
fixes I; impose R^T R=I and detR=1. Each CX acts as an exact signed permutation
on Pauli strings. For example Xc→XcXt, Yc→YcXt, Zc→Zc, Xt→Xt,
Yt→ZcYt, Zt→ZcZt; products retain Pauli multiplication signs. This recursion
encodes actual shared-template reachability, not arbitrary operators with
matching support. It automatically preserves Hermiticity, A²=I and trace.

Impose cycle-invariance of the resulting coefficient vector as exact linear
equalities in those coefficients. The coefficient entries are polynomials in
the SO3 variables, so this gives a concrete necessary feasibility system for
the frozen skeleton. Schmidt minors/cut caps are extra derived constraints
or useful diagnostics, not replacements for the recursion. Feasibility of
one observable is not sufficient for four commuting output axes or a full
diagonalizer. An exact infeasibility certificate would exclude that skeleton;
no such certificate, symbolic elimination or schedule selection was attempted.
All3024 six-CX survivors and global interval6<=Cmin<=13 remain unchanged.
Any future symbolic run needs a separately reserved finite resource budget
and review of actual SO3/CX recursion before execution.

## Relaxed family concurrence

Literal chronology is L0 CX01 A CX21 B CX31 L1 F02 L2 F13 L3 CX12
L4 F01 L5 F23 L6, where each local layer is four arbitrary SU2 gates and
Fct(a,b)=exp(i a XX+i b YY). Nine local layers give108coordinates and four
F blocks eight angles, total116. A=B=identity embeds each run385 checkpoint
exactly by inserting24zero coordinates while retaining all84local and8F
coordinates. This is a385 source embedding, not an assumed377 mapping.
Verifier should check the full native and compiled matrices before fitting.
Free A/B remove the earlier literal parity-prefix interpretation; no fixed
physical Hamming-parity tag is assumed. Native cost4CX+4F=8entanglers,
compiled upper12CX via the existing exact two-CX identity per F. No known
diagonalizing member or family exclusion follows from failed restricted385.

Adversarial review: strengthened firstcross bound correctly requires single
wire support; globalZ defeats rank-only lower claims; actual-template SO3
variables retained rather than arbitrary commutant coefficients; zeroA/B
mapping and116count explicit. No confirmed defect. Main unresolved tasks
are a genuine reachability obstruction and a passing relaxed-family member.
Report/chat files are intentional retained evidence; coordinator owns commit.

## Independently reviewed joint necessity

Coordinator proposed a stronger joint statement, independently reviewed by
this researcher and verifier. Rank1 across all three2|2 cuts implies complete
tensor factorization: write the coefficient tensor as psi01 tensor phi23;
its Schmidt rank across02|13 is rank(psi01)*rank(phi23), so both factors
must themselves factor. A Hermitian involution then has normalized local
factors that are I or unit Pauli axes. Cyclic invariance forces proportional
local factors; balance (trace zero) excludes I. Thus such an observable is
±Q tensor Q tensor Q tensor Q, with Q a unit Hermitian Pauli axis.

For two axes Q,R, QR has eigenvalues exp(±i theta), where cos(theta) is
their Bloch-vector inner product. Global commutation requires (QR)^tensor4
to equal (RQ)^tensor4. Tensor eigenphase sums include±2theta, so exp(4i theta)=1;
theta in[0,pi] gives parallel, antiparallel or orthogonal axes. Conversely
those cases commute globally. Three mutually orthogonal axes give global
observables with product relation Q^tensor4 R^tensor4=S^tensor4, since
(±i)^4=1. Therefore at most two independent commuting cycle-invariant fully
factorized observables exist, independence taken modulo scalar identity.

The four true U†ZjU observables are independent: every nonempty product is
unitarily conjugate to a nontrivial output Z string and has trace zero, so
cannot equal a scalar identity. Consequently at least two output images must
have operator Schmidt rank>1 across at least one2|2 cut. This is a proved
nonvacuous joint necessity beyond the individual globalZ counterexample.
It still excludes no six-CX skeleton until actual template reachability is
shown to force too many fully factorized images; no such forcing was proved.
