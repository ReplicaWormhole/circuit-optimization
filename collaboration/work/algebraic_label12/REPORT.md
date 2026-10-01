# Spectral-label escape: structural guidance

Hypothesis26; structural math and workflow, no optimization or certification
run. Gate model: ancilla-free, four wires, arbitrary one-qubit gates and
directed all-to-all CX; candidates in this round have12CX. Qubit0 is most
significant, gates chronological, and V is the right cycle. The source is
`collaboration/work/numerical12/deletion_slot1_deep.json`. It is an invalid
numerical lead, not an exact circuit.

## Exact multiplicities

V permutes computational basis strings. Thus trace(V^k) counts strings fixed
by rotation by k: 16 for k=0, two constant strings for k=1 or3, and four
period-dividing-two strings for k=2. Since V^4=I and V is unitary, its
eigenvalues belong to {1,i,-1,-i}. For label j in {0,1,2,3}, the spectral
projector is (1/4) sum_{k=0}^3 i^(-jk) V^k. Its trace gives multiplicities
(6,3,4,3). In particular every valid target D_rr=i^label_r must have those
counts. Nearest-root classification without this constraint need not be a
valid spectrum.

The assignment may be reordered arbitrarily. Within each root eigenspace,
every orthonormal basis is permitted. If U0 has target D0, every diagonalizer
has the form U=W P U0, where P is an output permutation and W commutes with
D=P D0 P†. Conversely every such unitary diagonalizes V. The block sizes of W
are6,3,4,3. This parametrizes mathematical freedom, not a claim that arbitrary
P or W is free in the specified circuit gate model.

## Deterministic finite target proposal

For unitary U, put B=UVU†. Right multiplication by U† preserves Frobenius
norm, hence ||UV-DU||_F²=||B-D||_F². The latter is
sum_{r!=s}|B_rs|² + sum_r |B_rr-D_rr|². At a fixed U, choosing the closest
admissible D therefore reduces to a finite capacity-constrained assignment
with label capacities(6,3,4,3), costs c(r,j)=|B_rr-i^j|². Resolve equal-cost
assignments by lexicographic order of their sixteen-label vector. This is a
selection rule to be executed only inside the numerical researcher's reserved
run, not an algebraic search already performed here.

Starting from that assignment, consider pairs a<b with distinct labels. Swap
their labels, retaining all other labels. The normalized target-loss increase
is exactly

delta(a,b)=2 Re[(conj(D_aa)-conj(D_bb))(B_aa-B_bb)]/16.

My initial suggestion was to rank these pairs by delta, then a, then b. The
coordinator instead froze the actual round's rule to score
|B_ab|²+|B_ba|², rank descending then a then b, and choose four unequal-label
pairs. Each target gets100L-BFGS iterations followed by at most three SVD
corrections, with120seconds total and one thread. The larger coupling rule
targets a conflicting offdiagonal block; it too is a heuristic. At an
exact closest assignment, all delta are nonnegative. Floating roundoff may
give tiny negative values, so this property is a diagnostic, not a numerical
proof. Swapping equal labels does not create a distinct target and should be
excluded. Neither ranking proves that a useful escape is among its four
selected swaps. The initial ranking suggestion was not executed.

The earlier thirteen-CX breakthrough used a targeted rows10/14 label swap,
followed by refinement and truncated-SVD Newton. It motivates a different
initialization mechanism here but supplies neither twelve-CX success nor a
preferred pair for this different saved matrix.

## Scope and independent review

The numerical researcher owns actual assignment selection, fitting, ledger
reservation and serialized lists. Report both fixed-target residual and the
label-free offdiagonal residual, as well as gate count and unitarity. A small
fixed-target residual is numerical evidence only. Failed finite frozen-label
fits exclude neither other assignments nor degenerate eigenspace bases nor
the twelve-CX topology. Independent exact certification remains necessary
before any new incumbent promotion. The current13CX exact incumbent and
global interval6<=Cmin<=13 are unchanged by this derivation.

Verifier's direct reply independently checks the character traces and Fourier
multiplicities and agrees on keeping both objective reports. Its independent
matrix/list checks are a separate task. No exact-certification work was
duplicated here.

Adversarial review against the repository protocol: no confirmed defect.
Checked signs by expanding both swapped diagonal costs; checked spectral
projectors against k=0 and total multiplicity16; kept gate-model cost separate
from arbitrary output-basis freedom. Remaining risk is heuristic target
selection, not the multiplicity or norm identity. This report and actual chat
record are intentional retained artifacts; coordinator owns shared commit.

## Unscheduled next-round review

Run375's four targets failed numerical validation. Coordinator asked whether
the joint swap(2,6)(3,7) would merely duplicate free local output relabeling.
With binary row convention q0q1q2q3, this permutation flips q1 only when q0=0
and q2=1, independent of q3. It is a negative-control Toffoli tensored with
identity on q3, rather than a whole-wire X. Appending it has no demonstrated
free CX cost. Selecting its target D is a new initialization proposal, not
appending that gate to a saved circuit.

There is a stronger target-level obstruction to free local equivalence.
From the saved base labels in numerical_label12/targets.json, partial trace
of D over wires0,2,3 has eigenvalues {1+i,1-i}. The joint target has
{3-i,-1+i}; the two individual targets(2,6) and(3,7) have respectively
{3+i,-1-i} and{1-i,1+i}. For every tensor-product output unitary W,
Tr_other(W D W†)=W1 Tr_other(D) W1†, so these eigenvalues are invariant.
The joint target is therefore not related to the base or either individual
target by any product-local output conjugation. This distinguishes it from
those three targets, not every possible circuit symmetry or all four tested
targets. A new reserved joint-target fit could be justified; no such fit was
performed and no different-basin or success claim follows.
Verifier independently reproduced all four reduced spectra by direct
Gaussian-integer label sums and confirmed the binary conditional action.
