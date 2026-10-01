# Round22 independent static review

Board141, verifier10. Read updated AGENTS.md, README.md, docs/STATUS.md,
collaboration/README.md and algebraic11nearhit/EXACT11_PROPOSAL.md. The accepted
source is run411's30-gate11CX word, SHA557c55cddefb697bcd4e038db300117d0aaaf6f6af97e851e0831baf257dbe4a.
This collaborator prepared and reviewed source/metadata without executing
matrices, RNG, AD, an optimizer or an exact checker. Root later executed the
reserved415 preflight and417 endpoint audit; their runtime results are
reviewed read-only in REPORT.md.

## Algebraic10 suffix identity and narrow obstruction

With s=sin(pi/8), c=cos(pi/8), B=sX+cZ and Nprime=-cX+sZ,
independent Pauli multiplication gives B Nprime=-iY and
B Z Nprime=(s^2-c^2)X+2scZ=(Z-X)/sqrt2=A. Nprime(v;y)
commutes through T_x and CH(u;x), whose wires are disjoint. The resulting
target-only suffix sectors are F00=I,F10=Z,F01=-iY,F11=A, retaining phases.

For chronological u- then v-controlled target interactions with arbitrary
unconditional target locals, U11=U01 U00^-1 U10; reversing control order
reverses the last two factors. If first three sector matrices are monomial,
their product is monomial. Left sectorwise target monomial output gauges
retain this property for I,Z,-iY; no such gauge can make A monomial because
all four A entries are nonzero. Thus no two-interaction target-only suffix
in that fixed computational-sector architecture gives the required bases.

The fixed-prefix branch operators also agree independently: K00=CZxy,
K10=(-i)^y Zx,K01=-(i^x)Zy and
K11=[(i+1)I/2-(i-1)Xy/2]Zx, since A Z A=-X. The x=1 part of uv00
is nondegenerate and forces the target monomial restriction there despite
x=0 degeneracy. This is a restricted fixed-prefix/unchanged-CH/target-only
obstruction; control-label permutations, changed earlier gates, general
degenerate-eigenspace mixing and the full numerical10 family are excluded
from its scope. It is not a10CX lower bound or family exclusion.

## Numerical10 preparation review

The frozen numerical10delete family has11local layers, four wires and three
free SU2 coordinates each,132total, around10fixedCX. All11chronological
deletions and every retained direction/order are represented. Source collapse,
principal SU2/global-phase conversion, U3 singular chart and Torch sinc/
Pauli rotation formula were inspected. All initial perturbations are allocated
first from PCG64seed10012201, eleven132coordinate normal draws sigma0.025.
Initial/best complete54-gate lists, bases and histories are preserved.

Fit guard checks exact canonical config bytes against ledger, input hashes,
accepted source, code snapshot, dependencies, owner numerical10/board140,
parent411/targetV4 and one-shot marker. Case8-second stops and total105-
second BaseException deadline are distinct, with preserved best words/vectors.
Max80iterations/max4000calls/no restarts and105internal/115external seconds
are within the stated finite budget. No additional confirmed static defect was found after correcting the target
guard as recorded below. Reserved415 subsequently passed independent matrix/
serialization/gradient consistency; static review alone was not validation.

## Independent preflight provenance and ownership

New verifier10/preflight adapts the earlier independent verifier11router
primitive/SciPy-expm implementation. Its independent references do not import
shared checkers or family matrix helpers; numerical10's copied family functions
are audit subjects. It checks121U builds across11cases,11AD and22FD losses.
Syntax passed. Numerical10 separately reviewed its source/config/counts and
reported no blocker without running matrices.

The coordinator accidentally claimed142, and board-ledger ownership prevented
run413 linking before computation. That frozen metadata is preserved in
verifier10/superseded142. Run414 subsequently caught a historical target guard
string cycle versus actual CLI targetV4 before the execution marker, RNG,
matrices or AD. The unsupported string was a confirmed metadata/guard defect;
both numerical10's fit and this preflight were corrected and refrozen.
Actual reserveCLI choices V4/V4+S2 were checked directly. The144 packet is
preserved inverifier10/superseded144. Fresh145 is owned by verifier10 with
matching ledgeragent,targetV4,parent411; no budget expansion or silent retry
occurred. Current code/config/PLAN hashes are in preflight/MANIFEST.json.
Root owns reservation,
dispatch, endpoint-audit budget, closure and commit. All previous directories,
shared checker code and incumbent11 are preserved.

Adversarial closeout found no unresolved source/output issue. Reserved415 and
417 passed their finite consistency scopes. The remaining limitation is the
deliberately narrow algebraic obstruction and finite numerical diagnostic,
which establish no10CX exclusion or optimality claim.
Source packets are retained evidence; caches are ignored. No remote action,
dependency installation or commit was performed by this collaborator.
