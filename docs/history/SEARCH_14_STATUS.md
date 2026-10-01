# Four-cycle search record: accepted exact13 incumbent

Historical chronology retained from the root-level status file. Paths in the
earlier entries were written for the former flat root. Historical computational
files now live in `archive/legacy_root/`. For the reviewed current state, see
`docs/STATUS.md`; manuscripts now live in `writing/drafts/`.

## Latest reviewed round: runs387–388

Projected24-coordinate Hessian and independent finite differences fail the
frozen negative-curvature trigger; zero escape fits. No minimum or family
exclusion follows. Exact control-wire absorption reduces the relaxed family
to an equivalent98-coordinate description with boundary remapping, not fewer
CX. Native4CX+4F compiles to12CX but no passing member is known. Exact13 and
global6<=Cmin<=13 remain. See `collaboration/ROUND_13.md` and actual chats.
Next full coupled98-coordinate diagnostic is recommended, not dispatched.


## Latest reviewed round: run386

The116-coordinate relaxed-prefix family embeds both385 checkpoints exactly,
independently checked beforefit. Both unperturbed warmstarts immediately stay
invalid; tiny added-coordinate gradients do not exclude the family. Native
4CX+4XXYY=eight entanglers compiles to12CX, but no passing member. A new
joint observable lemma is independently reviewed; it excludes no sixCX
skeleton yet. Exact13/global6<=Cmin<=13 remain. See `collaboration/ROUND_12.md`
and actual chats. Next projected-curvature diagnostic is recommended, not
dispatched, and requires its own frozen finite budget.


## Latest reviewed round: runs383–385

The orbit rewrite saved six local gates but no CX:143 total gates,67 CX.
Independent exact384 certifies the changed list. Restricted92parameter
4CX+4XXYY family compiles to12 CX; two frozen starts in385 both fail
independent checks (max off-diagonal0.943862/0.707107). No promotion or
family exclusion; exact13 and6 <= C_min <= 13 remain. See
`collaboration/ROUND_11.md` and actual researcher chats. Next round is a
newly frozen relaxation/warm-start hypothesis, not dispatched; support-only
observable elimination cannot establish a stronger lower bound.


## Latest constructive round: runs381–382

Bounded deterministic beam381 finds exact orbitP with7CX+3CCX, compiledCX
upper25, independently checked all16inputs. SharedGray phase compiler yields
mixing42 and fullnative67CX, exact_check382 passes conductor64 UV4=DU and
rootcounts6/3/4/3. This is highercost than accepted13, no incumbentpromotion
or globalboundimprovement. See `collaboration/ROUND_10.md`, native gate list,
exact output and actual chats. Next phase-aware joint simplification not dispatched.

## Latest architecture round: runs379–380

Two unordered-pair mutations fail independent checking. Run379actually used
376source despite parent377 metadata; frozenarchive preserved with explicit
correction. Corrective380 independently verifies377source path/hash before
fitting,120seconds total2x200iterationbudget one thread; actual4.57seconds.
Invalid12CX outputs maxoff0.511528 and0.890495; secondoptimizer abnormal
termination retained without hidden retry. No improvement or exclusion.
Orbit-address Fourier family and Clifford-only obstruction reviewed, but
no smaller exact circuit or newCXbound. See `collaboration/ROUND_9.md`.

## Latest curvature and residual-label rounds: runs377–378

Run377 converges47label-free iterations; Hessianleastvalue-3.34e-8 does not
cross-1e-7, so no perturbedfit executes. Independent finite differences do
not resolve the nearzero sign. Run378 swaps10/14 in a separately reserved
fixedtarget fit; maxoff0.143016 improves while offloss0.018068 worsens versus
377offloss0.009423. Both12CX lists independently invalid. Runs failed,
board33–39 terminal; no certification or incumbent promotion. Next recommended
changed unordered-pair round is not dispatched. Exact13 andglobal6<=C_min<=13
remain unchanged; optimum unproved. See `collaboration/ROUND_7_8.md`.

## Latest joint-label round: run376

Two composed targets,200L-BFGS iterations and five SVD corrections each,
120seconds total budget, one thread; actual8.64seconds. First saved12CX list
improves maxoff to0.15172779155 and offdiagonal loss to0.0094231083920;
second maxoff0.70710679183. Both independently invalid. Board29–32 terminal,
run376 failed; no new certificate or promotion. The inherited halfshift-support
consequence supplies no new survivor exclusion. See `collaboration/ROUND_6.md`
and actual chats. Next recommended label-free curvature round is not dispatched.
Exact13 and global6 <= C_min <= 13 remain unchanged; goal incomplete.

## Latest bounded spectral-label round: run375

Four fixed multiplicity-preserving label swaps around the saved12CX lead were
fitted with full168local freedom,100L-BFGS iterations and at most three SVD
corrections each,120seconds total cap, one thread. Actual14.85seconds. All
four saved lists independently fail, maxoff0.258819,0.258819,0.279401,0.279401.
Run375 is failed; board26–28 are terminal. No exact certificate was duplicated
or new incumbent promoted. See `collaboration/ROUND_5.md` and actual exchanges
in `collaboration/CHATS.md`. Next recommended joint-swap round is not dispatched.
The global interval remains6 <= C_min <= 13; optimum unresolved.

## Continued optimum goal: runs373–374

The persistent goal remains active and incomplete. No smaller exact circuit
was found in these two bounded runs; the global interval remains
`6 <= C_min <= 13`. Runs373–374 and hypotheses20–25 are terminal.
Actual chats and the explicit signed-radical 73-gate circuit are linked from
`collaboration/ROUND_3_4.md`.

Run373 tests two changed12CX schedules; maxoff0.4488 and0.6729. Run374 follows
two true13CX braid precursors through commuting local constraints and a
reviewed13-to-12 collapse. Independent complete matrices match the collapse
to5.36e-16 and3.56e-16 up to phase, but final12CX residuals0.6105 and0.6674
still fail. No topology exclusion or exact-certification target follows.

A joint transverse-input-axis lemma and full product-local CNOT centralizer
are independently reviewed. The former requires a schedule-specific forcing
argument before it can strengthen the global lower bound. The latter does
not turn bounded fitting failures into exclusions of the entire gate family.

## Current state after coordinated round2

The accepted CX incumbent is
`collaboration/work/algebraic13_new/canonical_algebraic_ansatz.json`:
13CNOTs and60local rotations with exact signed radical metadata. Acceptance
and independent evidence review: `collaboration/EXACT13_ACCEPTANCE.json` and
`collaboration/work/verifier_round2/REPORT.md`. Exact runs370–371 report256
zero intertwining entries and60local normalizations. Their full-matrix
implementation is the same; separate tower/sign/scalar and normalized100-digit
checks supply independent corroboration. No second full exact implementation
or algebraic degree32 claim is made. Preserve the exact14 baseline and proofs.

In the ancilla-free, all-to-all, arbitrary-single-qubit plus directed-CX model,
`6 <= C_min <= 13`; optimality is open. No lower-bound transfer to other
entangler models, total-spin, strong Schur or novelty claim follows.

Round2 reconciled371terminal runs and16terminal board hypotheses before
assignments17–19. Run372 tested one changed12CX pair order from two seeded
full-local starts (max200iterations each,120seconds total cap, one thread),
finishing in8.052seconds. Both native gate lists fail independent NumPy and
shared-checker validation, maxoff0.7561 and approximately1. These are bounded
failures, not topology exclusions. Algebraic work excludes only fixed-angle
terminal deletion and proposes a decorated three-wire rewrite. See
`collaboration/ROUND_2.md` for outcomes and an unscheduled next-round budget.

## Historical chronological search notes

The following sections preserve earlier snapshots, including then-current
bounds, incumbents and next-step recommendations. They are superseded by the
current-state section above; no old failed experiment or certificate is erased.

The target is an ancilla-free circuit $U$ for the four-qubit right shift
$V_4|x_0x_1x_2x_3\rangle=|x_3x_0x_1x_2\rangle$, with at most 14 CNOTs,
arbitrary one-qubit gates, and all-to-all connectivity. The objective is
achieved by `topology14_exact_matchgate_rational.json`: a 14-CNOT gate list
with rational-$\pi$ angles. `exact_check.py` and the independent
`delete14_integer_audit.py` each prove $UV_4=DU$ over the 96th cyclotomic
field, and `algebraic14_raw_matchgate_certificate.py` checks the four-block
construction and its exact equality to the gate list. Qiskit independently
simulates the gate list. The circuit does not diagonalize total spin; no
optimality claim is made. See `topology14_exact_matchgate_derivation.md`.

`experiments.sqlite3` is the authoritative attempt ledger. Run
`python3 experiment_log.py summary` and `python3 experiment_log.py audit`
before assigning new work. Every bounded search must be reserved before
execution and finished even when unsuccessful. The ledger records code and
candidate hashes. The findings below describe **restricted searches and
specific fixed operators**, not a global 15-CNOT lower bound.

## What has been tested

- Five-CNOT parity-prefix scans over the rational-$\pi/8$ phase models
  (runs 45, 61, 68) tested 824 eligible endpoints; the best completed tail
  has ten CNOTs, for 15 total. The 54 five-mask models produce four distinct
  valid 15-CNOT topologies, saved in `fivemask15_topology_0.json` through
  `_3.json` (run 100). Reordering commuting Fourier butterflies tested 2,472
  further tail syntheses and still bottomed out at 15 (run 79). This is a
  bounded Qiskit synthesis result, not an optimality proof.
- Fixed-circuit PyZX optimization, Qiskit seed/basis scans, output diagonal
  gauges, short output CNOT permutations, wire permutations, and input
  dihedral symmetries did not yield 14 (runs 44, 49, 51, 56, 59, 66, 73).
- Numerical deletion, rerouting, and local SU(2) optimization of 14-CNOT
  topologies have not passed the diagonalization check (runs 47, 50, 52,
  57, 60, 62, 64, 67, 70, 71, 74, 80, 81, 83, 113, 114, 117). The best
  reported maximum off-diagonal entry in
  the first mutation family is about 0.1267 (runs 50 and 67). A zero-angle
  serializer bug affected saved candidates for runs 81 and 83; corrected
  copies and notes are in runs 109 and 111. Other audited candidate losses
  were rechecked.
- A fixed global-parity-gauged four-wire boundary has operator Schmidt rank
  16 across the 01|23 cut for **every** phase angle: its exact realignment
  determinant is 1/256 (run 89,
  `topology14_boundary_schmidt_derivation.md`). It therefore needs at least
  four CNOTs crossing that cut and cannot replace its current four-CNOT
  implementation with three. This applies only to that fixed boundary.
- The unchanged final five-CNOT Fourier suffix cannot use four CNOTs. Its
  exact ranks are 16 across 02|13 and 03|12, but 2 across 01|23; the cut
  requirements contradict every four-edge interaction graph (runs 97 and
  102, `topology14_suffix_rank_derivation.md`). A changed intermediate basis
  or circuit split is required for a 14-CNOT solution by this route.
- Output eigenspace Givens and two-level gauges, including coordinated pairs,
  have not produced the required three-CNOT fixed boundary or a 14-CNOT
  circuit in their bounded scans (runs 76, 90, 91, 92, 93, 94, 96).
- A direct orbit-coordinate construction found two reversible quadratic
  shears that place the three length-four shift orbits in parallel affine
  planes (run 108), but a conventional separately synthesized preprocessor
  takes 21–22 CNOTs (run 115). Relative-phase shear variants still took at
  least 16 CNOTs before the Fourier block (run 118), whose chosen exact
  branch construction itself needs at least four by cut ranks (run 119).
- Uniform one-qubit Clifford input symmetries and rotation-invariant
  entangling diagonal Clifford input gauges preserve the shift, but bounded
  Qiskit scans returned at least 15 CNOTs (runs 127 and 131). These are
  compiler observations for the tested representatives and seed.
- The qubit shift and standard fermionic translation have different sign
  holonomy on a two-cycle, so no diagonal gauge alone conjugates one to the
  other (run 123). A parity-conditioned Fock Fourier basis **does** exactly
  diagonalize the qubit shift (run 125). Its stand-alone parity twist has an
  exact shortest closed parity tour of eight CNOTs in the tested phase-gadget
  model (run 132). This factorization reproduces the original eight-CNOT
  baseline prefix and ten-CNOT Fourier circuit up to a commuting gauge and
  output row relabeling/phases (run 135); it is not a new eigenspace basis.
- One reversible shear plus two CNOTs exposes conserved parity and two
  parity-controlled affine three-bit actions (run 124). A separate exact
  preprocessor takes 14 CNOTs, or 11 with relative-phase CCX gates before
  correcting orbit-dependent phases (run 129). The canonical Fourier block
  has a six-CNOT cut-rank lower bound (run 128). The same bound held for
  4,000 sampled orbit-row relabelings and phase-corrected variants (run 133).
  These bounds apply to the specified separate blocks, not jointly
  synthesized circuits or arbitrary eigenspace bases.
- A shared linear normal form for the two parity-sector affine actions uses
  two CNOTs and leaves a one-CNOT conditional discrepancy (run 136). The
  resulting controlled-affine permutation has a numerically checked
  13-CNOT Qiskit compilation (run 145), but generic monolithic synthesis of
  the subsequent orbit Fourier target is much larger. Structured joint
  synthesis remains open.
- On the exact15 final-five-CNOT suffix, all 52 single three-cycle output
  eigenrow permutations and 516 pairs across distinct eigenvalue sectors
  fail the four-CNOT cut-rank graph filter (runs 139–140). On catalog
  topology 2's earlier nine-CNOT split, four of the 52 three-cycle gauges
  retain four rank-compatible five-CNOT graphs, all among the ungauged
  eight (run 144). These are necessary rank screens for fixed prefixes;
  compatible graphs are not synthesized circuits.
- Run 146 found a near-hit with the exact six-CNOT prefix followed by four
  two-CNOT pair blocks on `(02,13,01,23)`. Two independent refinements
  (runs 148 and 150) saved 14-CNOT gate lists passing the numerical checker
  and independent Qiskit simulation at off-diagonal error near $10^{-15}$.
  At that stage this was **numerical evidence**, not yet an exact circuit
  certificate. The four pair blocks had numerically identified Weyl coordinates
  $(\pi/4,\pi/8,0)$ for `02,13` and $(\pi/8,\pi/8,0)$ for `01,23`
  (runs 149 and 151). A nullspace gauge search moved the arbitrary
  degenerate-eigenspace mixing angle while increasing the number of
  near-$\pi/16$ local angles from 18 to 64 of 108 (run 153). The rational
  construction in runs 167–171 subsequently closed the exactness gap.

## Earlier hypotheses before the exact construction

1. Exactify the numerically passing 14-CNOT topology from runs 148/150:
   combine its four rational-Weyl pair blocks with exact local frames,
   exploit common-$R_Z$ gauge freedom, and seek a rational-$\pi$ gate list.
   The pairwise block structure is in run 149 and the nullspace gauge search
   in run 153.
2. Jointly change the intermediate eigenbasis and the suffix topology. An
   earlier split with nine prefix CNOTs and a five-CNOT suffix passes the
   necessary cut-rank filter (run 107), although the first local optimizations
   failed (runs 113–114). Moving the split earlier and changing several
   interaction edges is still open.
3. Search multiple-edge 14-CNOT topology changes with stable SU(2)
   parametrization and label-free loss, using distinct 15-CNOT seeds from
   the topology catalog. This differs from the exhausted single-delete and
   single-reroute family.
4. Interleave reversible orbit-coordinate operations with the orbit Fourier
   blocks, rather than synthesizing a two-shear preprocessor and controlled
   Fourier block separately. The exact branch actions are recorded in runs
   112, 116, and 119.
5. Use the parity-sector affine maps from run 124 directly: seek a shared
   short CNOT conjugator and interleave its Fourier rotations with the shear
   instead of paying for a separate orbit Fourier block.

The run-148/150 topology was exactified in runs 167–171. Its full rational
gate list passes `check_circuit.py`, two independent exact gate-list checks,
an exact four-block identity check, and Qiskit. Earlier numerical residuals
were useful search evidence; the exact certificates establish the result.

## First search below 14 CNOTs

- Run 172 deleted each of the eight CNOTs in the four matchgate blocks of
  the exact 14-CNOT circuit and refit arbitrary local SU(2) gates once per
  topology. None passed; the best 13-CNOT numerical candidate has maximum
  off-diagonal entry 0.30348. This is a bounded negative result, not a
  lower bound.
- Runs 173 and 175 tested output-CNOT gauges on the fixed terminal pair
  blocks. Exact magic-basis invariants rule out a one-CNOT resynthesis of
  either isolated block. Exact operator-Schmidt ranks rule out a three-CNOT
  synthesis of the fixed terminal four-qubit suffix after any single output
  CNOT. These statements do not cover a shifted synthesis boundary.
- Run 174 proves a global four-CNOT lower bound by a spectral-path
  obstruction for any local Pauli axis whose image stays one-body. See
  `global_lower_bound4_spectral_path.md`. The established optimum interval
  is therefore `4 <= C_min <= 14`.

The next construction hypotheses at that stage were a joint seven-CNOT
realization after moving the synthesis boundary across the center local
layer and a first-layer matchgate, and rerouting a neighboring matchgate
edge after deleting one CNOT. Bounded versions of both were subsequently
tested below.

## Subsequent 13-CNOT tests and stronger global bound

- Runs 176 and 177 tested all 12 single-output-CNOT gauges and 132
  nontrivial ordered pairs of output CNOTs with fixed-seed Qiskit
  recompilation. Every valid result retained at least 14 CNOTs. Run 181
  tested all eight invertible binary circulant input symmetries of the
  four-cycle, using shortest CNOT preparations and the same compiler;
  again no 13-CNOT result. Run 183 applied three PyZX modes directly to
  the exact 14-CNOT circuit; basic optimization retained 14, full
  optimization did not support its non-Clifford+T angles, and graph
  extraction gave 45. These are bounded compiler observations.
- Runs 178, 179, and 186 fitted several seven-CNOT joint-tail and
  delete/reroute topologies with arbitrary local SU(2) layers. None
  passed the diagonalization checker. Run 184 used exact modular
  operator-Schmidt ranks to exclude a five-CNOT resynthesis of a fixed
  shifted-boundary suffix under eight cross-pair output-CNOT gauges.
  These results do not exclude other intermediate eigenbases or global
  topologies.
- Runs 185, 187, and 188 establish the global lower bound
  `C_min >= 5`: an exact spectral-path Gram calculation excludes a final
  CNOT on adjacent physical wires in a hypothetical four-CNOT circuit,
  while the opposite-pair case conflicts with the `(10, 6)` eigenspace
  multiplicities of `V4^2`. See `global_lower_bound5_spectral_gram.md`.
  The established interval is `5 <= C_min <= 14`; optimality remains open.
- Run 189 tested three altered five-CNOT parity prefixes followed by the
  exact 14-CNOT circuit's eight-CNOT matchgate-tail topology, with two
  independent local fits per prefix. No candidate passed; the best
  13-CNOT gate list has maximum off-diagonal entry 0.40710. The next
  construction hypothesis is a five-CNOT prefix whose endpoint mismatch
  lies on the first matchgate pair `(0,2)` or `(1,3)`, so the first pair
  block might absorb the mismatch during joint synthesis.
- Run 190 classified all 7,776 chronological five-CNOT interaction-pair
  schedules. Degree and full-support conditions leave 1,008 schedules;
  the terminal spectral-path lemma excludes 664 of these, leaving 344.
  See `global_lower_bound6_terminal_filter.md`. This is a restricted
  necessary-condition screen and does not establish a six-CNOT bound.

## Later bounded searches below 14

- Run 191 tested all 1,452 noncancelling three-output-CNOT eigenrow
  gauges with fixed-seed Qiskit recompilation; the minimum remained 14.
  Run 192 fitted three substantially rewired 13-CNOT suffix topologies;
  none passed. Runs 214 and 217 used a continuous Ising entangler in
  place of one CNOT and followed its strength toward zero while
  reoptimizing local gates; no valid 13-CNOT endpoint was found.
- Run 193 proves that an arbitrary same-pair two-qubit diagonal output
  gauge cannot reduce an isolated terminal XX/YY matchgate from two
  CNOTs to one. Runs 197, 201, and 202 show that every nonlocal
  **single** cross-pair parity phase, at any real angle, makes the fixed
  two-pair terminal suffix require at least five CNOTs. Run 205 found
  the same lower bound for 1,296 sampled two-mask gauges. These bounds
  apply to that fixed suffix; products of masks at arbitrary angles and
  changed boundaries remain open.
- Runs 194, 196, 198, and 200 exhaust the 62 catalogued four- and
  five-mask rational phase models against five-CNOT prefixes. No
  eligible prefix has a first-pair endpoint that absorbs the omitted
  mask locally. Runs 204, 206, 210, 213, and 216 then rule out moving
  that omitted phase solely into the fixed second layer, into a
  disjoint 03/12 reroute, or into any four-CNOT final layer: the
  transported target has balanced-cut ranks `(5,16,16)` and needs at
  least six CNOTs as a fixed suffix. These are phase-model and
  fixed-boundary obstructions, not a global 13-CNOT lower bound.
- Runs 209 and 211 shift the boundary to the full eight-CNOT tail.
  Several two-mask gauges pass necessary seven-CNOT cut-graph screens,
  but bounded fits on two compatible topologies failed. Run 219
  screens simple same-eigenvalue row swaps and Givens rotations: none
  allows a three-CNOT fixed final layer, while all retain some
  seven-CNOT-compatible full-tail graph. Compatibility is only
  necessary; synthesis remains open.
- Runs 195, 199, 203, 207, 208, 212, 215, and 218 sharpen the
  five-CNOT topology filter with exact spectral-path and rank-one
  certificates. Of the 344 schedules surviving run 190, 88 remain
  unexcluded. Exact selected-axis witnesses survive in two of the
  remaining classes, so the global lower bound is still five.

The current constructive priorities are a changed intermediate
eigenbasis or first-layer operator with a seven-CNOT tail, and a
multi-entangler homotopy that can pass the single-entangler loss barrier.
The lower-bound priority is a joint constraint on multiple commuting
input observables for the 88 surviving five-CNOT schedules.

## Current lower bound and later experiments

The preceding chronology records what was known at each stage. Run 223
proves the global bound `C_min >= 6` for ancilla-free circuits with arbitrary
one-qubit gates and all-to-all CNOTs. A half-shift trace obstruction closes
the opposite-pair cases left by the earlier spectral-path tests. The exact
rank-one certificates close the other five-CNOT chain shapes; the exhaustive
schedule checker accounts for all 7,776 pair schedules. An independent agent
rebuilt the spectral Gram matrices and audited the finite classification.
See `global_lower_bound6_complete_filter.md`. The current established
interval is `6 <= C_min <= 14`; the optimum remains open.

Runs 220 and 222 tested a subset of output eigenrow gauges with a full
seven-CNOT tail. Several graphs passed necessary rank screens, but bounded
local fits found no 13-CNOT circuit. Run 221's two-entangler homotopy likewise
found no valid 13-CNOT endpoint. These are limited searches, not a global
lower bound of 14.

## Subsequent six-CNOT screens and 13-CNOT searches

- Run 224 enumerated all 46,656 ordered unoriented six-CNOT pair schedules
  and applied the already proved degree, full-support, and selected-axis
  tail obstructions. Exactly 3,576 schedules survived. Runs 229 and 231
  proved that every shift eigenbasis contains a vector of Schmidt rank at
  least three across each adjacent balanced cut, so each such cut needs at
  least two crossing CNOTs. All 3,576 schedules already meet that test.
- Runs 237 and 240 established a stronger exact spectral-path obstruction:
  no nonzero traceless Hermitian image supported on an adjacent physical
  pair can occur. Run 242 applied it to the six-CNOT schedules and excluded
  368 more, leaving **3,208**. An independent agent rebuilt the 120-by-120
  Gram matrix using integer arithmetic, reran the rational Gröbner
  certificate, and independently reproduced the schedule counts. See
  `global_lower_bound7_adjacent_pair_support.md`. This is still a
  necessary-condition reduction, not a seven-CNOT lower bound.
- Runs 243, 245, and 246 examined opposite-pair and three-wire Pauli-image
  spans for degree-three wires. The exact 24-dimensional three-wire Gram
  matrices have nullities 45, 40, and 64 for the three ordered physical
  shapes. Their nullspaces do not by themselves exclude an entire shape;
  the normalized rank-one and involution systems remain open. No additional
  six-CNOT schedule was excluded. See
  `global_lower_bound7_degree3_centered_gram.md`.
- Runs 225, 227, 232, 233, 235, 241, and 244 tested paired and phased
  eigenspace gauges, interleaved seven-CNOT tails, and distinct five-CNOT
  prefix endpoints. Their exact fixed-operator rank results and bounded
  numerical fits found no valid 13-CNOT diagonalizer. See
  `docs/history/search13/search13_joint_gauge_summary.md`; numerical failures do not rule out
  other topologies or gauges.
- Runs 226 and 228 tested 12 rank-compatible changed matchgate block
  graphs and four refinements. Runs 230, 234, and 238 deleted each prefix
  CNOT, perturbed the closest near-hit, and rerouted nearby prefix edges.
  All checked 13-CNOT candidates failed. The closest candidate by maximum
  off-diagonal entry was about 0.267 (run 230); the perturbations returned
  to the same local loss. These are bounded searches, not no-go proofs.

The established global interval remains `6 <= C_min <= 14`. Current next
routes are a three-wire selected-axis certificate for the remaining
six-CNOT schedules and a 13-CNOT construction with a five-CNOT prefix
whose endpoint stays close to the exact six-CNOT prefix while its first
pair layer is jointly synthesized.

## Later seven-tail gauges and six-CNOT schedule classes

- Runs 247, 249, 250, 255, and 257 screened altered five-CNOT prefixes,
  cross-pair reroutes, and an adaptive beam of 169 distinct 13-CNOT
  schedules. Bounded local fits found no valid circuit; the best maximum
  off-diagonal entry remained about 0.267. Exact cut-rank tests allowed
  all 32 directed reroutes of the fixed residual, so they did not exclude
  that family. See `docs/history/search13/search13_endpoint_guided_summary.md` and
  `search13_adaptive_topology_result.json`.
- Run 251 derived an exact modular commutant necessary condition for the
  residual after the fixed six-CNOT prefix. It excludes 75,012 of 279,936
  ordered unoriented seven-CNOT tail pair schedules for that prefix,
  leaving 204,924. This is a fixed-prefix screen, not a global 13-CNOT
  obstruction. Runs 254 and 258 explored continuous output gauges and
  fitted two compatible seven-CNOT tail schedules; none gave a valid
  13-CNOT diagonalizer.
- Runs 264–266 turned the numerical gauge lead into an exact
  same-eigenvalue unitary gauge over `Q(zeta_96)`. Its gauged full tail has
  operator-Schmidt rank eight across cut 03|12, exactly, while still
  diagonalizing the shift. Thus a gauge-independent rank-greater-than-eight
  obstruction on that cut is false. It remains an eight-CNOT tail; a
  seven-CNOT synthesis is open. See
  `search13_fulltail_invariant_exact_witness.py` and its result JSON.
- Runs 248, 253, 256, and 259–263 studied selected images on six-CNOT
  schedules. Exact one-axis and commuting two-axis witnesses satisfy the
  tested involution and spectral-path conditions for one representative,
  but no common six-CNOT circuit was constructed. Run 262 partitions the
  3,208 surviving schedules into 401 size-eight physical dihedral orbits
  (159 degree-3333 classes). No further schedules were excluded. See
  `global_lower_bound7_joint_aa_pair.md` and
  `global_lower_bound7_dihedral_orbits_result.json`.

The established global interval is still `6 <= C_min <= 14`; the optimum
is open. The exact rank-eight output gauge and the 401 orbit representatives
are concrete starting points for the next constructive and lower-bound
experiments.

## Exact mixed-cut screens and broader gauge tests

- Runs 267 and 270 applied the exact rank-eight output gauge to seven-CNOT
  tails after the fixed six-CNOT prefix. Ordinary cut ranks and the exact
  reverse-commutant test leave 109,944 of 279,936 ordered unoriented pair
  schedules. Three numerical fits failed independent validation. Run 273's
  PyZX reductions of the exact 14-CNOT gate list also found no shorter
  valid list; the basic rewrite retained 14 CNOTs and passed an exact
  checker.
- Runs 278, 280, 283, 284, and 286 developed mixed input/output tensor
  cuts. An exact 5-by-5 minor first excluded one promising seven-CNOT
  order. Across the 124 schedules obtained by deleting one CNOT from the
  exact eight-CNOT tail and changing at most one remaining pair, ordinary
  cuts exclude 94, mixed wire-bond cuts 22 more, and splitting each CNOT
  into a rank-two control-target bond excludes the final eight. This is
  an exact no-go for that **fixed target and 124-schedule neighborhood**;
  see `docs/history/search13/search13_rank8_neighborhood_exact_obstruction.md`.
- Runs 287–289 extended the exact modular-rank and CNOT-factor mincut
  screen to all 279,936 seven-CNOT pair schedules for this fixed output
  gauge and prefix. Ordinary cuts exclude 152,256; the earlier commutant
  excludes 17,736 more; all 127 complementary mixed cuts exclude another
  108,368. **1,576 schedules survive necessary tests**, saved in
  `search13_rank8_global_allcuts_result.json`. An independent agent
  reproduced the count and all modular ranks with a separate max-flow
  implementation. No valid seven-CNOT tail follows from survival.
- Runs 285 and 291 fitted selected candidate 13-CNOT topologies. The
  interrupted run 285 is retained as a superseded numerical artifact;
  run 291 fitted six distance-three survivors with process, direct, and
  structured warm starts. All checked 13-CNOT circuits failed
  independent validation.
- Run 279 (after a superseded run 276) expanded the exact five-CNOT
  parity-prefix search to *arbitrary continuous* diagonal input phases
  constant on shift orbits.
  Integer-character certificates exclude all 312 first-pair-local
  endpoint/mask patterns. This applies to parity-`Rz` prefixes with the
  fixed first-pair blocks, not arbitrary local gates. Run 282 tested one
  distinct nonmonomial five-CNOT prefix with a jointly fitted eight-CNOT
  tail; its two fits returned to the previous invalid local basin.
  Run 292 jointly optimized output gauges toward rank eight on both
  balanced cuts 01 and 03. Three bounded starts reached their iteration
  caps without a near-zero simultaneous solution; this numerical failure
  does not exclude another gauge.
- Runs 268, 271, 272, 275, 277, 281, and 290 sharpened one six-CNOT
  schedule's selected-axis test. For a fixed canonical first-axis image,
  exact trace moments rule
  out every commuting second-axis image in the full 12-dimensional
  spectral-path/involution family. An exact sparse scan found that all
  signed equal-weight involutions with at most three Pauli terms in the
  first-axis span are canonical symmetry variants. Unequal-weight or
  denser first axes remain open. **No six-CNOT schedule was excluded;
  3,208 remain.** See `global_lower_bound7_joint_aa_family_trace.md`.

The best valid circuit still has 14 CNOTs, and the proved global lower
bound is six. The optimum remains open. The next constructive route is
to change the output gauge or prefix beyond the fixed rank-eight tail;
the next lower-bound route is to classify dense first-axis images in the
surviving six-CNOT schedule classes.

## New exact identity-gauge screen and proof manuscript (runs 293--310)

- The standalone six-CNOT lower-bound draft is
  `four_qubit_cycle_lower_bound_paper.tex` (compiled PDF alongside it).
  It consolidates the exact spectral-path, half-shift trace, algebraic
  certificate, and exhaustive five-CNOT schedule arguments. A separate
  agent reran the five Gram generators and four certificate/classifier
  scripts and found no mathematical gap conditional on those certificates;
  this is internal review, not external peer review.
- Runs 293 and 297 numerically optimized an output eigenspace gauge to
  make two mixed cuts rank eight. Run 300's numerical full-cut profile
  leaves 1,428 seven-CNOT pair schedules. This was only candidate
  selection. Runs 305--306 found a simpler **exact** fact: the original
  eight-CNOT tail with the identity output gauge has rank four at both
  decisive cuts and an exact rank profile at all 127 complementary mixed
  cuts. Run 309's exact necessary screen leaves the same 1,428 schedules
  out of 279,936; survival does not give a seven-CNOT implementation.
  Runs 302 and 310 fitted selected survivors against the numerical and
  exact identity targets, respectively; all checked 13-CNOT candidates
  failed independent validation. These bounded fits are not no-go proofs.
- Runs 295 and 299 excluded simultaneous rank-eight behavior for fixed
  exact-G output gauges formed by disjoint row swaps or one/two
  complex-Hadamard row blocks. Runs 303 and 307 tested the larger
  within-eigenspace row-permutation family: all 622,080 permutations
  have modular rank at least 11 at mask 83, so none gives rank eight for
  that fixed G. These restricted-gauge results do not apply to the
  identity gauge.
- Run 298 exactly excluded all 27 three-CNOT pair orders for a fixed
  three-qubit middle block of the 14-CNOT circuit. Runs 294 and 304
  classified real selected-axis images with at most four Pauli terms
  for one six-CNOT survivor; any feasible first image there needs at
  least five terms. No complete six-CNOT schedule was excluded.

The established interval remains `6 <= C_min <= 14`; the optimum is open.

## Output-row gauges and expanded seven-tail neighborhoods (runs 311--321)

- Run 311 completed the pairwise-regrouping check begun in run 298 for
  the exact 14-CNOT matchgate tail. The four overlapping three-wire
  middle blocks are each excluded from a three-CNOT synthesis by exact
  mixed-cut ranks over all 27 three-pair schedules for that fixed block.
  This does not exclude
  a jointly changed intermediate basis or a general 13-CNOT circuit.
- Runs 312--315 and 317, 320 screened the 124 seven-pair tail schedules
  obtained by deleting one of the exact eight tail CNOTs and rewiring at
  most one retained pair. Exact modular mixed-cut ranks exceed the
  CNOT-factor mincut bound after every output row gauge in
  $\mathrm{GL}(4,2)$ (20,160 maps), every single output-row swap, and
  every pair of disjoint output-row swaps. Thus none of these
  gauge/schedule combinations synthesizes its **fixed target** with
  arbitrary one-qubit gates. See
  `docs/history/search13/search13_output_row_gauge_neighborhood.md` for the exact scope.
- Run 316 allowed two retained-pair rewires. Of 1,632 schedules across
  all 20,160 linear output gauges, 6,096 gauge/schedule pairs survive
  the necessary exact rank screen, representing only 13 distinct pair
  orders, all at distance two. Survival is not a circuit construction.
  Run 319 fitted the nine newly revealed schedules with their shortest
  compatible output gauges, and run 321 fitted 18 more diverse
  identity-gauge survivors. No checked candidate passed the independent
  diagonalization test. Those are bounded numerical failures, not
  exclusions of the surviving schedules.
- Runs 322--323 extended the finite output gauge to all 1,120 oriented
  three-cycles of output rows. None of the 124 one-rewire tail schedules
  survives the exact mixed-cut screen. With two rewires, 416
  gauge/schedule pairs survive over 104 gauges, but their four distinct
  pair orders were already among the identity-gauge survivors. This
  screen adds no new compatible pair order in that neighborhood.
- Run 324 optimized the full continuous output gauge inside the four
  shift-eigenspace sectors using block Procrustes updates alternating
  with seven-CNOT local fits. Four bounded trials on two schedules
  failed independent circuit verification; the best maximum
  off-diagonal error was about 0.534. See
  `docs/history/search13/search13_procrustes_gauge_summary.md`. This numerical failure does
  not rule out a continuous gauge or a 13-CNOT circuit.
- Run 326 fitted the two three-cycle gauges with the lowest selected
  mixed-cut-rank sums against representative pair orders 506 and 922,
  from deletion-warm and random starts. All four checked 13-CNOT
  candidates failed, with best maximum off-diagonal error about 0.992.
  The rank score is a heuristic; the fit does not exclude either order.
- Run 325 changed both the prefix and suffix topology of the exact
  14-CNOT circuit: one prefix deletion, one other prefix reroute, and
  both CNOTs of one matchgate pair rerouted. Four joint 168-parameter
  fits with continuous output-gauge continuation failed independent
  validation; the best 13-CNOT candidate had maximum off-diagonal
  error about 0.387. Run 327 deepened that candidate three times from
  perturbations of size 0, 0.02, and 0.08 for up to 800 iterations;
  all returned to the same invalid loss, about 0.078733. This records
  one local basin, not a no-go for the topology. See
  `docs/history/search13/search13_joint_prefix_orbit_summary.md` and
  `search13_joint_prefix_orbit_refine_result.json`.

- Run 328 screened 143 single directed-CNOT replacements around the
  joint-prefix candidate and fitted six replacements on distinct slots.
  All six failed validation; the best maximum off-diagonal error was
  about 0.3864, with loss 0.078797. The raw replacement score did not
  predict a successful fit. See `docs/history/search13/search13_joint_prefix_one_step_escape_summary.md`.

- Run 330 fitted 17 diverse one-slot changes independently of raw-loss
  ranking: a disjoint-pair replacement at every slot and four direction
  reversals. All failed validation, but slot 8 changed from 0->3 to 1->2
  reached a new basin: loss 0.06019279537 and max offdiagonal error
  0.25491651759. Run 332 refined that topology from four starts up to
  1000 iterations; three returned to loss 0.060192795058 and error
  about 0.25491354, while the largest perturbation was worse. No
  validated 13-CNOT circuit was produced.
- Run 329 stopped on a finite-difference sanity check before fitting.
  Run 331 repeated with a smaller finite-difference step and computed
  the full 168-parameter Hessian at the earlier run327 basin. No
  eigenvalue was below -1e-6 (smallest about -3.67e-8, gradient norm
  1.13e-7). Six signed weak-positive-curvature escapes failed; the
  best returned to the old invalid basin. This is numerical local
  evidence, not a certified minimum or topology exclusion.

- Run 333 changed the objective while retaining the cycle eigenspaces:
  H_t = Re(V) + t Im(V), with t = 1/4, 4, -1/4, -4. Each Hermitian
  target has four distinct eigenvalues. Four weighted fits followed by
  direct-cycle polish generated eight invalid candidates; all polished
  fits returned to loss about 0.060192795058. The equivalence of exact
  diagonalizers does not imply small residuals transfer between targets.
- Run 334 stopped before fitting on a duplicate-proposal guard. Run 335
  corrected the selection and fitted twenty new two-slot changes from
  the run332 source. No raw-loss ranking was used. All candidates failed;
  the best returned to the source basin. Complete proposals are frozen
  for replay, since exclusion depends on the contemporaneous matching
  search13*.json files and does not cover every historical topology.
- Run 336 used large seeded local-angle perturbations (Gaussian std0.8)
  on the two lowest fitted-loss run335 topologies. All four candidates
  failed independent checks and had worse losses, about 0.328--0.453.
  These bounded trials do not exclude their topologies.

- Run 337 broadened the constructive search beyond the exact14-derived
  prefix/suffix neighborhood. Two complete thirteen-CNOT pair orders
  (round-robin and chain) were fitted with two direction variants and
  two independent Gaussian local-axis starts each. No local chunks
  were inherited. All eight candidates were invalid, but the chain
  seed33821 reached loss0.001215834838 and maxoff0.04303923. Run339
  deepened it from three starts; all converged to loss0.001215834228
  and maxoff about0.0430346. This is the best archived thirteen-CNOT
  numerical lead, not a validated upper bound. See
  `docs/history/search13/search13_fresh_schedules_summary.md`.
- Run338 tested positive output-row weights on the two older source
  topologies, followed by ordinary-cycle polish. All sixteen stage
  candidates failed; polishing returned to the older0.0601928 basin.
  The full weight vectors, source hashes and schedules are saved.
- Run340 rerouted each of the thirteen slots in the new chain source
  to edge03 with alternating directions. All13 full-local fits failed
  and were worse than the source. These are one-start tests per
  specified topology, not exclusions of all-to-all circuits.
- Run341 tested twelve-CNOT circuits by deleting either the first or
  last CNOT of the new chain lead, with warm and noisy starts. All four
  failed independent checks; the best maxoff was about0.5688734. No
  smaller exact circuit or new lower bound resulted.

- Run342 inferred nearest root labels on the fresh chain source and
  checked the correct multiplicities (6,3,4,3). Two automatic-Jacobian
  Gauss-Newton solves of U V = D U converged to nonzero squared residual
  about0.001216502002; both13CNOT candidates failed. Maximum offdiagonal
  entry was about0.04299887, without improving the ordinary-loss basin.
- Run343 used three dynamic positive-row-weight updates from each of
  three starts on the fresh chain source. All12 stage/final candidates
  failed. A transient stage reached maxoff0.03305318 at a worse ordinary
  loss; all final ordinary-loss polishes returned to0.001215834228.
  The transient weighted point is not a valid circuit or an improved
  ordinary-loss basin. Full stage vectors are retained.
- Run344 smoothly relocated one interaction to edge03 or02 through
  powered-CNOT products. Intermediate gates are outside the counted
  model. All four native endpoint/polish candidates had13CNOTs and
  failed; both paths ended worse than the source.
- Run345 tested fresh full13 ring and star pair orders, two direction
  variants and two independent starts each, no inherited chunks.
  All eight candidates failed. One optimizer line search terminated
  abnormally, explicitly recorded. These are two unordered architecture
  families, not four; no topology is excluded by these bounded fits.

- Run346 appended one of three output-CNOT basis permutations to the
  fresh chain source, deleted earlier slot0,4 or8, and refitted all
  local layers. All nine native13CNOT candidates failed and were worse
  than the source. Output permutation alone preserves diagonalization;
  deleting an interaction is not assumed to preserve it.
- Run347 proposed65 one-slot other-pair mutations, excluded only the13
  known run340 closing-edge unordered schedules, and short-fit all52
  remaining proposals for40 iterations. The best six by fitted loss
  were deepened for400 iterations. All58 saved candidates failed
  independent checks. Ranking and frozen schedules were independently
  verified. Proposal0 changes first edge01 to02; its maxoff0.03519944
  is lower than the ordinary chain source's, but its ordinary loss
  0.001455266622 is worse. It is an alternative numerical starting
  point, not a valid circuit or a better ordinary-loss basin.

- Run348 refined the alternative opening02 circuit from four perturbation
  scales (0, .05, .2, .6), with 650 iterations per fit. All four native
  13-CNOT candidates failed independent checks. The first three return
  to loss 0.001455266622, so these starts did not escape that basin.
- Run349 tested all eleven interior deletions from the original chain
  source, short-fit for 80 iterations, then deepened the best three for
  500 iterations. All fourteen native 12-CNOT candidates failed. Best
  maximum off-diagonal entry is 0.25881905382. Root independently
  reconstructed all matrices and checked deletion schedules and ranking.
  Together with run341, every single deletion slot has now been tested
  from this source; this does not exclude other starts or 12-CNOT circuits.

The best validated circuit still has 14 CNOTs and the proved global
lower bound is six. The optimum remains open. Next mutate one interaction
of the alternative opening02 source, excluding saved original-chain
one-slot schedules. This explores two-slot changes relative to the original
chain; rank by short fitted ordinary loss before deeper optimization.

- Run350 computed the full numerical Hessian at the fresh chain source
  and tested six signed-pi escapes along low positive-curvature modes.
  All six native 13-CNOT candidates failed independent checks; four
  returned to the source ordinary-loss basin and two reached a worse one.
  Floating-point curvature is not a local or global minimum certificate.
- Run352 computed a numerical polar spectral-sector correction and
  retained its full256-term Pauli generator. Independent expm reconstruction
  diagonalizes the numerical source to about1e-14. Its largest coefficients
  have weight3; weight4 contributions also occur. This correction has no
  assigned CNOT cost and is not a valid improvement. Gauge freedom prevents
  interpreting its support as a universal synthesis obstruction.
- Run351 was interrupted after partial short fits. The frozen60 schedules
  and saved partial candidates are retained; recovery must inspect actual
  process status before relaunching. Its ledger running flag alone is not
  evidence of a live computation. No complete tournament result exists yet.

- Run353 tested twelve fresh all-six-pair schedules with independent Haar
  SU(2) local starts,120 short iterations, and500 iterations for the three
  best fitted losses. All fifteen native13-CNOT candidates independently
  failed; best ordinary loss0.240446671 is worse than the chain lead.
  Frozen exclusions, seed conversion, saved schedules and ranking passed
  independent checks. Run351 recovery is recorded separately as run354.

- Run354 recovered the seven preserved run351 endpoints and fitted the
  remaining53 frozen schedules, then deepened six selected by fitted loss.
  All66 endpoints (59 new gate lists) independently failed. Best maximum
  off-diagonal entry is0.10922791848, best ordinary loss0.009117300939,
  worse than both source leads. Root checked all saved matrices, schedules,
  ranking and recovered-file hashes. Recovered optimizer metadata is
  explicitly unavailable; those seven entries are not claimed as new fits.
  Checkpoint hash/config binding is a known resilience limitation for a
  future interrupted resume, not evidence of an error in this completed run.

Bounds remain6–14. The original chain lead is the best ordinary-loss source.
New searches should change architecture substantially or exploit spectral
structure; repeating these same local fits is not supported by these results.

## New thirteen-CNOT numerical success (runs355,358,359)

A targeted output-eigenlabel swap at rows10 and14 escaped the old chain
basin. Run355 reached1.69e-7 maximum off-diagonal entry; run358 reached
4.74e-8. One truncated-SVD Newton step in run359 reached1.997e-15 and
passed the numerical checker with13 native CNOTs. Independent NumPy
UV-DU error is1.58e-15 and Qiskit off-diagonal error1.85e-15; Qiskit
transpilation retains13 CNOTs. Source/evidence are in
`collaboration/work/root/labelswap_svd_refine_candidate.json` and its
result/review files. This is numerical evidence, not an exact certificate.
The exact14 incumbent and proved6 lower bound remain authoritative;
board hypothesis5 requests exact recognition/gauge fixing of this source.
The first run355 summary claiming no better basin was incorrect and is
explicitly corrected in the review. Saved numerical measurements are intact.

Runs356–357 establish a general opposite-pair support obstruction for
input-axis images. Exact rational nullspace constraints force swap symmetry;
traceless Hermitian involutions then contradict the proper-support halfshift
trace identity. Root independently checked all six squared-difference
functionals against all45 rational nullvectors and reviewed the trace proof.
The necessary six-CNOT screen excludes184 additional schedules, leaving3024.
This does not establish a seven-CNOT lower bound or circuit witnesses.

Run361 froze108 Euler coordinates to the pi/8 grid through12 bounded
fits, preserving numerical13-CNOT diagonalization (maxoff2.24e-14).
The all-coordinate pi/8 rounded list fails numerical and exact checks;
no exact13 certificate resulted. Root independently checked both lists.
Original valid source is preserved. A finer-grid or analytic construction
requires a distinct bounded hypothesis; no failed rounding excludes13.

Runs362–364: pi/24 gauge freezing prescribed108 Euler coordinates and
preserved numerical13-CNOT diagonalization, but fullgrid rounding fails
(maxoff0.255 and exact conductor96 check). Root's canonical commuting-axis
conversion replaces168 local Euler coordinates by64 physical rotation
coordinates (60 after removing final outputZ phases), retaining13 CNOTs.
Matrix equivalence upglobalphase error1.12e-15 and independent Qiskit
maxoff2.27e-15 pass; generic angles remain. This representation provides an
analytic-recognition source, not an exact certificate. Run363 tests all13
singleCNOT deletions from the passing source: all16 saved short/deep
12-CNOT candidates independently fail. Best maxoff0.2588190370.

Next distinct constructive direction: exploit collective input rotations
M tensor4, which commute with the cycle, to remove three intrinsic source
parameters before recognizing the reduced canonical circuit. This has not
been tested in these runs. No passing exact13 or12 construction is claimed.
