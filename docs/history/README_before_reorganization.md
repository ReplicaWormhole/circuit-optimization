# Four-qubit cycle circuit search

Archived copy of the former full README. Its commands and paths were written
for execution from the former flat repository root. Historical computational
files now live in `archive/legacy_root/`; prefix their former root paths with
that directory. Mathematical notes now live in `docs/research/`. The current
entry point is `README.md`, and the current reviewed summary is
`docs/STATUS.md`.

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


## Accepted exact 13-CNOT result

The CX incumbent is
`collaboration/work/algebraic13_new/canonical_algebraic_ansatz.json`: thirteen
CNOTs and sixty local rotations, with exact signed radical metadata. The
acceptance record is `collaboration/EXACT13_ACCEPTANCE.json`, reviewed in
`collaboration/work/verifier_round2/REPORT.md`. Runs370–371 certify all256
entries of `UV4-DU` as zero with Fraction arithmetic. Independent tower,
sign/scalar and normalized100-digit checks support the result; the isolated
full-matrix reproduction uses the same implementation, not a second verifier.
Floating `theta` fields alone do not define the exact circuit.

In the ancilla-free, all-to-all, arbitrary-one-qubit plus directed-CX model,
`6 <= C_min <= 13`; optimality remains open. This does not establish a strong
Schur transform, total-spin diagonalization or novelty. The six-CNOT lower
bound is preserved from earlier proof work, not re-proved this round.

`collaboration/ROUND_2.md` records the reconciled bounded round. One changed
12-CNOT topology with two starts failed numerical validation; no topology
exclusion or lower bound follows. Native gate models remain separate: the
older ten-entangler CX+independently-tunable-XX/YY baseline compiles to at most
fourteen CNOTs and has no direct exact native-JSON certificate.

The continued optimum goal is active. `collaboration/ROUND_3_4.md` records two
further bounded runs, neither producing a smaller diagonalizer. Actual recorded
researcher exchanges are in `collaboration/CHATS.md`. The complete readable
exact gate prescription is in
`collaboration/work/verifier_round3/EXPLICIT_CIRCUIT.md`, including signed
radicals, all 73 chronological gates, and the output-label root mapping.

The latest bounded spectral-label round is in `collaboration/ROUND_5.md`.
Run375 tests four frozen unequal-label swaps on the best saved12CX lead;
all four lists fail independent checking. The exact13 incumbent and interval
6 <= C_min <= 13 remain unchanged. A joint-swap follow-up is recommended,
with a new finite budget required before execution.

`collaboration/ROUND_6.md` records the subsequent composed-swap run376. It
improves an invalid12CX lead to maxoff0.151728; neither output passes. The
new full gate list is `collaboration/work/numerical_label12/joint/target0.json`.
No new exact incumbent or stronger global lower bound follows. A bounded
label-free curvature experiment on this new source is recommended next.

`collaboration/ROUND_7_8.md` records runs377–378: no curvature escape passes
the frozen threshold, and a dominant-pair label swap lowers maxoff to0.143016
while worsening total offdiagonal loss. All new lists remain invalid. Separate
metric leads, actual chats, and a bounded unordered-pair mutation recommendation
are preserved; no smaller exact circuit or stronger global bound is claimed.

`collaboration/ROUND_9.md` records two changed interaction-pair hypotheses,
the independently detected run379 source-provenance error and corrective380
with preflight source verification. Both corrected12CX outputs fail. An explicit
orbit-address/Fourier family and a Clifford-only obstruction are reviewed;
neither improves the exact13 incumbent or global lower bound.

`collaboration/ROUND_10.md` constructs an exact orbit encoder with7CX+3CCX
(compiled upper25), then exports a complete67CX rational-pi orbit/Fourier
diagonalizer, independently exact-certified in run382. It is a higher-cost
research baseline for joint simplification, not a replacement of exact13.

## Preserved exact 14-CNOT baseline

`topology14_exact_matchgate_rational.json` is an ancilla-free exact
14-CNOT diagonalizer of the four-qubit right shift. Its chronological gate
list contains only rational-π one-qubit rotations and CNOTs. The construction
and two-CNOT XX/YY block identity are in
`topology14_exact_matchgate_derivation.md`; the generator is
`topology14_exact_matchgate_candidate.py`. The count improves the paper's
18-CNOT explicit circuit by four CNOTs and the earlier result below by one.

```bash
python3 check_circuit.py topology14_exact_matchgate_rational.json
python3 exact_check.py topology14_exact_matchgate_rational.json
python3 delete14_integer_audit.py topology14_exact_matchgate_rational.json
python3 algebraic14_raw_matchgate_certificate.py
python3 qiskit_crosscheck.py topology14_exact_matchgate_rational.json
```

The two gate-list exact checks independently prove `U V4 = D U` in the
96th cyclotomic field; the third exact certificate independently assembles
the four XX/YY blocks and matches the gate list up to global phase. Qiskit
simulates the original gate list with maximum off-diagonal error
`4.83e-16`. Its output eigenvalue multiplicities are `(6, 3, 4, 3)` for
`(1, i, -1, -i)`. This circuit diagonalizes the cycle, but does not
diagonalize total spin and is not claimed optimal. Ledger runs 167–171
record the exact construction and audits.

The current global lower bound is six CNOTs, proved in
`global_lower_bound6_complete_filter.md` by an exact half-shift trace
obstruction, spectral-path Gram certificates, and an exhaustive
five-CNOT schedule classification (ledger run 223). Its intermediate
five-CNOT theorem is in `global_lower_bound5_spectral_gram.md`.
Thus the established interval is `6 <= C_min <= 13`; the optimum is open.
Run359 supplied the numerical source for the subsequently accepted exact
radical construction in runs366–371. Other failed fits and fixed-suffix
screens retain their restricted scopes.

## Historical numerical source of the exact13 construction

The passing gate list is
`collaboration/work/root/labelswap_svd_refine_candidate.json` (run359).
`collaboration/work/root/canonical13_stripped_candidate.json` (run364)
expresses it with60 single-axis rotation angles and13 CNOTs by pushing
commuting local rotations forward and removing four final diagonal phases.
Independent matrix comparison and Qiskit preserve numerical diagonalization.
Neither list has an exact certificate. Pi/8 and pi/24 gauge-freezing attempts
(runs361–362) preserved numerical solutions but full grid rounding failed.
All single-CNOT deletions of run359 were tested in run363 and failed their
bounded local fits; this does not exclude12-CNOT circuits.

A general opposite-pair support obstruction (runs356–357) now excludes184
additional six-CNOT schedules, leaving3024 necessary survivors. This does not
raise the global lower bound. See
`global_lower_bound7_opposite_pair_support_certificate.md`.

## Exact 15-CNOT result

`topology16_15_exact_candidate.json` is an ancilla-free, exact 15-CNOT
diagonalizer of the four-qubit right shift. It uses only one-qubit rotations
by rational multiples of pi and CNOTs. The six-CNOT parity prefix uses a
shift-invariant diagonal gauge; two modified Fourier butterflies then use
two CNOTs each. `topology16_15_exact.py` generates the gate list, and
`topology16_15_derivation.md` explains the construction. This improves the
paper's 18-CNOT explicit circuit by three CNOTs. No optimality claim follows.

```bash
python3 check_circuit.py topology16_15_exact_candidate.json
python3 exact_check.py topology16_15_exact_candidate.json
python3 delete16_integer_audit.py topology16_15_exact_candidate.json
python3 qiskit_crosscheck.py topology16_15_exact_candidate.json
```

The two exact checks independently establish `U V4 = D U` over the
32nd cyclotomic field; the Qiskit check is numerical. The output eigenvalue
multiplicities are `(6, 3, 4, 3)` for `(1, i, -1, -i)`. The circuit does not
diagonalize total spin, so it is a cycle diagonalizer, not a verified strong
Schur transform. Ledger runs 36, 38, 39, and 40 record the discovery,
prefix certificate, exact circuit, and independent integer-ring audit.

## Exact 17-CNOT result

`algebraic_boundary_candidate.json` is an exact 17-CNOT diagonalizer of the
four-qubit right shift. The gate derivation is in
`algebraic_boundary_derivation.md`; ledger run 6 records the construction.
The construction removes a full-parity phase commuting with the shift and
fuses the remaining parity correction with the first Fourier CZ into an
eight-CNOT closed parity network. The rest of the Fourier circuit uses nine
CNOTs. This improves the paper's 18-CNOT upper bound; it does not establish
that 17 is optimal.

```bash
python3 check_circuit.py algebraic_boundary_candidate.json
python3 exact_check.py algebraic_boundary_candidate.json
python3 delete_search_exact_audit.py algebraic_boundary_candidate.json
python3 qiskit_crosscheck.py algebraic_boundary_candidate.json
```

The first command is a numerical screen. The next two independently check
`U V4 = D U` with exact cyclotomic arithmetic and count 17 CNOTs; ledger runs
9 and 10 record those checks. Qiskit independently reports maximum
off-diagonal error `3.93e-16`. The output eigenvalue multiplicities are
`(6, 3, 4, 3)` for `(1, i, -1, -i)`. The circuit does not diagonalize total
spin, so no strong Schur-transform claim is made.

The target was an **ancilla-free, exact** diagonalizer of the four-qubit
right shift with **at most 17 CNOTs**, arbitrary one-qubit unitaries, and
all-to-all connectivity. The shift is

`V4 |x0 x1 x2 x3> = |x3 x0 x1 x2>`.

This is the task behind the 18-CNOT entry in Table 1 of *Exact ancilla-free
diagonalization of qubit cycles*. The paper's 18 is an upper bound from its
mixed-radix construction (10 CNOTs for the fermionic Fourier part and 8 for
the parity correction), not a lower bound. Its optimality result covers only
three qubits. The other local paper studies a three-qubit Schur sampling
circuit. A circuit that diagonalizes `V4` is enough for a destructive four-cycle
measurement with classical outcome labels, but may not satisfy all conditions
of a strong four-qubit Schur transform.

## Candidate format and check

Write a JSON file with chronological gates, for example:

```json
{
  "n": 4,
  "gates": [
    {"gate": "h", "qubit": 0},
    {"gate": "cx", "control": 0, "target": 3},
    {"gate": "u3", "qubit": 1, "theta": "pi/2", "phi": 0, "lam": "-pi/4"}
  ]
}
```

Supported one-qubit gates: `h`, `x`, `s`, `sdg`, `rx`, `ry`, `rz`, and `u3`.
Angles may be real numbers (radians) or simple multiples of pi such as
`-3*pi/4`. The only two-qubit gate is `cx`. There are no implicit wire swaps.
The example is a format illustration, not a valid diagonalizer.

Run `python3 check_circuit.py candidate.json`. The checker composes the full
unitary and reports CNOT count, unitarity error, maximum off-diagonal entry of
`U V4 U†`, fourth-root eigenvalue error, and the output eigenvalue for each
computational-basis result. It also reports whether the total-spin Casimir is
diagonal in the same basis; this is a necessary check for the stronger Schur
task, but does not by itself establish a strong Schur transform. It accepts
**any** eigenvalue order and **any** orthonormal basis inside each degenerate
eigenspace. A zero exit code means
numerical validation passed at the specified tolerance; it is not a proof of
exactness. `python3 -m unittest test_check_circuit.py` checks the checker.

The paper's baseline is in `baseline_18.json`, generated by
`python3 build_baseline.py` from its Appendix A and B gate formulas. It passes
the checker with 18 CNOTs. Optional Qiskit cross-checking and transpilation:

```bash
python3 -m pip install qiskit
python3 qiskit_crosscheck.py baseline_18.json --output optimized.json
python3 check_circuit.py optimized.json
```

Qiskit uses the reverse numerical qubit index, which the adapter handles.
Qiskit's standard optimization level 3 retains 18 CNOTs for this baseline
(Qiskit 2.5.2); this is an observation, not an optimality result.

## Search protocol

1. Use the checked 18-CNOT baseline. Its gate list calibrates qubit order,
   multiplication order, Fourier labels, and fermionic signs.
2. Give separate agents bounded approaches: simplify the complete Fourier
   and parity circuit; remove each of the 18 CNOTs in turn and reoptimize
   nearby one-qubit gates; search other 17-CNOT topologies while using the
   full freedom inside degenerate eigenspaces. Each agent submits a JSON gate
   list and a derivation, never only a circuit drawing or a low loss value.
3. Run the checker on every proposal and rank only passing gate lists by CNOT
   count. Use Qiskit as an independent simulator and optional transpiler.
4. Require every proposed winner to supply the full gate list, numeric residuals,
   qubit-order convention, and a derivation of exact angles. Recheck from the
   gate list independently, then verify the matrix identity symbolically or
   with an analytic proof. Numerical success alone does not establish an
   exact circuit or a lower bound.

The baseline supports local simplifications, while the wider search uses the
full freedom in the degenerate eigenspaces. Searching only the paper's fixed
Fourier eigenbasis could miss a shorter diagonalizer.

The 14-CNOT target described above has now been met. For any future search
below 14 CNOTs, agents should run bounded, non-duplicating searches in parallel.
After every failed or inconclusive run, finish its ledger entry, state the
specific hypothesis tested and its limit, and propose a distinct next
hypothesis using the observed failure. Reassign an available agent to that
hypothesis; do not retry an unchanged method or infer a lower bound from a
restricted search. A numerical candidate becomes the incumbent only after
an exact gate specification and independent exact verification. Continue
this loop until a better exact circuit is found or the experiment reaches a
stated resource limit or genuine blocker.

## Experiment ledger and agent handoff

`experiments.sqlite3` is the authoritative local record of attempts. Use
`experiment_log.py` before and after every bounded search run. SQLite prevents
simultaneous agents from overwriting one another. An exact duplicate of the
target, method, configuration, parent, and search-code hash is rejected before
work starts. Include the random seed, optimizer, iteration or time limit, and
topology or topology-generation rule in the configuration. Changing one of
these makes a distinct run.

```bash
python3 experiment_log.py summary
python3 experiment_log.py list
python3 experiment_log.py show 1
python3 experiment_log.py audit
python3 experiment_log.py export > attempts.jsonl

python3 experiment_log.py reserve --method delete-one-cx --agent agent-a \
  --description 'Remove baseline CNOT 0 and optimize local gates' \
  --config '{"deleted_cx":0,"seed":1,"optimizer":"scipy-least-squares","max_nfev":5000}' \
  --code search.py --parent 1
# Use the returned run ID after the computation finishes:
python3 experiment_log.py finish 3 --status inconclusive \
  --candidate best_candidate.json --notes 'No exact diagonalizer found.'
```

The final command illustrates the format; use the actual returned ID and
outcome. Finish unsuccessful runs too, including the best residual or the
reason for failure. Each submitted candidate is copied to `candidates/` under
its SHA-256 content hash. The ledger stores the ordered CNOT-topology hash,
candidate hash, CNOT count, residuals, validity flags, method, code hash,
agent, parent, and UTC timestamps. Search code is saved under
`code_snapshots/` by hash. `audit` rehashes the code and candidates and
rechecks numerical matrix properties. Identical topology hashes identify the
same ordered CNOT pattern;
different hashes do not prove mathematically different circuits.
The optional `evidence_level` field on a run is the submitting agent's claim;
`audit` verifies hashes and numerical matrix properties, not an analytic proof.

Runs 1 and 2 record the 18-CNOT paper baseline and Qiskit's level-3
transpilation. Both are numerical checks. A local Trackinizer issue (#1)
tracks the overall objective; the SQLite ledger and saved candidates remain
the source of truth for scientific attempts. Agent-session capture is off.
