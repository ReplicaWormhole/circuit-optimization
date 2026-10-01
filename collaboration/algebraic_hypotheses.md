# Algebraic team's first collaboration round

Date: 2026-09-26. This round inspected existing derivations and search records; it did not run a new optimizer or re-prove the incumbent. Candidate counts below are analytic consequences of the documented block construction, pending the team's independent gate-list verification.

## Useful native-gate identity and handoff

Use four wires, most-significant wire 0, chronological gates, and right shift `V4 |x0x1x2x3> = |x3x0x1x2>`.

The exact14 circuit is a six-CNOT prefix followed by local layers and four blocks `F(a,b)=exp(i a XX+i b YY)`. The first two disjoint blocks act on 02 and 13 with `(a,b)=(pi/4,pi/8)`; the last two act on 01 and 23 with `(pi/8,pi/8)`.

An independently tunable XX/YY native family therefore gives 10 two-qubit gates. It can be uniform: `F(pi/4,0)` is locally CNOT-equivalent, so it also realizes each of the six prefix CNOTs with one native gate and free single-qubit changes of basis. This changes the elementary gate model; it does not improve the CNOT bound.

A single-parameter exchange family `E(t)=exp(i t(XX+YY))` must not silently substitute for unequal first-layer blocks. The exact commuting-Pauli identity is

`F(a,b)=E((a+b)/2) X_c E((a-b)/2) X_c`.

Conjugating by `X_c` changes the sign of YY and preserves XX, proving the formula. Chronologically apply `X_c`, `E((a-b)/2)`, `X_c`, `E((a+b)/2)`. Thus the unequal first blocks require this two-exchange construction with angles pi/16 and 3pi/16; each equal terminal block is one exchange at pi/8. The hybrid CX plus variable-exchange family gives a 12-entangler upper bound. This is an upper bound on the decomposition cost, not a proof that each block requires precisely that many exchange gates.

The last prefix CNOT is 3->2, so neither first pair 02 nor 13 shares its pair; there is no immediate same-pair boundary fusion. This observation and the native-family distinction were communicated directly to the native-gate researcher.

Evidence: `topology14_exact_matchgate_derivation.md`, first 13 gates of `topology16_15_exact_candidate.json`, and `topology16_15_derivation.md`.

## H1: Nine independently tunable XX/YY native gates

Hypothesis: permitting the unequal-angle native family throughout allows a nine-gate diagonalizer after jointly changing the prefix and the degenerate-sector basis.

A concrete bounded next experiment is eight starts, at most 300 optimizer iterations each, on nine-gate architectures obtained by deleting one native gate from the 10-gate construction, using disjoint deletion locations. Fit all one-qubit frames and both native angles jointly. First run exact baseline construction and verify native-versus-CNOT matrix equivalence, then perform the bounded deletions; do not freeze the prefix target or output rows.

The objective is off-diagonal norm of `U V4 U†`, with any eigenvalue order accepted. Record topology, two angles per native gate, seeds, per-start residual and runtime. A numerical success requires independent full-unitary simulation and unitarity verification. An exact winner additionally needs rational-angle reconstruction or another exact specification and an independently verified `U V4 = D U` identity. Native 9 versus CNOT 14 must remain distinct registry entries; compilation back to CNOT is reported separately.

If all eight fail, conclude only that these deletion architectures and starts failed. A useful handoff is the best residual and Jacobian null directions, allowing the numerical researcher to propose a changed pair schedule rather than rerun identical starts.

## H2: A genuinely nonmonomial five-CNOT prefix

Hypothesis: a five-CNOT prefix containing non-diagonal one-qubit gates can compensate jointly with the existing eight-CNOT tail while escaping the phase-polynomial obstruction.

The old prefix screens exhausted five-CNOT parity-Rz paths for fixed first-pair endpoint conditions. Run 279 strengthened this to arbitrary continuous diagonal input phases constant on shift orbits. This excludes the stated monomial architecture; it does not exclude arbitrary local frames and a jointly changed first pair. A previous nonmonomial trial (run 282) failed numerically and is a single restricted construction.

A distinct bounded test is four prefix schedules with interaction orders that are not the incumbent prefix with just one deletion, each with two starts and a cap of 300 iterations. Include an Rx or Ry freedom before and after each prefix CNOT, and fit all tail local layers jointly; the input and intermediate bases must remain free. Use the actual diagonalization loss, not process fidelity to the fixed incumbent matrix. Before running, the numerical agent should compare the schedules and seeds against the ledger to avoid duplicating prior experiments.

Required evidence is the complete 13-CNOT gate list plus independent diagonalization check; an exact claim requires exact gate specification and an independent exact identity. A near-hit should be handed to the algebraic agent with its final gate list and parameter sensitivity, not only a scalar loss. Failure is numerical evidence about eight starts, not a no-go theorem.

## H3: Joint degenerate-sector rotations across the middle boundary

Hypothesis: a rotation in a three-or-more-dimensional repeated-eigenvalue sector enables a changed intermediate basis and a seven-CNOT tail, beyond finite row swaps and isolated Givens tests.

The incumbent has eigenspace dimensions 6, 3, 4, 3. Row rotations confined to any one sector preserve diagonalization. Existing fixed-output row-permutation and small-Givens screens, and exact three-wire regrouping obstructions, have restricted scopes; neither proves that the full continuous product U(6) x U(3) x U(4) x U(3) cannot help. Existing continuous gauge searches also have to be distinguished by the specific target and boundary they optimized.

A bounded next test should choose one three-row block in the six-dimensional eigenvalue-1 sector, vary a full SU(3) rotation, and fit a seven-CNOT tail with the six-CNOT prefix fixed. Compare exactly two independently chosen interaction schedules surviving the existing mixed-cut screen, with two starts each and 300 iterations maximum. Selection requires checking the ledger against prior continuous gauge fits; do not repeat their gauge parameterization and schedule together. This explicitly changes the fixed target, whereas many old rank exclusions assume a fixed target.

Promising gauge angles and the full resynthesized tail go to algebraic reconstruction. Exact rank survival is only a necessary filter, and no number of unsuccessful optimizer starts proves a global lower bound.

## Scope of established obstructions

- `algebraic13_terminal_diagonal_gauge.md`: exact no-saving statement for an isolated terminal same-pair diagonal gauge. Cross-pair gauges and jointly changed earlier boundaries are outside its scope.
- `docs/history/search13/search13_firstpair_prefix_derivation.md` and `docs/history/SEARCH_14_STATUS.md` run 279: catalogued rational phase-polynomial models, then continuous orbit-constant diagonal phase models, with fixed first-pair endpoint conditions. Arbitrary nonmonomial prefixes are outside scope.
- `docs/history/SEARCH_14_STATUS.md` runs 298 and 311: fixed overlapping three-wire middle operators cannot be replaced by three CNOTs. A changed intermediate basis is outside scope.
- `docs/history/search13/search13_output_eigenrow_gauge_derivation.md`: finite swaps and single real Givens rotations with the stated boundaries; its failed bounded numerical fits are not exact exclusions of other gauges.
- README records the global ancilla-free all-to-all arbitrary-single-qubit CNOT interval `6 <= C_min <= 14`. This round did not independently reproduce its certificates, and it does not transfer automatically to native XX/YY or exchange families.
