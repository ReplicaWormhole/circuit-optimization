# Reordered full98 Bell-prefix static review

Board hypothesis 103. This review is source/config-level only; no coordinates were generated and no matrix, objective, derivative, fit, or search was evaluated here.

## Static findings

The declared chronological word is `L0 F02 L1 F13 L2 CX01 A1 CX21 B1 CX31 L3 CX12 L4 F01 L5 F23 L6`. The matrix evaluator and serializer in the reviewed family enumerate the same chronology; the four directed CX are `(0,1),(2,1),(3,1),(1,2)`, and native XX/YY pairs are `(0,2),(1,3),(0,1),(2,3)`. The compiler emits two CX per native XX/YY gate, giving four native CX plus eight compiled CX, or twelve compiled CX total.

The chart covers 98 entries with slices `local=z[:84].reshape(7,4,3)`, `ab=z[84:90].reshape(2,3)`, and `f=z[90:].reshape(4,2)`. Source loops visit all seven layers and all four wires; both three-coordinate A1/B1 rotations on q1 are inserted between their designated CX gates; all four two-coordinate F gates are used. The optimizer passes the full 98-vector with no fixed-coordinate mask. There are seven 4-wire local layers (28 U3s) and the two A1/B1 U3s, totaling 30 one-qubit gates in the serialized family.

The Bell decoder initializer places its two input locals in L0, negative-angle F02/F13 in the first two interaction slots, and X0/X1 in L1/L2. The algebraic review independently checked the SO(3) chart vector for `Ry(pi/2) S†`, the negative interaction signs, and the nonduplicating chronological schedule. The serialized base prefix consists of the first fourteen gates (L0, F02, L1, F13, L2); because the later CX suffix remains present in the full family, only this prefix should be compared against direct `E_pair†`.

The initial perturbation is `default_rng(9261801).normal(0,0.1,size=98)` added to the Bell base, with one start and no restarts. Fit limits are 1000 iterations, 20,000 optimizer objective/gradient evaluations, one thread, 110 seconds internal and 120 seconds external. The runner documents at most one final best-point objective recomputation outside SciPy's maxfun count. It writes the initial native and compiled gate lists before fitting; the run's best lists and result are saved after the optimizer or cooperative internal stop. The checker uses a one-shot execution marker.

## Freeze issue and corrected snapshot

The first reserved source snapshot 401 was not executable: `search.py` read `config["noise_sigma"]`, but 401 `config.json` only contained `initialization.noise_sigma`; because the runner writes its one-shot marker before generating the initial point, this would fail with a missing-key error and block retry. Root closed 401 without execution and preserved it. The corrected, distinct source snapshot is `collaboration/work/numerical98bellprefixv2/`; it adds top-level `noise_sigma: 0.1`. Its config points to matching search/family hashes, and all five declared dependency hashes match at review time. Both Python sources pass AST parsing. The native-aware checker fix is present: native lists use `native_gate_check`, compiled lists use `check_circuit`.

The independent bounded preflight checker is being prepared separately. It is limited to three parameter points: exact Bell base, the frozen seeded start, and a deterministic all-nonzero probe. It will report six full-family matrix builds (Torch and independent NumPy per point), six serialized-list matrix builds (native and compiled per point), and two extra products for the serialized Bell prefix/direct `E_pair†` comparison. No preflight has yet been executed; it must bind the actual reserved corrected fit workspace and its hashes before root authorizes it.
