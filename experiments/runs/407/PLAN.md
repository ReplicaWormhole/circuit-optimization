# Independent exact rational-ring verification

Board128, numerical_sol. The coordinator owns reservation, execution, ledger
finish and commit. No circuit evaluation is permitted during preparation.

The frozen source is verifier_sol/candidate12/candidate.json. This checker
imports no repository checker, NumPy or SymPy. Its arithmetic consists of eight
Fraction coefficients in Q[z]/(z^8+1), z=exp(i*pi/8), with reduction z^8=-1
and exact conjugation z->z^-1. All frozen one-qubit angles must be rational-pi
strings with phases and half-angles represented by integer powers of z.
Unsupported angles are rejected. Candidate phases and global phases are kept.

It constructs chronological U by independent exact one-qubit row updates and
CX row permutations, with qubit0 most significant. Independently building
V4|x0x1x2x3>=|x3x0x1x2>, it checks all256 entries of UV4=DU, discovering
one unique root in(1,i,-1,-i) for each row. It checks all256 entries of UU†=I,
exactly12CX and root multiplicities(6,3,4,3). The full exact matrix is retained.
Passing is an exact certificate for this literal gate list; optimality, total
spin and strong Schur claims are outside scope.

Frozen resources: seed0, no randomness, one thread, one deterministic90-second
invocation, no optimizer. The root must reserve the exact config.json object
with method "exact Fraction polynomial arithmetic in Q[z]/(z^8+1)", then link
board128. Matching active code/config/method/board ownership are checked before
circuit work. An exclusive execution marker prevents silent retries. Runtime
Python, Fraction source and math-extension hashes are checked. All config keys
are audited for membership and fixed values or content hashes.

Root copies certify.py,candidate.json,config.json,PLAN.md into the run workspace
created with experiment_log.py workspace. Authorized command:
`timeout 90s python3 experiments/runs/ID/certify.py`.
Internal90-second timer includes preflight and certificate work from main entry.
Finish failed/incomplete runs and preserve failure.json/execution marker too.
