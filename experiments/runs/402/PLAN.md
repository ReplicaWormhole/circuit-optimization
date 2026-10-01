# Reordered Bell-decoder full98 fit preparation

Board 106. This is the corrected preparation snapshot replacing the config
that was reserved as run 401; run 401 must remain preserved and unexecuted.
Structural motivation and exact prefix formulas are in
`collaboration/work/algebraic98bellprefix/REPORT.md` (hash pinned in config).
The target word is

`L0 F02 L1 F13 L2 CX01 A1 CX21 B1 CX31 L3 CX12 L4 F01 L5 F23 L6`.

`family.py` implements that chronology in a 98-coordinate chart and serializes
the same word to native `u3`, `cx`, and `xx_yy` gate lists. The one start is
the report's Bell-decoder prefix plus deterministic Gaussian noise on every
coordinate: `default_rng(9261801).normal(0,0.1,size=98)`. All coordinates remain
free. The original normalized off-diagonal objective is minimized with
PyTorch autograd and SciPy L-BFGS-B, at most 1000 iterations/20000 optimizer
evaluations, one thread, and 110-second internal/120-second external limits.
There are no restarts. A separately counted final best-point objective
recomputation is at most one and does not increase the SciPy evaluation count.

`search.py` saves initial native/compiled lists before fitting and best
native/compiled lists and `result.json` after the optimizer returns or the
cooperative internal budget stops. The search wrapper derives the run ID from
its reserved numeric workspace directory. Code and dependency hashes are
frozen in `config.json`; the config's own digest is to be recorded externally
by the reservation/ledger. The script itself does not reserve or run itself.

This directory is preparation only. Root owns ledger reservation, workspace
linking, verifier preflight approval and execution. Static syntax checks are
not matrix or serializer validation. A failed fit is bounded numerical
evidence only; any candidate needs independent full gate-list validation and
an exact certificate before promotion.

Run402 reserved and linked board106 before computation. Independent preflight must pass before this one command.

```bash
timeout --signal=TERM --kill-after=5s 120s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/402/search.py > experiments/runs/402/stdout.txt 2> experiments/runs/402/stderr.txt
```
