# Frozen-check plan

Board hypothesis: 116. The root coordinator must reserve the bounded
numerical and exact check before execution, then copy `check_baseline.py` and
`config.json` into the corresponding `experiments/runs/<ID>/` workspace.

The one authorized command is:

```bash
python3 experiments/runs/<ID>/check_baseline.py
```

The script refuses to run from this preparation directory. In a reserved
workspace it verifies the ledger status and code hash, pinned checker hashes,
Python/NumPy/SymPy versions, one-thread setting, exact gate vocabulary, and
the 19-CX count before writing `candidate.json`. It runs the NumPy checker
once with a 30-second timeout and the exact cyclotomic checker once with a
60-second timeout. It preserves checker stdout/stderr and writes `result.json`
even when either checker times out. There are no optimizer, derivative,
search, retry, or restart paths.

The planned input is an explicit rational-pi 19-CX comparator, not a search
start and not a candidate for 12-CX promotion. Its complete ordered gate
builder is frozen in `check_baseline.py`; the root copies the exact config and
records the script/config hashes in the ledger workspace.
