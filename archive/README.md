# Historical computation archive

`legacy_root/` preserves the former flat root namespace for historical search
scripts and saved outputs. The 1,398 moved files have unchanged bytes; the
archive manifest records their SHA-256 hashes. Shared helper code is copied
into the archive at its migration version, and links to the active ledger,
snapshots, collaboration records, and accepted baseline preserve old relative
paths.

From the repository root, verify the move with:

```bash
python3 archive/verify.py
PYTHONDONTWRITEBYTECODE=1 python3 archive/legacy_root/search13_opening02_recovery_verify.py
```

The second command independently rechecks the 66 saved candidate gate lists
for historical run 354. It writes its validation JSON back to the archive.
For other old commands that named a root script or JSON file, prefix that file
with `archive/legacy_root/`. Reserve and bound any **new** computation in the
live ledger and collaboration board before running it.

The accepted native baseline records `topology14_exact_matchgate_candidate.py`
as its generator. A root-level link keeps that recorded path resolvable while
the generator and its companion inputs live in the archive. Qiskit is required
to run that historical generator; it is optional for the active exact checks.

The archive is a research record. The active commands, incumbent and current
status remain at the repository root and in `collaboration/` and `docs/`.
