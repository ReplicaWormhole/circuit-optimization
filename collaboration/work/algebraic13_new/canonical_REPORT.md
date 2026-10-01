# Canonical13 endpoint and radical recognition

Run366 / hypothesis11 uses source run365's stripped canonical physical rotation representation (60 axis rotations and13 CNOTs). Original source and complete coordinate/gate/wire mapping plus source SHA256 are preserved in `canonical_source_map.json`.

Command: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 collaboration/work/algebraic13_new/canonical_endpoint.py`.

Budget: at most8 constrained fits,400 evaluations each, one thread. Used2 fits: recognized constants (39 fixed,1 evaluation), then initial RY q1=0 and RY q2=π together (41 fixed,5 evaluations). The remaining-coordinate Jacobian had no usable nullspace direction, so no additional fits ran. This paired endpoint succeeds numerically: maxoff5.652e-16, eigenvalue-root error1.776e-15. The endpoint is a tested hypothesis, not an assumption inherited from the original source.

All19 remaining angles recognize rational multiples of π or clean radicals. Five previously complicated angles have sin² values:

| Coordinate | sin²(theta) |
|---|---|
|27|(252−45√2)/367|
|38|(3−√2)/7|
|39|5/12−5√2/48|
|44|(184+72√2)/367|
|45|(19+3√2)/28|

Other sin² values are2/3 (14),3/5 (15,48),1/6 (21),1/3 (42,54),2/5 (58),5/6 (59). Additional rational angles are7π/8 (12),π/2 (13,19),−5π/8 (18),0 (25,31).

`canonical_algebraic_ansatz.json` supplies exact signed sin/cos and branch-correct half-angle entries for every rotation, with rational-π metadata where available. The radical prescription matches the numerical endpoint trig values within3.33e-16. This recognition alone is not an exact UV4=DU certificate. Coordinator independently checked it at100 decimal digits in run367; that remains numerical evidence. Separate run368 / hypothesis12 performs bounded exact projective verification.

Full π/24 rounding remains invalid (maxoff0.09474), demonstrating why the algebraic radicals must be retained.

Adversarial review: source reset, actual paired wire/coordinate mapping, fixed constants, rotation signs and 13-CNOT count inspected. No original work, checker code, native-gate code or shared status documents overwritten. Python syntax compilation passes. All assigned artifacts intentionally retained. Run366 completed with numerical evidence; no incumbent promotion or local commit by collaborator. Exact verification outcome must be taken from run368, not this recognition report.
