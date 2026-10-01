# Run 400 FD freeze: source verification

This is a prepared, not executed, ten-matrix finite-difference check. The checker is generic for a numeric reserved run directory and refuses to run unless the frozen `run_id` equals that directory name. The root should copy `fd_check.py`, `fd_manifest.json`, and `fd_config.json` unchanged into the reserved run400 workspace and authorize one execution only after the board link exists.

The run399 source-audit mismatch was a metadata-only scope-string edit. I reconstructed the original JSON by changing only that scope field in a separate file, `experiments/runs/399/source_point_audit_expected.json`; its SHA-256 exactly matches the historical expected `9ac0678de91ca3eeed3ab799a1a1267ed96cf498302259e4694f442338c9e7f4`. The then-current `collaboration/work/verifier98curvature/PREFIT_POINT_AUDIT.json` was not overwritten. An immutable byte copy is bound in this freeze as `source_point_audit.json`, with the same exact hash.

Before freezing, all nine source entries in `fd_manifest.json` were independently re-hashed; every digest matched. The frozen source result is run398 `result.json`, SHA-256 `8ff1a48ab94146837f3756eb86e36bad87fc9dd2a0bbc78748371eddb651d1c2` (as recorded by its existing audit/manifest). The manifest binds the 98-coordinate point, the two full 98-component direction columns and hashes, two matching objectives, fixed labels copied from the run398 base assignment, and the steps `0.01` and `0.001`.

The source gate lists and matrix/evaluation plan are unchanged from the prior frozen design: two separate center calls and eight signed perturbation calls, for exactly ten NumPy coordinate-matrix evaluations. Perturbations use only the matching objective and retain the run398 fixed labels; no assignment is recomputed. The checker uses one thread, has a 30-second external cap in the plan, performs no optimization, and does not evaluate AD gradients or Hessians. It compares finite differences to the already saved run398 directional gradient/eigenvalue records.

Static validation only: Python AST parse passed; config-to-script and config-to-manifest hashes matched; run ID is 400; manifest point length is 98; both direction vectors have length 98 and their frozen f64le hashes match the manifest; four steps/objective-direction pairs imply the recorded ten calls. No circuit matrix or objective was evaluated while preparing this freeze.

Frozen SHA-256 values:

- `fd_check.py`: `d0c1f42c51cd46fc00f3671bf43daddc902dc57e867967acd2410369ea968706`
- `fd_manifest.json`: `7b395a9e7b92e15b8a4a704d537b05fb78b896f895eb7a691ab17a951475aea1`
- `fd_config.json`: `e09ee0e7480f86a4e3f298e2a2641dea3d3bcf812409698f0edaa2aac00addbe`
- `source_point_audit.json`: `9ac0678de91ca3eeed3ab799a1a1267ed96cf498302259e4694f442338c9e7f4`
