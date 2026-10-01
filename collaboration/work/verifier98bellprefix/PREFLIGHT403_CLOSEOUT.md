# Run403 preflight closeout

The one authorized PLAN command exited 1 after 1.004 s. It stopped at the pinned-runtime version guard with `pinned NumPy/SciPy/Torch runtime differs`. The guard runs before base construction, RNG generation, and the matrix loop; therefore the run performed zero coordinate-point evaluations and zero matrix builds. No `preflight_result.json` was created.

The mismatch is a version-string suffix: frozen/runtime metadata reports NumPy 2.3.5, SciPy 1.16.3, Torch 2.11.0, while imported module versions are NumPy 2.3.5, SciPy 1.16.3, Torch 2.11.0+cu130. The preflight config was generated from package metadata, while the checker compares the exact `torch.__version__` string. This is a guard/configuration mismatch; it provides no circuit, initializer, or serializer validation.

No source inputs were edited and the failed command was not retried. Full command, hashes, stdout/stderr, versions, and zero-call accounting are in `PREFLIGHT403_EXECUTION_AUDIT.json`. Root owns ledger/board closeout and any separately reserved corrected preflight.
