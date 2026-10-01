# Run404 preflight closeout

The single frozen PLAN command exited 0 in 1.729 s. It evaluated exactly three coordinate points: the Bell-decoder base, the seed-9261801 / sigma-0.1 initial point, and a deterministic all-nonzero probe. Accounting is 6 full-family matrix builds, 6 native/compiled serialized-list builds, and 2 additional prefix/direct-decoder products, for 14 top-level matrix builds total. The maximum reported matrix error is 1.44e-15; the serialized 14-gate Bell prefix agrees with direct `E_pair^dagger` up to phase with error 1.87e-16.

For all three points, Torch and independent NumPy full-family matrices agree up to phase, native serialization agrees with the independent chronological evaluator, and compiled lists agree with native lists. Each native circuit has 4 CX plus 4 XX/YY gates, and its compiled list has 12 CX. Source/config hashes and the seed construction passed the checker guards. Static slice review confirms all 98 coordinates remain in the family.

This validates the initialization and representation path at three fixed points only. It evaluates no target objective and ran no optimizer, derivative, or fit; it does not validate a diagonalizer or candidate. The full raw preflight result and command/hash record are saved in `preflight_result.json` and `PREFLIGHT404_EXECUTION_AUDIT.json`. Root owns ledger/board closeout and authorization of fit402.
