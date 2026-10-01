# Run402 bounded reordered full98 fit

Frozen seed9261801, exact Bell-decoder base plus default_rng Gaussian sigma0.1
on all98 coordinates; all98 remain free. Original offdiagonal loss; SciPy
L-BFGS-B with Torch AD gradient, max1000iterations/20000optimizer evaluations,
110sec internal/120sec external/one thread, no restart. Run404 independent
three-point preflight passed before the one fit command.

Command exited0:113iterations127objective-gradient calls plus one explicitly
counted final best-point loss recomputation,1.364887827sec. Bestloss is
0.3750000000000021, gradientnorm3.92257e-8. Full best native and compiled lists
are retained; compiledlist has12CNOTs, maxoff0.9999999999999969 and is invalid.
Source/config, complete initial/base/noise/best vectors, history, versions,
stdout/stderr and best candidate are retained. Ledger402 failed; board106
inconclusive. Independent endpoint audit passed: initial/best coordinates and all four full
lists agree upglobalphase below8e-16, unitarity below1.2e-15; independent
bestloss0.375000000000002 and maxoff0.9999999999999974 confirm invalidity.
See POSTFIT402_AUDIT.json. No exact certificate,
incumbent change, minimum proof or family exclusion. One small perturbation
cannot determine representability; next investigate constructive Bell-label
exchange suffix initialization before further unchanged starts.
