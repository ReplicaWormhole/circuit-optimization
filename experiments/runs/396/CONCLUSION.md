# Run396 conclusion

All eight frozen fresh starts completed in42.767seconds within the one-thread
110-second internal/120-second external batch cap. Seed9261707 reached the
1000-iteration limit; the others returned earlier. No nearhit occurred.
Every seed has its initial/best98 vectors, selectedassignments, objective,
gradientnorm, optimizercounts/termination, fullnative/compiled lists and hashes.
Frozen sources/inputs/limits/environment are in PLAN.md/config.json/result.json.

Best seed9261703: labeledloss0.1464466094067266,
originaloffdiagloss0.1357233047033634, maxoffdiag0.35355339341816694.
Best candidate: best_compiled.json (12 chronological CNOTs).
Independent labeledloss0.1464466094067264 and original0.13572330470336316.
The endpoint returns to the same loss plateau as run395 within roundoff.
Ledger396 is terminalfailed. No exactcertificate/incumbentchange.

Independent audit validates all eight RNG initializations, all16endpointlists
and globalbest lists, full98 matrix/native/compiler equivalence, unitarity,
counts, hashes, both residuals, and spectralmultiplicity in assignmenthistory.
See POSTFIT_AUDIT.json and VERIFIER_REPORT.md. All endpoints faildiagonalization.
Distinct-root labels changed during fitting, so the newobjective did explore
different assignmentbranches, without finding a solution in theseeightstarts.

Limitations: onlyeight frozen initializations; piecewisesmooth assignment can
affect quasi-Newton curvature. Optimizercalls exclude the one postfit scalar
original-loss diagnostic per retainedseed and gate-list checks; these were
declared beforeexecution and within the sametimecap. The tinydeadlinewindow
described in PLAN.md did notoccur. No failure implies localminimality, family
exclusion or global lowerbound.

Next hypothesis: one separatelyfrozen negativecurvature diagnostic at the
repeated lowplateau point, rather than more unchanged randomstarts. The earlier
run389 Hessian concerned a different point/objective. Any escape fit would
require a separatelyreserved budget and independentlyvalidated direction.
