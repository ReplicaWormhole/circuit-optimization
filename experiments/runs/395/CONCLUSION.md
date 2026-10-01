# Run395 conclusion

The one frozen spectrum-assignment fit from run394 seed9261620 returned after
one iteration and three objective-gradient evaluations in0.06673seconds.
Termination: relative reduction tolerance. Initial gradient norm1.34347e-8;
best gradient norm9.93068e-9. This is not a local-minimum proof.

Best labeledloss0.1464466094067263; originaloffdiagloss0.1357233047033631;
maximumoffdiag0.3535533925341737. No meaningful escape from the startingpoint.
Best compiled12-CNOT list is best_compiled.json; ledger395 is terminalfailed.
All98 coordinates were free, with no random draws or restarts. Complete frozen
inputs, code/dependency hashes, optimizerlimits, selectedassignment history,
and environmentversions are recorded in PLAN.md/config.json/result.json.

The independent audit checked allfive endpoint vectors and their ten full
native/compiled gate lists, counts, phaseequivalence, unitarity, hashes, both
objectives and correctspectrumassignment. See POSTFIT_AUDIT.json and the
verifier's report. No endpointdiagonalizes V4; no exactcertificate or promotion.
The incumbent13-CNOT and exact14baseline remain unchanged.

The ledger's finish-time note records that the postfit audit was thenpending;
the completedaudit is nowretained in thisworkspace without rewriting the row.
The failure neither excludes thisfamily nor establishes a lowerbound.
Next distinct hypothesis: use the validatedobjective from freshinitializations,
rather than only from an already stationarypoint. Separatelyreserved run396
tests eight suchstarts with a new finitebudget.
