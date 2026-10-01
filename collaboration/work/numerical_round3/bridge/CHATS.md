# Actual bridge exchanges

## numerical_round3 → algebraic_round3

Root dispatch round4 board23 bridge: free13 motif[01,12,01], then hardcommuting layersk+1/k+2 (wire0Rz,wire1Rx), then collapse12 identity,12,02 and move layerk+2 afterC02 merge nextLk+3@Lk+2, middlelayeridentity. Please confirm chronological identity including arbitrarywire2/3 local decorations before I freeze plan.

## algebraic_round3 → numerical_round3

Bridge hypothesis24 structural constraints confirmed: for motif slots k,k+1,k+2=[C01,C12,C01], BOTH interior layers L[k+1],L[k+2] must commute C01. Sufficient local SU2 rotations: wire0 exp(-i z Z/2), wire1 exp(-i x X/2), spectators2,3 arbitrary; no condition on exteriorL[k],L[k+3]. Chronological collapse Lk,C01,L1,C12,L2,C01,L3 -> Lk,L1,C12,C02,L2,L3. Fixed-slot reducedtemplate [None,12,02]: retainLk,L1; set middlelocal L[k+2]=I; replace nextlocal L[k+3] with L3@L2 independently on each wire. Products of SU2 remain SU2; Euler conversion may change globalphase only, acceptable for diagonalization, but validate exported matrix. No exact-angle claim; floating commutant fits need independent checks and exact recognition.

## algebraic_round3 → numerical_round3

Confirmed including arbitrary spectatorwire2/3 decorations: C01 acts only0/1 so spectatorlocals commute it. Original matrix L3 C01 L2 C12 L1 C01 Lk = L3 L2 C01 C12 C01 L1 Lk = L3 L2 C02 C12 L1 Lk. Thus your reducedslot template middlelayerI and nextL3@L2 is correct, spectators included. This requires no commutingL1/L2 with C12 or C02; only each with C01.

## verifier_round3 → numerical_round3

Claimed board 25. Please retain native bridge13 and collapsed12-before-polish JSON for each motif, in addition to final12 endpoints. I will independently compare their complete matrices up to global phase and check final native gate lists. Ensure both interior layers constrain wire 0 to Rz, wire 1 to Rx.

## numerical_round3 → verifier_round3

374reserved/link23 launched session43826. bridge/run.py saves case0/1_free13,constrained13,collapsed12,final12.json. Both interior layers hardmask0[x,y],1[y,z] at k+1/k+2; collapsed middleI nextL3@L2. Matrixequivalence assertion1e-10 upglobalphase before polish. Will notify completion.

## numerical_round3 → verifier_round3

Bridge374 done10.649s, all8savedgateJSONs ready in collaboration/work/numerical_round3/bridge/. Native13→collapsed12 equivalence residual5.836925621952896e-16/3.3306690738754696e-16. Final12invalid maxoff0.6104871878877347/0.6674475006279953. Please independently check allsavedgates and phaseequivalence; no exact certification warranted.

## verifier_round3 → numerical_round3

Independent bridge comparison passed on complete native matrices: case0 difference5.36049202688126e-16, case1 difference3.55715613020066e-16 up to global phase; counts13->12. Final12 checker results INVALID maxoff0.6104871878877347/0.6674475006279953. No exact certificate needed. Comparisons saved in verifier_round3/bridge.

## numerical_round3 → algebraic_round3

Bridge374terminalfailed; all8savedlists retained. Free13 alreadyinvalidmaxoff0.628056/0.665093, constrained13invalid0.610487/0.666605; collapsephaseequivalent~1e-16 independently confirmed;final12invalid0.610487/0.667448. Thus source-to-new13topology bridge did not find passing basin, not merely collapseconstraint cost. Next recommendation fresh structured initialization rather than longer same warmstarts; no furtherfit launched.
