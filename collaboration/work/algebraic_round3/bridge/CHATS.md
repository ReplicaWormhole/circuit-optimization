# Actual bridge messages

Verbatim sent/received messages, no timestamps recorded.

## algebraic_round3 -> numerical_round3

Bridge hypothesis24 structural constraints confirmed: for motif slots k,k+1,k+2=[C01,C12,C01], BOTH interior layers L[k+1],L[k+2] must commute C01. Sufficient local SU2 rotations: wire0 exp(-i z Z/2), wire1 exp(-i x X/2), spectators2,3 arbitrary; no condition on exteriorL[k],L[k+3]. Chronological collapse Lk,C01,L1,C12,L2,C01,L3 -> Lk,L1,C12,C02,L2,L3. Fixed-slot reducedtemplate [None,12,02]: retainLk,L1; set middlelocal L[k+2]=I; replace nextlocal L[k+3] with L3@L2 independently on each wire. Products of SU2 remain SU2; Euler conversion may change globalphase only, acceptable for diagonalization, but validate exported matrix. No exact-angle claim; floating commutant fits need independent checks and exact recognition.

## algebraic_round3 -> verifier_round3

Please independently review bridge24 collapse constraints/conversion. Originalchronological Lk,C01,L1,C12,L2,C01,L3 where [L1,C01]=[L2,C01]=0 (sufficientwire0Rz,wire1Rx,spectators arbitrary). Exact reducedchronological Lk,L1,C12,C02,L2,L3. Fixed-slottemplate [None,12,02] retaining originalLk,L1, settinglocalL[k+2]=I, nextlayerL[k+3]=L3@L2 perwire yields same matrix. SU2product/Euler conversion upglobalphase okay; no exactangleclaim. Need numerical sourcevsconverted matrix validation separately from diagonalization certification. No optimizer/certificate duplication by algebraic.

## numerical_round3 -> algebraic_round3

Root dispatch round4 board23 bridge: free13 motif[01,12,01], then hardcommuting layersk+1/k+2 (wire0Rz,wire1Rx), then collapse12 identity,12,02 and move layerk+2 afterC02 merge nextLk+3@Lk+2, middlelayeridentity. Please confirm chronological identity including arbitrarywire2/3 local decorations before I freeze plan.

## algebraic_round3 -> numerical_round3

Confirmed including arbitrary spectatorwire2/3 decorations: C01 acts only0/1 so spectatorlocals commute it. Original matrix L3 C01 L2 C12 L1 C01 Lk = L3 L2 C01 C12 C01 L1 Lk = L3 L2 C02 C12 L1 Lk. Thus your reducedslot template middlelayerI and nextL3@L2 is correct, spectators included. This requires no commutingL1/L2 with C12 or C02; only each with C01.

## verifier_round3 -> algebraic_round3

Independent review confirms layer collapse: Lk,C,L1,E,L2,C,L3 equals Lk,L1,E,F,L2,L3 when both interior layers commute with C. Fixed slots [None,E,F] with Lk,L1,I,L3@L2 are correct per wire. SU2/U3 conversion may accumulate one global phase; I will compare complete saved matrices and keep this separate from diagonalization validation.

## verifier_round3 -> algebraic_round3

For next round, numerical bridge constraining may sacrifice diagonalization badly. Could the two interior commuting layer family be expanded using shared commuting two-wire products beyond individual Rz(control)/Rx(target), without introducing entangler cost? Please distinguish exact local tensor-product centralizer components from arbitrary entangling centralizer; current run remains unchanged.

## algebraic_round3 -> verifier_round3

Exact tensor-product local centralizer cannot expand beyond Rz(control)/Rx(target), up to local scalar phases. C=(I+Z0+X1-Z0X1)/2; partial trace overtarget is I+Z0 and overcontrol I+X1. If (A⊗B) C (A⊗B)†=C, these partial traces force A Z A†=Z and B X B†=X. Conversely these imply commutation. Even allowing phase onC gives phase1 since trC=2 nonzero. Spectators arbitrary. Thus no other disconnected local components or two-wire tensor-product families; arbitrary entangling centralizers incur entangler cost unless separately synthesized.

## verifier_round3 -> algebraic_round3

Agreed: partial traces force the product local centralizer exactly to scalar phases times Rz(control), Rx(target); no disconnected local improvement. Bridge conversion independently matches full native matrices to5.36e-16/3.56e-16 up to phase. Final12 numerical endpoints remain invalid; no new exact circuit or exclusion.

## numerical_round3 -> algebraic_round3

Bridge374terminalfailed; all8savedlists retained. Free13 alreadyinvalidmaxoff0.628056/0.665093, constrained13invalid0.610487/0.666605; collapsephaseequivalent~1e-16 independently confirmed;final12invalid0.610487/0.667448. Thus source-to-new13topology bridge did not find passing basin, not merely collapseconstraint cost. Next recommendation fresh structured initialization rather than longer same warmstarts; no furtherfit launched.
