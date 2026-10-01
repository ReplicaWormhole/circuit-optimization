# Independent bridge verification

Board25, no optimizer or duplicate exact certificate. All artifacts are intentional retained evidence.

## Structural identity

Chronological Lk,C01,L1,C12,L2,C01,L3 collapses to Lk,L1,C12,C02,L2,L3 if both interior layers commute with C01. Product-layer layout [None,12,02] with local layers Lk,L1,I,L3@L2 is correct. The sufficient control-Rz/target-Rx restriction is also the complete tensor-product local centralizer, up to scalar phases: partial traces of CNOT = (I+Z0+X1-Z0X1)/2 force A Z A†=Z and B X B†=X under (A tensor B) conjugation. Its nonzero trace excludes nontrivial projective centralizing phase. Spectator locals remain arbitrary. Entangling centralizers would require their own gate cost.

## Independent native matrices

`compare.py` constructs complete native U3/CX matrices independently using explicit bit columns and local Kronecker products. Commands compare `case0_constrained13.json` with `case0_collapsed12.json` and corresponding case1 in numerical_round3/bridge. Saved output `comparison0.json`, `comparison1.json` binds source hashes. Both assertions pass, counts13->12. Maximum differences up to global phase: 5.36049202688126e-16 and 3.55715613020066e-16. This numerical comparison validates conversion to floating precision; the symbolic commutation identity justifies the structural collapse. It does not certify diagonalization.

## Final endpoints

`python3 check_circuit.py collaboration/work/numerical_round3/bridge/case0_final12.json` and corresponding case1 both return expected exit1 INVALID, each12CX. Maximum offdiagonal errors:0.6104871878877347 and0.6674475006279953; unitarity errors6.661338147750939e-16 and1.7763568394002505e-15. No exact certification is warranted and no topology exclusion follows. Accepted13CX remains unchanged; optimum remains unproved.

## Adversarial review

Checked collapse order, both commutations, local multiplication order, global-phase convention, native counts, source binding, and invalid endpoint scope. No confirmed issue in bridge artifacts. Root found the prior readable circuit omitted interpretation of numeric eigenlabels; exporter and regenerated artifact now explicitly define D_rr=i**labels[r]. That confirmed presentation defect is corrected without changing any gate. Remaining risk: failed optimization is bounded evidence only and no global optimum proof exists. Coordinator owns shared files and commit.

All eight saved gate lists also evaluated through shared native checker; `all_saved_checks.json` retains results. All are invalid. Counts are13 for free/constrained bridges and12 for collapsed/final cases; constrained and collapsed offdiagonal residuals agree to floating precision. The centralizer proof describes the exact local tensor-product class; failed100-iteration fits do not exclude this class or any topology, and the hard coordinate mask is the particular parameterization optimized in run374.
