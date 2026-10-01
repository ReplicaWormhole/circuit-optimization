# Run368 exact projective attempt: timeout, not certificate

Command: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 collaboration/work/algebraic13_new/canonical_projective_verify.py`.

One bounded exact symbolic run, one CPU thread, hard wall timer300 seconds. Input: `canonical_algebraic_ansatz.json`, parent run366. The program exactly verified all60 local `sin²(theta)+cos²(theta)=1` identities and counted13 CNOTs. The independent coordinator review `collaboration/work/root/radical13_local_exact_review.json` establishes exact realness and nonzero projective representatives. An elementary written proof is `canonical_local_unitarity_proof.md`.

The projective representation `(1+cos(theta))I−i sin(theta)P` avoids half-angle radicals and preserves the homogeneous target equation; `cos=-1` uses pureP. Exact expression expansion nevertheless accumulated large scalar factors. Completed milestones: gate35 at86.6s (2614 operation nodes maximum entry), gate40 at177.2s, gate45 at251.1s (1349 maximum nodes). The300-second alarm stopped assembly before the complete matrix or any global UV4=DU row checks. `canonical_projective_progress.json` records `stage:timeout`, `exact_diagonalizer:null` and the limitation.

This is partial exact local evidence, not a thirteen-CNOT global certificate or disproof. Run368 and hypothesis12 are closed inconclusive. No automatic retry occurred. A distinct next algorithm could remove per-gate global factors via normalized projective `I−i tan(theta/2)P`, simplify the tangent algebra before assembly, and exploit repeated-angle/radical relations. That algorithm was not run here.

Adversarial review: no floating-point zero tests in the symbolic verifier; timeout records incomplete work honestly; source scripts and initial guesses preserved; no shared checker/native-gate code or status documents edited. Local unitarity proof and coordinator independent review cover the nonzero factor condition. Remaining obligation is every exact global residual entry. All artifacts retained intentionally; coordinator owns commits.
