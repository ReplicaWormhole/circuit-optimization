# Bridge round4: run374, hypothesis23

Two truebraid13-CX topologies replace source[01,32,12] by[01,12,01] at locations0/6. Canonical warmstart plus0.08 Gaussian noise seeds926401/926402. Stage1 all168 local axis coordinates free,100iterations. Stage2 each interiorlayer wire0[x,y]=0 andwire1[y,z]=0,100iterations. Collapse to[None,12,02], interiorsecondlayerI andnextlayerL3@L2. Stage3 all168 coordinates free,200iterations. Caps120seconds permotif240total,one compute thread. Reservation374 andlink23 preceded fits. Actual10.649seconds. Full source/dependency/code hashes frozen inplan.json. Eight gate lists retained.

|Motif|Free13maxoff|Constrained13maxoff|Final12maxoff|13to12matrixequivalence|
|---|---|---|---|---|
|0|0.6280560871346066|0.6104871880414621|0.6104871878877347|5.84e-16|
|6|0.6650929662143774|0.6666053704526018|0.6674475006279953|3.33e-16|

Independent verifier compared complete exported matrices up to globalphase, residual5.36e-16/3.56e-16, counts13->12; finalgate checks rejected both. No exact certification needed. Structural identity independently checked by algebraic24; numerical equivalence preserves invalid diagonalization, not a passing construction. Optimizer stopping does not exclude topologies. CurrentCX interval remains6<=Cmin<=13.

Reproduce `python3 collaboration/work/numerical_round3/bridge/run.py`; this reruns fits so reserve a NEW ledger attempt first. No furtherfit dispatched. Next useful hypothesis needs fresh structured initialization because unrestricted13 fits were already invalid, not merely constraint degradation.

Adversarial review: reviewed chronological multiplication, hardmask, SU2 conversion/globalphase, nativegate counts, reservationtiming andfitbounds. No confirmed defect; residual risk is finite basin exploration and numerical rather than exactangle evidence. Own generatedJSONs/chatnotes/scripts intentional evidence; coordinator owns commit, shared checker unchanged.
