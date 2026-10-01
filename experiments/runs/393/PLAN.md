# Run 393: one fresh full98 fit

Board75; structural73; independent verifier76; coordinator77. Seed9261501, exactly numpy.default_rng(seed).normal(0,.3,size=98), full literal vector in config.json. Locals84+A1/B1 six+F eight; ALL98 free. Native chronology L0 CX01 A1 CX21 B1 CX31 L1 F02 L2 F13 L3 CX12 L4 F01 L5 F23 L6. Four CX and four F=exp(i aXX+i bYY); compiled12CX. Qubit0MSB, gateschronological; loss=||offdiag(UV4Udagger)||F^2/16, unrestricted eigenvalue order/eigenspace bases.

ONE L-BFGS-B start with Torch AD, maxiter300/maxfun10000/maxls30/ftol1e-15/gtol1e-10; internal110second deadline savesbest, external120second timeout, onecompute thread. ZERO restarts/followups. Code SHA256 b19f1aaae55d705fb2487ecb04ca1f8e878688404c2e9a03133edb4e73b305aa; config file SHA256 d347c2aeff25828c476354e10baeb0df184022eb8c1ecfc859e7bcc2e7b89b52; complete dependency hashes and frozenvector in config.json and ledger. No input is an inherited385 checkpoint.

Independent preflight must approve exactcode/config and RNG/initialmatrix/serializer/compiler BEFORE executing search. Execute from repository root:
`timeout --signal=TERM --kill-after=5s 120s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/393/search.py > experiments/runs/393/stdout.jsonl 2> experiments/runs/393/stderr.txt`

Save result/history/bestparams and initial/best native/compiled fullgatelists; independent NumPy/SciPy checks of allfourlists afterrun. Bestcandidate selected by lowest objective among evaluatedpoints, including linesearchsamples. Numerical success is NOT an exact certificate; no incumbent change absent independent exact certification. Failure is one boundedstart, not family/topology exclusion. Never rerun this workspace; reproducibility needs a fresh separately reserved run.
