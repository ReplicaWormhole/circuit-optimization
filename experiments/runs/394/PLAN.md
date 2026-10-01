# Run394: full98 broad batch32

Board80. New seeds9261600..9261631. All98 free, literal chart/schedule in config. Local/interior90 normal sigma [.3,1,2,3] cyclic by index; F8 independentuniform[-pi/2,pi/2]. One L-BFGS-B fit perstart max800iterations/maxfun20000/maxls30/ftol1e-16/gtol1e-12. Internal285s/external300s total; one thread. Atmost32starts. Earlystop normalizedFrobeniusloss<1e-16 is a promisingnearhit, NOT exactsuccess; independently check beforepolishing/certification. No hidden extra restart. Codehash c8e57295372068ad9ac01535c337d831c49bda0fdc64fcaf5ea31e1285cb7bcd; configfilehashcc6334176092ff9fd128aa1bb7710a064d297f74f07eb8ef8007204ffb8e3481; source/dependencyhashes inledger. Reserve precedes independentpreflight, which mustprecede optimizer.

Command: timeout --signal=TERM --kill-after=5s 300s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 experiments/runs/394/batch.py > experiments/runs/394/stdout.jsonl 2> experiments/runs/394/stderr.txt

Save perseed initialization, bestpoint, counts, caps, completeendpointlists and independentchecks; result.json checkpoint aftereachfit. Close even ifunsuccessful. Full12CXcompiledspec targetUV4=DU arbitraryeigenbasis; zerooptimizerloss needs independentgatevalidation and exactcertificate before incumbentchange. No failure exclusions. Goalcontinueuntilsuccess remainsactiveafter a failedbatch.
