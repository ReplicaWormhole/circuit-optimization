# Static review of frozen full-local deletion-refit pack

Board130; no circuit matrix, objective, gradient or optimizer executed during
review. Frozenfit f90819b38d17156c1f20bbf4b21274f4023a4f39829f52e524d8a7426cc5593c;
family b7361e44000a2a0221a61fdef1173c6d68f509995aa279c3b9446f9b3dd04e0b;
config277e8509a1eb57531fe800f79c646b1c1e6564fbc3c66c6641576cf3cb7f9748.

Source406candidate exacthash58a862a4acef4e768c4613e9b142504592f5237a68539941148653a7f083eea2.
The topology generator deletes each chronological CX ordinal0..11 once and
preserves remaining edges/order. Each word has11CX and12layers with4independent
SU2 rotationvectors perlayer; all144real coordinates remain free, including
control/spectator wires. The collapse multiplies chronological local gates on
the left and converts each local matrix up to a global determinant phase; the
product of discarded local phases is one global fullword phase. Serializer
outputs48U3+11CX=59gates with all four wires in each layer.

OnePCG64(seed10012101) generator preallocates12successive144coordinate Gaussian
draws withsigma.025 before preparation/optimization. All cases therefore retain
fixed noise allocations if a deadline later interrupts. Base/noise/initial
vectors, complete source-deleted/base/initial/best words and local matrices are
saved. L-BFGS-B takes the entire144vector with no bounds/frozen coordinates;
Torchautograd supplies the full gradient. maxiter80,maxfun4000,maxls20,
ftol1e-15,gtol1e-10; one8second bounded fit percase,102secondglobal optimizer
deadline,110secondinternalalarm,120externalcap,one thread,zero restarts.

The exclusive marker forbids repeat execution. Code/family/checker/source/
topology/config/runtime/dependency hashes and active raw-config ledger plus
board ownership are guarded before RNG/matrix construction. The BaseException
TotalDeadline bypasses percase recovery, disables the alarm and preserves
available bestvectors/words without a new fullword matrix evaluation. Serializer
still builds one-qubit matrices for checkpointing; it does not evaluate a
fullword objective. Ordinary percase exceptions are retained as terminations.

No static blocking issue found. The independently reserved preflight133 is
required before root dispatch to verify runtime matrix, phase, serializer and
gradient consistency at all12base/initial points. Fullyfree locals can break
u/v preservation, so ROUTER_REPORT's literal-tail obstruction does not exclude
this distinct fit hypothesis. Every numerical pass needs independent saved-list
validation and exact certification. Any failed bounded fit excludes neither a
word globally nor the11CX goal. Adversarial review used the local protocol;
concurrent work and accepted12 unchanged. Root owns ledger execution/commits.
