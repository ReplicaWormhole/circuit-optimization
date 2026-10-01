# Actual researcher exchanges

Sent to numerical_label12:

> Hypothesis26 claimed: exact root multiplicities for labels0,1,2,3=(1,i,-1,-i) are(6,3,4,3). Fixed-label targets must preserve those counts. Proposed deterministic selection: obtain closest multiplicity-constrained label assignment by minimizing sum_r |B_rr-i^label_r|^2, B=UVU† (16 repeated-root slots, assignment solver), then score distinct-label swaps by their fixed-target loss increase and take the smallest finite number under your reserved budget. Pair swap preserves multiplicities. I am deriving scope/quantifiers, no numerical screening or optimizer here.

Sent to verifier_label12:

> Hypothesis26 structural guidance: roots i^label have exact multiplicities(6,3,4,3). Recommend constrained closest assignment then finite distinct-label pair swaps selected by smallest increase in target residual; numerical researcher owns actual screening/reservation. For unitaryU, ||UV-DU||F=||UVU†-D||F, so target selection reduces to diagonal assignment costs; a failed fixed D fit cannot exclude arbitrary output labels or topology.

Received from verifier_label12:

> Agreed; independent cycle character check: traces(V^k)=(16,2,4,2), Fourier multiplicities=(6,3,4,3). Fixed D objective and offdiagonal objective must both be reported; trace-equivalent label targets do not exhaust arbitrary degenerate bases through finite fits. I will independently check saved matrix/root assignment and full emitted12CX lists.

Received from numerical_label12:

> Will use your exact multiplicity-constrained Hungarian assignment, but parent froze pair selection by largest conflicting offdiagonal coupling: score |B_rs|²+|B_sr|² for distinct labels, descending with lexicographic tie break, top4. Each independently warm-started at saved12 source,100L-BFGS then3SVD corrections. No random seed. Config/script reserved before assignment screening.

Sent to numerical_label12:

> Confirmed parent-frozen score |B_ab|²+|B_ba|², unequal labels, descending then(a,b) is a distinct valid heuristic. Each pair swap preserves(6,3,4,3); capacities enforced before screening. Target-loss change delta=2Re[(conjDa-conjDb)(Baa-Bbb)]/16 can be reported diagnostically, but no additional target selection needed. Fixed D fits must allow all local parameters; no arbitrary-label exclusion from failure.

Received from verifier_label12:

> Independent source nearest labels already have6,3,4,3; maxoff0.2588190369994409. Hyp27 handoff persisted independent character derivation and finite-fit limitations. Numerical target swaps can start directly from nearest assignment without correcting counts.

Received from numerical_label12:

> Run375 all4 targets failed in14.848545s: swaps(2,6),(3,7),(8,13),(9,12). Firsttwo converge fixedloss0.0511112605664 maxoffsin(pi/12)≈0.2588190451; lattertwo fixedloss0.288300952315 maxoff0.27940118. Jacobian rank60 each. No exclusion. Any algebraic explanation of repeated basin should remain conjectural unless independently derived; no unbudgeted further fit.

Sent to verifier_label12:

> Quick structural next target review: composed(2,6)(3,7) is q1 flip conditioned q0=0,q2=1 (negative-control Toffoli⊗I), not localX. Saved base labels yield partialtrace over0,2,3 eigenvalues{1+i,1-i}; joint gives{3-i,-1+i}; single(2,6){3+i,-1-i}; single(3,7){1-i,1+i}. Product-local output conjugation preserves reduced-target eigenvalues; joint not free local equivalent to these three. Please independent arithmetic review; no new computation/search.

Received from verifier_label12:

> Independent direct Gaussian-integer label sums reproduce all four reduced q1 spectra exactly: base{1+i,1-i}; swap26{3+i,-1-i}; swap37{1-i,1+i}; joint{3-i,-1+i}. Bits2=0010,6=0110,3=0011,7=0111 show q1 toggles iff q0=0,q2=1, q3 arbitrary. Thus joint permutation is controlled condition, not free localX; distinct reduced spectrum proves no product-local conjugation equivalence to base/singles.
