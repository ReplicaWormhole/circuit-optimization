# Parity-prefix native helper and independent verification

Hypothesis54. Owned native_helpers.py adds differentiable F(a,b)=exp(i aXX+
i bYY) kernel, embedding and exactlytwo-CX serialization. Chronological
Rx_c,t(-pi/2),CXct,Rx_c(-2a),Rz_t(-2b),CXct,Rx_c,t(pi/2) implementsF;
this is an upper compilation count, not minimum for special parameters.
Existing matchgate identity supplies analytic justification. No sharednative
checker modification or new dependency. Focused2tests pass15angle/wirecases,
including nonadjacent/reversed embedding and AD gradientcheck.

Before reservation/fit, helper SHA992aa07290b35bcc0ad0860e1f7053f790b0d72e45c477e12ebd13a902f734b3,
scriptSHA dd7a738bce56d137f2f1ab9bdb67c4a90137935f7964a34065e341da3ed0e1dc
and configSHA74e348841156ad1a29117cb43b81944754df60e4bfa9efa9bcb936ca371beb0f
were independently validated and directapproval persisted. Restrictedfamily
has sevenfree local layers,84axis+8F coordinates, nointeriorprefixlocals;
L0 CX01,CX21,CX31 L1 F02 L2 F13 L3 CX12 L4 F01 L5 F23 L6.
ArbitraryL0 means prefix parity is currentbit parity, not originalHamming
weight; downstreamF neednotpreserve tag. This is a hypothesis, not exact
construction. TwoGaussian.3 seeds9261101/02,500iterations each240seconds
totalone-thread frozen;385 executed2.228609seconds.

IndependentSciPyexpm XX/YY generators and explicitbit embeddings verify
native/compiled fullmatrices differ3.86e-16/4.17e-16. Native4CX4F costs
8entanglers in preciselydefined hybridset; compiledlists12CX. Both invalid:
maxoff0.9438619418015138/0.7071067866282266; offloss0.5922239996302485/
0.5000000000000019. Nearestroot counts fail required6,3,4,3, consistent
with non-diagonalization. No hybrid8gate or12CX upperbound promoted.
Exactseedinitialization and fullnative chronology independently parse;
saved92axis/F coordinates reconstruct matrices within7.47e-16/6.66e-16
upglobalphase. Independentunitarity below1.89e-15. No exact certification
needed for rejection; parentcertificates not duplicated.

Commands test_helpers.py, parameter_audit.py, check_candidates.py with both
native/compiledpairs pass. One checker import-sideeffect found during audit
was fixed using main guard and affectedchecks rerun. Adversarial review no
remaining confirmeddefect; failure excludes neither restrictedfamily nor
broader12CXarchitectures. Exact13incumbent and6<=Cmin<=13 unchanged. New
source/tests/JSON/report/chats intentionalretained; coordinatorownscommit,
no push.
