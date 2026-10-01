# Round23 preparation adversarial review

Artifact types: code, computation protocol, workflow. Reviewed using
/home/lev/.codex/adversarial_review.md. No confirmed static preparation issue
remains. Root clarified no-ranking as no adaptive allocation/deepening/restarts;
the current frozen pack adds deterministic global min-loss closeout summary
after both fixed-budget case records. Earlier unreserved no-summary pack is
preserved in superseded_no_global_summary; it was never numerically evaluated.

Actual metadata-only checks: AST parse fit.py/family.py; canonical raw config;
all eleven current pack hashes and all immutable origin hashes; source best132
unchanged; literalnewbase retains all old U3 local gates and has prescribed
54gate/10CX topology; seed/sigmas/two-start counts and250iter10000calls,
15/35/45/60second budget values; V4/parent416/agentnumerical10reorder/board148;
focused deterministic-summary diff. All passed. Runtime dependency metadata
was read by importing fit.py only, with no family import/evaluation or main.
No matrix, RNG, AD, optimizer, scientific run reservation or shared code edit.

Static optimizer inspection: two noise draws use one PCG64 generator and are
allocated before either start. Both complete initial states are written before
optimization. Each case retains132free coordinates and single frozen budget;
no restart/deepening loop. Guard/source/config/board/code/runtime checks precede
numerical work. Exception salvage retains current best endpoint histories.

Residual risk: syntax and static review do not establish numerical family,
serializer, gradient or deadline correctness. Independent reserved preflight
and by-hand algebraic review remain mandatory before root dispatch. The fixed
classical-sector/computational-x obstruction does not exclude full-local words.
No exactness, topology exclusion, lower bound or incumbent change follows.

Cleanup: new directory and superseded freeze are intentional research evidence;
Python caches are ignored. Root owns execution, shared files and commits.
