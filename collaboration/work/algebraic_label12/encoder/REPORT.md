# Fixed orbit encoder: ANF and actual reversible program

Hypothesis46; one fixed construction, no optimization or symbolic search.
Definition-level verification is `python3 collaboration/work/algebraic_label12/encoder/derive.py`.
Artifact types are math, computation and a runnable verification script.
Qubit0 is MSB; P|x>=|p(x)> with the original orbit-to-address convention.

The exact sixteen-row mapping, input order0..15, is
[0,4,7,8,6,2,11,12,5,9,3,13,10,14,15,1]. It has cycles
(1,4,6,11,13,14,15), (2,7,12,10,3,8,5), with0 and9 fixed.

Write inputs as a=q0,b=q1,c=q2,d=q3, arithmetic over F2. The output ANFs are

* q0'=cd+bc+bcd+ad+acd+ab+abd+abc;
* q1'=a+b+c+d;
* q2'=c+cd+b+bc;
* q3'=c+cd+a+ac+acd+ab+abc.

These are Boolean identities, not independent in-place instructions. Assigning
one output overwrites input values needed by others. The explicit reversible
implementation below handles that issue without temporary or ancillary bits.

For each nontrivial cycle(a0,a1,...,a6), apply chronologically the row swaps
(a0,a1),(a0,a2),...,(a0,a6). This sends each cycle entry to its listed
successor. For each row swap(a,b), make a fixed Gray path changing differing
bits in q0,q1,q2,q3 order. If adjacent path swaps are e1,...,eh, apply
e1,...,eh,e(h-1),...,e1; this exchanges endpoints and restores all interior
rows. Each adjacent swap is an X on the changed bit controlled on fixed
values of the other three bits. `encoder.json` contains all44 chronological
three-controlX gates with exact positive/negative control values. The script
evaluates these gates on every input and reproduces P, rather than treating
the ANF as an executable assignment. It also evaluates every ANF on all16
inputs: all64 output-bit identities pass.

This fixed gate list has44 native three-controlX primitives. It is not a44CX
construction. Each three-controlX is local-H conjugate to a four-bit pi-phase
projector. Expanding that projector in commuting Z strings costs34CX using
parity ladders: six weight-two terms cost12, four weight-three terms cost16,
one weight-four term costs6; constants and local terms cost zeroCX. Thus this
explicit encoder has certified generic compilation upper1496CX. Combining
with the previous explicit conditional Fourier upper268 gives1764CX for the
full abstract diagonalizer. These are loose constructive upper bounds, not
minimal counts, an improved incumbent, or a serialized nativeCX circuit.

The parity output q1' suggests a useful alternative synthesis seed: first
compute it using CX0→1,CX2→1,CX3→1, then synthesize the other outputs while
preserving or restoring q1. That three-CX prefix is exact, but no remaining
short program is asserted here. The straight-line44gate program is the
completed constructive deliverable. No cancellation with conditionalF4/H2
was demonstrated; any optimized program needs its own frozen search budget
and independent full16input check. Matching phase-ladder parity networks at
the encoder/Fourier boundary is an actionable joint-synthesis hypothesis,
not a gate saving already obtained.

Adversarial review: cycle direction, Gray restoration and overwritten-input
issue explicitly checked; both direct truth-table and gate evaluations pass;
nativeCCCX and compiledCX counts are separated. No confirmed defect. Remaining
challenge is efficient synthesis and any phase-aware cancellation. No
exact13 certificate duplicated. Source/JSON/report/chat files are intentional
retained evidence, coordinator owns the shared commit.

Verifier's independent check_encoder.py reproduces all16 native-program
outputs and all64 ANF bits exactly, and confirms44*34=1496 and full bound1764.
The numerical researcher's separately reserved beam search is a distinct task;
this derivation neither executes it nor promotes approximate programs.

Subsequent direct peer reports for reserved run381 give an exact10step program,
7CX+3CCX with compiled upper25, independently evaluated on all16 inputs. That
separate numerical result supersedes this44CCCX program as an efficient
encoder baseline, not as a whole diagonalizer. Its first three CX operations
3→2,0→1,2→1 compute c+d, then b+a, then the required input parity on q1.
The generic Fourier upper still gives293CX; no comparison near13 is justified.
An unscheduled joint-synthesis idea is cheaper relative-phase CCX replacement,
but it requires tracking accumulated phases, proving they are constant on
cycle orbits or explicitly compensating them in the Fourier blocks. Pure
permutation truth-table agreement alone would not validate that modification.
