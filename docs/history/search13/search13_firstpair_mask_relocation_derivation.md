# Relocating one parity phase into the first pair block

This successor to the first-pair endpoint screen permits the omitted phase
mask to be implemented inside the first `02` or `13` two-qubit block and
permits the block's Cartan angles to change. It fixes the later matchgate
architecture of the exact 14-CNOT construction. The five-CNOT prefix's
linear endpoint may differ from the six-CNOT endpoint by **any** invertible
two-bit linear map on the selected pair, not just one CNOT.

Let the six-CNOT endpoint rows be `E6=(8,5,10,1)`. A five-CNOT path has
endpoint rows `E5`. A phase mask `m` omitted from the prefix becomes, after
that prefix, the Pauli `Z` string whose output-wire support solves
`m = XOR(E5[j] for j in support)`. For the phase to sit inside one first pair
block, the exact necessary condition is

`m in {E5[a], E5[b], E5[a] XOR E5[b]}`

for `(a,b)=(0,2)` or `(1,3)`. The first two values are local `Z` phases in
the block's input frame; the XOR value is a `ZZ` phase. The linear mismatch
must also leave the other pair untouched and map the selected pair to itself.
There are 11 distinct pair-local endpoints: six `GL(2,2)` maps on each pair,
with the identity counted once.

Runs 198 and 200 exhaust all `12^5=248832` directed five-CNOT schedules
and all 62 catalogued rational-`pi/8` phase models with four or five
nonlocal masks. The phase model's required support must have exactly one
mask not visited by the path. There are 3,574 paths at pair-local endpoints,
and 196 path/pair/model cases with exactly one missing required mask. None
of these 196 cases satisfies the pair-support condition above: zero `ZZ`
relocations and zero local-`Z` relocations. The omitted mask always involves
at least one wire outside the selected pair after the prefix.

This exact finite screen excludes the stated one-mask relocation ansatz with
the remaining tail held at its exact-14 architecture. Because no case passes
the support condition, changing the first block's Cartan angles cannot help
within that ansatz, and a numerical joint fit would test an impossible
decomposition. It does not constrain a 13-CNOT circuit that also changes the
center layer or the later pair blocks, or a prefix outside the catalogued
phase-polynomial models.

Reproduce the screens and ledger entries with:

```bash
python3 search13_firstpair_mask_relocation.py
python3 search13_firstpair_mask_local_screen.py
python3 experiment_log.py show 198
python3 experiment_log.py show 200
```
