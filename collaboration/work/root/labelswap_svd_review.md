# Thirteen-CNOT numerical candidate review

Artifacts: computation, math, native CNOT gate list.
Source chain order is(01,23,12) repeated four times followed by01, with alternating directions. All168 local coordinates are free; no ancillary wire or hidden entangler occurs. A targeted nearest-root eigenlabel swap at rows10/14 in run355 produced the new nearhit. Run358 reduced it to4.74e-8. Run359 used one truncated-SVD minimum-norm Newton step, retained rank56, cutoff1e-6, and reached numerical roundoff.

Independent NumPy chronological gate multiplication verifies UV-DU with maximum error1.58e-15. The local checker reports exactly13 CNOTs, unitarity error2.00e-15, maximum off-diagonal entry2.00e-15, and correct root multiplicities(6,3,4,3). Qiskit using the existing /tmp/circuit-qiskit-venv independently gives maximum off-diagonal entry1.85e-15; transpilation retains13 CNOTs and error3.17e-15. No package was installed. Default python3 lacks Qiskit; that initial failed invocation was superseded by the venv run and explicitly corrected on the board.

Confirmed reporting error fixed: the first summary of run355 incorrectly said all polishes returned to the old basin. Case4 instead reached a nearhit; case2 reached a worse basin. The run355 ledger note's phrase 'No better basin' is superseded by this correction and run359's passing numerical record; its archived gate list and measured1.69e-7 error were accurate.

Adversarial review: no confirmed mathematical or gate-model defect in the numerical result. Floating-point checks do not constitute an exact upper-bound certificate. The exact incumbent remains14 until symbolic certification; total spin is not diagonalized. Exact-certification hypothesis5 was handed to the algebraic agent, with a bounded12fit budget. Root numerical source and evidence are kept artifacts.
