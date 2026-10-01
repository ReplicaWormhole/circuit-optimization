# Targeted eigenlabel swaps — run 355

Command: `python3 search13_fresh_chain_label_swaps.py`.

Starting from run339 trial1, choose the six largest unordered off-diagonal transitions between different nearest-root sectors. Swap the eigenlabels of those two output rows, retaining multiplicities(6,3,4,3). Each proposed assignment receives a full168-parameter warm fit of UV=DU for400 iterations, followed by free-label ordinary off-diagonal polish for300 iterations. The frozen plan records source hash, schedule, labels, pairs and iteration limits.

All12 saved native13-CNOT candidates independently fail. Case4 (rows10,14) reaches ordinary loss1.47e-14 and maximum entry1.69e-7: a strong numerical nearhit, still above checker tolerance1e-9. Four other polishes return to the source basin and one reaches a worse basin. This candidate requires further numerical refinement and exact validation. Fixed-label and ordinary objective values are recorded separately; the free polish does not preserve the swapped assignment by assumption.

Root independently reconstructed every saved gate list, checked both Frobenius objectives within1e-11, maximum entries within1e-12, CNOT counts, label swaps and multiplicities. Syntax and ledger audit pass. Adversarial review: no confirmed defect. These six warm starts do not exhaust labels, eigenbasis freedom, starts or architectures, and provide no lower bound. Saved scripts, plan, candidates and result are intentional artifacts.
