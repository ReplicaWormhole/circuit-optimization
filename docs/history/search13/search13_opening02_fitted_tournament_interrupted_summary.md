# Interrupted opening02 tournament (run 351)

Original handle 1646 became unavailable; OS inspection found no script process. The frozen 65 proposals (five excluded, 60 retained) and seven saved short candidates, proposals 5 through 11, were preserved. No result file or optimizer iteration/convergence metadata survived. Every saved gate list independently failed `check_circuit.py` with 13 CNOTs. Saved-gate NumPy losses and hashes are in `search13_opening02_fitted_tournament_interrupted_validation.json`. Best recovered loss is 0.014839314687471914, maximum entry 0.16489445225373636 (proposal 7), worse than source.

Run 351 is inconclusive due to interruption, not a topology exclusion. A separate recovery run reuses the seven saved gate lists as short-fit endpoints, honestly reconstructing losses while leaving unavailable optimization metadata absent. It fits remaining 53 proposals and selects six by fitted direct loss. Original files remain untouched.
