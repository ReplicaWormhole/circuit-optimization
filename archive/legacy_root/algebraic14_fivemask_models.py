"""Enumerate five-mask phase-polynomial supports on the pi/8 angle grid."""

import json
from itertools import combinations
from pathlib import Path

from algebraic14_fourmask_models import model_for_support
from algebraic14_phase_support import NONLOCAL


ROOT = Path(__file__).parent


def main():
    models = []
    tested = 0
    for support in combinations(NONLOCAL, 5):
        tested += 1
        coefficients = model_for_support(support)
        if coefficients is not None:
            models.append({"support": support,
                           "coefficients_pi_over_8": coefficients})
    path = ROOT / "algebraic14_fivemask_models.json"
    path.write_text(json.dumps(models, indent=2) + "\n")
    print(json.dumps({"tested_supports": tested,
                      "satisfiable_supports": len(models),
                      "models_file": path.name,
                      "supports": [record["support"] for record in models]},
                     indent=2))


if __name__ == "__main__":
    main()
