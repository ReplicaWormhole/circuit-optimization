"""Bounded symbolic determinant audit for higher-body parity output gauges.

For masks of weight 3 or 4, set the odd-parity row phase to an indeterminate
t (the omitted even-row phase is an irrelevant nonzero global scalar).
Factor the two 16x16 operator-Schmidt realignment determinants of G(t)S.
The symbolic run has a per-determinant wall-time cap; missing factors are
reported as timeouts rather than interpreted as rank obstructions.
"""

import json
import signal
import time
from pathlib import Path

import sympy as s


ROOT = Path(__file__).resolve().parent
MASKS = (7, 11, 13, 14, 15)
CUTS = ((0, 2), (0, 3))
LIMIT_SECONDS = 20


class TimedOut(Exception):
    pass


def on_alarm(_signum, _frame):
    raise TimedOut()


def bit(value, wire):
    return (value >> (3 - wire)) & 1


def realign(unitary, cut):
    other = tuple(q for q in range(4) if q not in cut)
    result = s.zeros(16)
    for output in range(16):
        for input_ in range(16):
            oa = 2 * bit(output, cut[0]) + bit(output, cut[1])
            ia = 2 * bit(input_, cut[0]) + bit(input_, cut[1])
            ob = 2 * bit(output, other[0]) + bit(output, other[1])
            ib = 2 * bit(input_, other[0]) + bit(input_, other[1])
            result[4 * oa + ia, 4 * ob + ib] = unitary[output, input_]
    return result


def main():
    a = s.sqrt(2) / 2
    f = s.Matrix([[1, 0, 0, 0], [0, a, s.I * a, 0],
                  [0, s.I * a, a, 0], [0, 0, 0, 1]])
    suffix = s.kronecker_product(f, f)
    t = s.symbols("t")
    signal.signal(signal.SIGALRM, on_alarm)
    rows = []
    for mask in MASKS:
        gauge = s.diag(*[t if (basis & mask).bit_count() % 2 else 1
                         for basis in range(16)])
        u = gauge * suffix
        for cut in CUTS:
            start = time.monotonic()
            try:
                signal.alarm(LIMIT_SECONDS)
                determinant = s.factor(realign(u, cut).det(method="domain-ge"))
                status = "complete"
                factor = str(determinant)
            except TimedOut:
                status = "timeout"
                factor = None
            finally:
                signal.alarm(0)
            row = {"mask": mask, "cut": list(cut), "status": status,
                   "factor": factor, "seconds": time.monotonic() - start}
            rows.append(row)
            print(json.dumps(row, sort_keys=True), flush=True)
    path = ROOT / "algebraic13_crosspair_parity_symbolic_result.json"
    path.write_text(json.dumps({"row_phase_definition":
                                "even parity 1, odd parity t; physical t=e^(-2ic)",
                                "per_determinant_time_limit_seconds": LIMIT_SECONDS,
                                "rows": rows}, indent=2) + "\n")


if __name__ == "__main__":
    main()
