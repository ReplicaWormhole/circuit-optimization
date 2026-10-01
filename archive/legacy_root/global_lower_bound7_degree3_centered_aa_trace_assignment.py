"""Exact eigenvalue-count assignment check for the AA witness trace moments."""

import json
from itertools import product
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent


def main():
    witness = json.loads((ROOT / 'global_lower_bound7_degree3_centered_aa_witness_result.json').read_text())
    target = tuple(sp.sympify(v) for v in witness['trace_V_powers_1_2_3'])
    eigenvalues = (sp.Integer(1), sp.I, sp.Integer(-1), -sp.I)
    multiplicities = (6,3,4,3)
    pz = 1/sp.sqrt(3)
    assignments = []
    for left in product(*(range(n+1) for n in multiplicities)):
        if sum(left) != 8:
            continue
        difference = tuple(2*a-n for a,n in zip(left,multiplicities))
        moments = tuple(sp.simplify(pz*sum(d*lam**power
                                          for d,lam in zip(difference,eigenvalues)))
                        for power in (1,2,3))
        if all(sp.simplify(a-b)==0 for a,b in zip(moments,target)):
            assignments.append({'bit_zero_eigenvalue_counts':dict(zip(('1','i','-1','-i'),left)),
                                'bit_one_eigenvalue_counts':dict(zip(('1','i','-1','-i'),
                                                                     tuple(n-a for a,n in zip(left,multiplicities))))})
    result={'witness_trace_moments':[str(v) for v in target],
            'input_axis_z_component':str(pz),
            'matching_bit_partition_count':len(assignments),
            'example':assignments[0] if assignments else None}
    (ROOT/'global_lower_bound7_degree3_centered_aa_trace_assignment_result.json').write_text(
        json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,sort_keys=True))


if __name__=='__main__':
    main()
