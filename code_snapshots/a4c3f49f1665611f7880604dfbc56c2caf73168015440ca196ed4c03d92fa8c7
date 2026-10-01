"""Exact feasible commuting pair in the fixed-AA centralizer test."""

from itertools import permutations
import json
from pathlib import Path

import sympy as sp

ROOT=Path(__file__).resolve().parent
I=sp.eye(2)
X=sp.Matrix([[0,1],[1,0]])
Y=sp.Matrix([[0,-sp.I],[sp.I,0]])
Z=sp.diag(1,-1)
P={'I':I,'X':X,'Y':Y,'Z':Z}


def pauli(label):
    return sp.kronecker_product(*(P[a] for a in label))


def main():
    central=json.loads((ROOT/'global_lower_bound7_joint_aa_centralizer_result.json').read_text())
    a_labels=('ZXII','XIIZ','YXIZ')
    b_labels=('IIXZ','ZIIX','ZIXY')
    assert all(label in central['second_overapproximate_basis_labels'] for label in b_labels)
    a_num=sum((pauli(label) for label in a_labels),sp.zeros(16))
    b_num=sum((pauli(label) for label in b_labels),sp.zeros(16))
    assert a_num*a_num==3*sp.eye(16)
    assert b_num*b_num==3*sp.eye(16)
    assert a_num*b_num==b_num*a_num
    shift=sp.zeros(16)
    for state in range(16):
        shift[((state&1)<<3)|(state>>1),state]=1
    ms=[sum(((-sp.I)**(k*r)*shift**r for r in range(4)),sp.zeros(16))
        for k in range(4)]
    assert all(sp.simplify(ms[i]*a_num*ms[j]*a_num*ms[k])==sp.zeros(16)
               for i,j,k in permutations(range(4),3))
    assert all(sp.simplify(ms[i]*b_num*ms[j]*b_num*ms[k])==sp.zeros(16)
               for i,j,k in permutations(range(4),3))
    result={'schedule':central['schedule'],
            'first_image_pauli_labels':list(a_labels),
            'second_image_pauli_labels':list(b_labels),
            'normalization':'divide each sum by sqrt(3)',
            'first_involution':True,'second_involution':True,
            'commute_exactly':True,'first_spectral_path':True,
            'second_spectral_path':True,
            'interpretation':'feasible joint necessary system; not a shift diagonalizer'}
    (ROOT/'global_lower_bound7_joint_aa_pair_witness_result.json').write_text(
        json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,sort_keys=True))


if __name__=='__main__':main()
