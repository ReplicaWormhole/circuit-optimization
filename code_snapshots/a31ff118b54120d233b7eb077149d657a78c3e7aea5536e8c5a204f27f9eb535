"""Exact finite joint trace-moment assignment test for AA commuting pair.

For A=U p_0 U† and B=U p_3 U†, with p_j a local Bloch axis having
input-computational Z component z_j, diagonality of D=U†V4U gives
tr(V4**r A)=z_0 tr(D**r Z_0), similarly for B, and
tr(V4**r AB)=z_0*z_3 tr(D**r Z_0 Z_3).
Enumerate all four input-bit-pair blocks of size four with the exact
eigenvalue multiplicities (6,3,4,3).
"""

from itertools import product
import json
from pathlib import Path

import sympy as sp

ROOT=Path(__file__).resolve().parent
I=sp.eye(2)
X=sp.Matrix([[0,1],[1,0]])
Y=sp.Matrix([[0,-sp.I],[sp.I,0]])
Z=sp.diag(1,-1)
P={'I':I,'X':X,'Y':Y,'Z':Z}
EIGENVALUES=(1,1j,-1,-1j)
TOTAL=(6,3,4,3)
BLOCKS=((0,0),(0,1),(1,0),(1,1))


def pauli(label):
    return sp.kronecker_product(*(P[c] for c in label))


def gaussian(value):
    return sp.Integer(int(round(value.real)))+sp.I*sp.Integer(int(round(value.imag)))


def signed_moment(counts,bit_mask,power):
    value=0j
    for (b0,b3),row in zip(BLOCKS,counts):
        sign=(-1)**((b0 if bit_mask&1 else 0)+(b3 if bit_mask&2 else 0))
        value+=sign*sum(n*(lam**power) for n,lam in zip(row,EIGENVALUES))
    return gaussian(value)


def main():
    a_num=sum((pauli(label) for label in ('ZXII','XIIZ','YXIZ')),sp.zeros(16))
    b_num=sum((pauli(label) for label in ('IIXZ','ZIIX','ZIXY')),sp.zeros(16))
    shift=sp.zeros(16)
    for state in range(16):
        shift[((state&1)<<3)|(state>>1),state]=1
    target={
        'A':tuple(sp.simplify(sp.trace(shift**r*a_num)/sp.sqrt(3)) for r in (1,2,3)),
        'B':tuple(sp.simplify(sp.trace(shift**r*b_num)/sp.sqrt(3)) for r in (1,2,3)),
        'AB':tuple(sp.simplify(sp.trace(shift**r*a_num*b_num)/3) for r in (1,2,3)),
    }
    rows=[row for row in product(range(5),repeat=4) if sum(row)==4]
    counts_checked=0
    marginal_pass=0
    matches=[]
    for row00,row01,row10 in product(rows,repeat=3):
        row11=tuple(t-a-b-c for t,a,b,c in zip(TOTAL,row00,row01,row10))
        if any(n<0 for n in row11) or sum(row11)!=4:
            continue
        counts=(row00,row01,row10,row11)
        counts_checked+=1
        s0=signed_moment(counts,1,1)
        s3=signed_moment(counts,2,1)
        if s0==0 or s3==0:
            continue
        z0=sp.simplify(target['A'][0]/s0)
        z3=sp.simplify(target['B'][0]/s3)
        if z0.is_real is not True or z3.is_real is not True:
            continue
        if abs(float(z0))>=1 or abs(float(z3))>=1:
            continue
        if any(sp.simplify(z0*signed_moment(counts,1,r)-target['A'][r-1])!=0
               or sp.simplify(z3*signed_moment(counts,2,r)-target['B'][r-1])!=0
               for r in (2,3)):
            continue
        marginal_pass+=1
        if all(sp.simplify(z0*z3*signed_moment(counts,3,r)-target['AB'][r-1])==0
               for r in (1,2,3)):
            matches.append({'block_counts':{''.join(map(str,block)):list(row)
                                             for block,row in zip(BLOCKS,counts)},
                            'z0':str(z0),'z3':str(z3)})
    result={'target_trace_moments':{k:[str(v) for v in values] for k,values in target.items()},
            'block_row_options':len(rows),'valid_total_assignments_checked':counts_checked,
            'assignments_passing_single_axis_moments':marginal_pass,
            'assignments_passing_joint_product_moments':len(matches),
            'example_joint_assignment':matches[0] if matches else None,
            'all_joint_assignments':matches}
    (ROOT/'global_lower_bound7_joint_aa_trace_moments_result.json').write_text(
        json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='all_joint_assignments'},sort_keys=True))


if __name__=='__main__':main()
