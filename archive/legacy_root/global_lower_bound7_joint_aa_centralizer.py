"""Exact Pauli-space centralizer for two AA selected axes in one schedule.

The first image is the canonical exact run248/253 witness on wires 0,1,3.
The second selected wire is 3, with partners 2 then 0.  Without using a
second collective gauge, its image lies in the full 48D Pauli span having
nonidentity on wire 3 and arbitrary factors on wires 2 and 0.
"""

import json
from pathlib import Path

import sympy as sp

ROOT=Path(__file__).resolve().parent
P={('X','Y'):(sp.I,'Z'),('Y','X'):(-sp.I,'Z'),
   ('Y','Z'):(sp.I,'X'),('Z','Y'):(-sp.I,'X'),
   ('Z','X'):(sp.I,'Y'),('X','Z'):(-sp.I,'Y')}


def product(a,b):
    phase=sp.Integer(1)
    out=[]
    for u,v in zip(a,b):
        if u=='I':out.append(v)
        elif v=='I':out.append(u)
        elif u==v:out.append('I')
        else:
            factor,letter=P[(u,v)]
            phase*=factor
            out.append(letter)
    return phase,''.join(out)


def main():
    selected=json.loads((ROOT/'global_lower_bound7_joint_aa_schedule_select_result.json').read_text())
    representative=selected['examples_by_second_geometry']['j_center_adjacent_v']
    assert representative['schedule']==[[0,1],[2,3],[1,2],[0,1],[2,3],[0,3]]
    first=('ZXII','XIIZ','YXIZ')
    second_labels=[''.join((q0,'I',q2,q3)) for q0 in 'IXYZ'
                   for q2 in 'IXYZ' for q3 in 'XYZ']
    assert len(second_labels)==48 and len(set(second_labels))==48
    all_labels=[''.join(letters) for letters in
                __import__('itertools').product('IXYZ',repeat=4)]
    index={label:i for i,label in enumerate(all_labels)}
    comm=sp.zeros(256,48)
    for col,b in enumerate(second_labels):
        for a in first:
            phase,label=product(a,b)
            reverse,_=product(b,a)
            assert reverse==-phase or reverse==phase
            if reverse==-phase:
                comm[index[label],col]+=sp.simplify((phase-reverse)/(2*sp.I))
    null=comm.nullspace()
    assert all(comm*v==sp.zeros(256,1) for v in null)
    result={'schedule':representative['schedule'],
            'first_selected_wire':0,'second_selected_wire':3,
            'first_axis_pauli_labels':list(first),
            'second_overapproximate_basis_labels':second_labels,
            'commutator_rank':48-len(null),'centralizer_dimension':len(null),
            'centralizer_basis_coefficients':[[str(v[i]) for i in range(48)] for v in null]}
    (ROOT/'global_lower_bound7_joint_aa_centralizer_result.json').write_text(
        json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in
                      ('second_overapproximate_basis_labels','centralizer_basis_coefficients')},
                     sort_keys=True))


if __name__=='__main__':main()
