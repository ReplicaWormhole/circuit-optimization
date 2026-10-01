"""Exact swap-symmetry certificate and sound opposite-pair support screen.

For real Pauli coefficients, spectral-path zero implies swap symmetry via
six square differences in the exact Gram row space. A traceless involution
on this symmetric two-qubit space has trace(SWAP A) = +/-2: the antisymmetric
singlet is one-dimensional with eigenvalue +/-1. Thus trace(V4^2 A)=+/-4,
contradicting the proper-support half-shift trace lemma for an input axis.
"""
from collections import Counter
from itertools import product
import json
from pathlib import Path
import sympy as sp
import global_lower_bound7_schedule_filter as prior
import global_lower_bound7_adjacent_pair_support_filter as adjacent

ROOT=Path(__file__).resolve().parent

def cone_independent(schedule,j):
    first=next(t for t,e in enumerate(schedule) if j in e)
    reached=1<<j
    for edge in schedule[first+1:]:
        mask=sum(1<<w for w in edge)
        if reached & mask: reached |= mask
    return {w for w in range(4) if reached & (1<<w)}

def main():
    data=json.loads((ROOT/'global_lower_bound7_opposite_pair_general_gram_result.json').read_text())
    gram=sp.Matrix(data['gram']); null=gram.nullspace()
    assert len(null)==45
    labels=data['labels']; monomials=list(map(tuple,data['monomials']))
    certificates=[]
    for i,label in enumerate(labels):
        swapped=label[2]+label[1]+label[0]+label[3]
        k=labels.index(swapped)
        if k<=i:continue
        row=sp.zeros(1,len(monomials))
        row[0,monomials.index((i,i))]=1
        row[0,monomials.index((k,k))]=1
        row[0,monomials.index((i,k))]=-2
        assert all((row*v)[0]==0 for v in null)
        certificates.append({'labels':[label,swapped],
                             'square_difference_in_gram_row_space':True})
    assert len(certificates)==6
    counts=Counter(); examples={}; remaining=Counter()
    for schedule in product(prior.EDGES,repeat=6):
        degree=prior.degrees(schedule)
        if min(degree)<2 or not prior.all_forward_cones_full(schedule):continue
        if prior.tail_rejection(schedule,degree) is not None:continue
        cones=[adjacent.first_commuting_cone(schedule,j)[0] for j in range(4)]
        assert cones==[cone_independent(schedule,j) for j in range(4)]
        if any(len(c)==2 and tuple(sorted(c)) not in prior.OPPOSITE for c in cones):continue
        key=''.join(map(str,sorted(degree,reverse=True)))
        witnesses=[{'wire':j,'pair':sorted(c)} for j,c in enumerate(cones)
                   if len(c)==2 and tuple(sorted(c)) in prior.OPPOSITE]
        if witnesses:
            counts['excluded_opposite_pair']+=1
            examples.setdefault(key,{'schedule':[list(e) for e in schedule],
                                    'witnesses':witnesses})
        else:
            counts['survive']+=1;remaining[key]+=1
    assert sum(counts.values())==3208
    result={'gram_exact_rank':120-len(null),'swap_symmetry_certificates':certificates,
            'prior_survivors':3208,'partition':dict(counts),
            'remaining_by_degree':dict(remaining),'excluded_examples':examples,
            'scope':'Sound necessary screen; surviving schedules are not feasible circuit certificates.'}
    (ROOT/'global_lower_bound7_opposite_pair_support_certificate_result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,sort_keys=True))

if __name__=='__main__':main()
