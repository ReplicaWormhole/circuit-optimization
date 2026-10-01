"""Exact rank-one spectral-path and involution ideal in 12D centralizer."""

from collections import defaultdict
import json
from pathlib import Path

import sympy as sp

from global_lower_bound7_joint_aa_centralizer import product as pauli_product

ROOT=Path(__file__).resolve().parent


def primitive(expr,symbols):
    poly=sp.Poly(sp.expand(expr),*symbols,domain=sp.QQ)
    return None if poly.is_zero else poly.primitive()[1].as_expr()


def main():
    data=json.loads((ROOT/'global_lower_bound7_joint_aa_centralizer_gram_result.json').read_text())
    gram=sp.Matrix(data['gram'])
    monomials=[tuple(m) for m in data['monomials']]
    null=gram.nullspace()
    assert len(null)==9
    eliminated=[i for i in range(12)
                if all(v[monomials.index((i,i))]==0 for v in null)]
    assert [data['basis_labels'][i] for i in eliminated]==['IIIZ','IIZZ','ZIZX','ZIZY']
    remaining=[i for i in range(12) if i not in eliminated]
    symbols=sp.symbols('x:'+str(len(remaining)))
    x=[sp.Integer(0)]*12
    for i,v in zip(remaining,symbols):x[i]=v
    m=sp.Matrix([x[i]*x[j] for i,j in monomials])
    spectral=set()
    for expr in gram*m:
        p=primitive(expr,symbols)
        if p is not None:spectral.add(p)
    square=defaultdict(lambda:sp.Integer(0))
    for i,label in enumerate(data['basis_labels']):
        if x[i]==0:continue
        square['IIII']+=x[i]**2
        for j in range(i+1,12):
            if x[j]==0:continue
            phase,out=pauli_product(label,data['basis_labels'][j])
            if phase in (1,-1):square[out]+=2*phase*x[i]*x[j]
    square['IIII']-=1
    involution=set()
    for expr in square.values():
        p=primitive(expr,symbols)
        if p is not None:involution.add(p)
    equations=sorted(spectral|involution,key=str)
    gb=sp.groebner(equations,*symbols,order='grevlex',domain=sp.QQ)
    result={'exact_gram_rank':78-len(null),'exact_gram_nullity':len(null),
            'forced_zero_labels':[data['basis_labels'][i] for i in eliminated],
            'remaining_labels':[data['basis_labels'][i] for i in remaining],
            'spectral_equation_count':len(spectral),
            'involution_equation_count':len(involution),
            'groebner_basis':[str(p.as_expr()) for p in gb.polys],
            'has_complex_solution':not gb.contains(sp.Integer(1))}
    (ROOT/'global_lower_bound7_joint_aa_centralizer_certificate_result.json').write_text(
        json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,sort_keys=True))


if __name__=='__main__':main()
