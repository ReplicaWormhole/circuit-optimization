"""Bounded exact feasibility attempt for the adjacent/adjacent centered tail.

Use the actual image's common Bloch axis on the first partner k, fix the
last partner l axis to Z by a collective rotation, and fix the remaining
azimuth of the k axis by a collective Z rotation.  The resulting 14-real-
variable overapproximation is tested against the exact spectral-path Gram
and Hermitian involution equations.
"""

from collections import defaultdict
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
PAULI_PRODUCT = {
    ('X','Y'): (sp.I,'Z'), ('Y','X'): (-sp.I,'Z'),
    ('Y','Z'): (sp.I,'X'), ('Z','Y'): (-sp.I,'X'),
    ('Z','X'): (sp.I,'Y'), ('X','Z'): (-sp.I,'Y'),
}


def product_labels(a, b):
    phase = sp.Integer(1)
    result = []
    for u,v in zip(a,b):
        if u == 'I':
            result.append(v)
        elif v == 'I':
            result.append(u)
        elif u == v:
            result.append('I')
        else:
            factor, letter = PAULI_PRODUCT[(u,v)]
            phase *= factor
            result.append(letter)
    return phase, ''.join(result)


def primitive(expr, variables):
    poly = sp.Poly(sp.expand(expr), *variables, domain=sp.QQ)
    if poly.is_zero:
        return None
    return poly.primitive()[1].as_expr()


def main():
    data = json.loads((ROOT / 'global_lower_bound7_degree3_centered_gram_adjacent_adjacent_result.json').read_text())
    assert data['ordered_triple_j_k_l'] == [0,1,3]
    a = sp.symbols('a:3')
    b = sp.symbols('b:3')
    c = sp.symbols('c:3')
    d = sp.symbols('d:3')
    s,t = sp.symbols('s t')
    variables = (*a,*b,*c,*d,s,t)
    by_mu = {letter:i for i,letter in enumerate('XYZ')}
    x = []
    for label in data['basis_labels']:
        mu, k, _, ell = label
        i = by_mu[mu]
        if k == 'I':
            x.append(a[i] if ell == 'I' else c[i])
        elif k == 'X':
            x.append((b[i] if ell == 'I' else d[i])*s)
        elif k == 'Z':
            x.append((b[i] if ell == 'I' else d[i])*t)
        else:
            assert k == 'Y'
            x.append(sp.Integer(0))
    monomials = [tuple(pair) for pair in data['monomials']]
    m = sp.Matrix([x[i]*x[j] for i,j in monomials])
    gram = sp.Matrix(data['gram'])
    equations = set()
    for expression in gram*m:
        poly = primitive(expression, variables)
        if poly is not None:
            equations.add(poly)
    spectral_count = len(equations)
    square_coefficients = defaultdict(lambda: sp.Integer(0))
    for i,label in enumerate(data['basis_labels']):
        if x[i] == 0:
            continue
        square_coefficients['IIII'] += x[i]**2
        for j in range(i+1,len(x)):
            if x[j] == 0:
                continue
            phase, output = product_labels(label,data['basis_labels'][j])
            if phase in (1,-1):
                square_coefficients[output] += 2*phase*x[i]*x[j]
    square_coefficients['IIII'] -= 1
    involution = set()
    for expression in square_coefficients.values():
        poly = primitive(expression, variables)
        if poly is not None:
            involution.add(poly)
    equations.update(involution)
    equations.add(s*s+t*t-1)
    metadata = {'shape':'adjacent_adjacent','variables':[str(v) for v in variables],
                'spectral_equation_count':spectral_count,
                'involution_equation_count':len(involution),
                'total_distinct_equations':len(equations),
                'status':'groebner_running'}
    path = ROOT / 'global_lower_bound7_degree3_centered_aa_common_axis_result.json'
    path.write_text(json.dumps(metadata,indent=2,sort_keys=True)+'\n')
    print(json.dumps(metadata,sort_keys=True),flush=True)
    gb = sp.groebner(sorted(equations,key=str), *variables, order='grevlex', domain=sp.QQ)
    result = {**metadata, 'status':'groebner_finished',
              'groebner_basis':[str(p.as_expr()) for p in gb.polys],
              'has_complex_solution':not gb.contains(sp.Integer(1))}
    path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='groebner_basis'},sort_keys=True))


if __name__ == '__main__':
    main()
