"""Exact spectral-path Gram and shift traces on fixed-AA 12D centralizer."""

from itertools import combinations_with_replacement, permutations
import json
from pathlib import Path

import numpy as np

from global_lower_bound6_chain_gram import kron_label, rank_mod_prime

ROOT=Path(__file__).resolve().parent


def main():
    data=json.loads((ROOT/'global_lower_bound7_joint_aa_centralizer_result.json').read_text())
    original=data['second_overapproximate_basis_labels']
    labels=[]
    for vector in data['centralizer_basis_coefficients']:
        nonzero=[i for i,v in enumerate(vector) if v!='0']
        assert len(nonzero)==1 and vector[nonzero[0]]=='1'
        labels.append(original[nonzero[0]])
    assert len(labels)==12
    axes=[kron_label(label) for label in labels]
    shift=np.zeros((16,16),dtype=np.complex128)
    for state in range(16):
        shift[((state&1)<<3)|(state>>1),state]=1
    powers=[np.linalg.matrix_power(shift,r) for r in range(4)]
    ms=[sum(((-1j)**(k*r)*powers[r] for r in range(4)),
            np.zeros((16,16),dtype=np.complex128)) for k in range(4)]
    monomials=list(combinations_with_replacement(range(12),2))
    blocks=[]
    for a,b,c in permutations(range(4),3):
        per_triple=[]
        for i,j in monomials:
            block=ms[a]@axes[i]@ms[b]@axes[j]@ms[c]
            if i!=j:
                block+=ms[a]@axes[j]@ms[b]@axes[i]@ms[c]
            per_triple.append(block.reshape(-1))
        blocks.append(np.stack(per_triple))
    image=np.concatenate(blocks,axis=1)
    gram_float=image.conj()@image.T
    assert np.max(np.abs(gram_float.imag))==0
    assert np.max(np.abs(gram_float.real-np.rint(gram_float.real)))==0
    gram=np.rint(gram_float.real).astype(np.int64)
    assert np.array_equal(gram,gram.T)
    assert int(np.max(np.abs(gram)))<2**27
    a_num=sum((kron_label(label) for label in ('ZXII','XIIZ','YXIZ')),
              np.zeros((16,16),dtype=np.complex128))
    moment_b=[]
    moment_ab=[]
    for power in (1,2,3):
        row_b=[np.trace(powers[power]@axis) for axis in axes]
        row_ab=[np.trace(powers[power]@a_num@axis) for axis in axes]
        assert all(z.real.is_integer() and z.imag.is_integer() for z in row_b+row_ab)
        moment_b.append([[int(z.real),int(z.imag)] for z in row_b])
        moment_ab.append([[int(z.real),int(z.imag)] for z in row_ab])
    result={'basis_labels':labels,'monomials':[list(m) for m in monomials],
            'gram':gram.tolist(),'rank_mod_prime_1000003':rank_mod_prime(gram,1000003),
            'max_absolute_gram_entry':int(np.max(np.abs(gram))),
            'moment_B_rows_real_imag':moment_b,
            'moment_A_numerator_B_rows_real_imag':moment_ab}
    (ROOT/'global_lower_bound7_joint_aa_centralizer_gram_result.json').write_text(
        json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in
                      ('basis_labels','monomials','gram','moment_B_rows_real_imag',
                       'moment_A_numerator_B_rows_real_imag')},sort_keys=True))


if __name__=='__main__':main()
