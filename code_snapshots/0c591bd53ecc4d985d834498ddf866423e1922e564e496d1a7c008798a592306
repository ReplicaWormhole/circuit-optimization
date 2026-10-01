"""Independent int64 Gaussian-integer audit of stored opposite-pair Gram.

No floating point is used in projector-block multiplication or Gram sums.
"""
import json
from pathlib import Path
from itertools import permutations
import numpy as np
ROOT=Path(__file__).resolve().parent

def mul(a,b):
    return (a[0]@b[0]-a[1]@b[1],a[0]@b[1]+a[1]@b[0])

def main():
    d=json.loads((ROOT/'global_lower_bound7_opposite_pair_general_gram_result.json').read_text())
    I=np.eye(2,dtype=np.int64);Z=np.diag([1,-1]);X=np.array([[0,1],[1,0]],dtype=np.int64)
    z=np.zeros((2,2),dtype=np.int64)
    p={'I':(I,z),'X':(X,z),'Z':(Z,z),'Y':(z,np.array([[0,-1],[1,0]],dtype=np.int64))}
    axes=[]
    for label in d['labels']:
        a=(np.ones((1,1),dtype=np.int64),np.zeros((1,1),dtype=np.int64))
        for ch in label:
            b=p[ch];a=(np.kron(a[0],b[0])-np.kron(a[1],b[1]),np.kron(a[0],b[1])+np.kron(a[1],b[0]))
        axes.append(a)
    V=np.zeros((16,16),dtype=np.int64)
    for x in range(16):V[((x&1)<<3)|(x>>1),x]=1
    powers=[np.linalg.matrix_power(V,r) for r in range(4)]
    phases=[(1,0),(0,-1),(-1,0),(0,1)]
    ms=[(sum(phases[k*r%4][0]*powers[r] for r in range(4)),
         sum(phases[k*r%4][1]*powers[r] for r in range(4))) for k in range(4)]
    blocks=[]
    for a,b,c in permutations(range(4),3):
        images=[]
        for i,j in d['monomials']:
            out=mul(mul(mul(mul(ms[a],axes[i]),ms[b]),axes[j]),ms[c])
            if i!=j:
                alt=mul(mul(mul(mul(ms[a],axes[j]),ms[b]),axes[i]),ms[c])
                out=(out[0]+alt[0],out[1]+alt[1])
            images.append(np.concatenate([out[0].ravel(),out[1].ravel()]))
        blocks.append(np.stack(images))
    image=np.concatenate(blocks,axis=1)
    maximum=int(np.max(np.abs(image)))
    assert maximum<=128
    # Every product and sum fits int64: 12288 * 128^2 < 2^31.
    gram=image@image.T
    assert np.array_equal(gram,np.array(d['gram'],dtype=np.int64))
    result={'gaussian_integer_gram_matches':True,'max_block_coordinate':maximum,
            'summation_length':image.shape[1],'bound_on_dot_product':image.shape[1]*maximum**2}
    (ROOT/'global_lower_bound7_opposite_pair_support_audit_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
if __name__=='__main__':main()
