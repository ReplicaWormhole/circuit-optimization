"""Exact finite trace obstruction for all real B in fixed-AA centralizer.

The run277 Groebner basis implies, for every real spectral-path involution
B in the 12D centralizer, that s=x6-x5=epsilon/sqrt(3), epsilon=+-1,
u=x2=t/sqrt(3), |t|<=1, and x5=-epsilon*t*t/sqrt(3).  Hence
tr(V B)=2*i*epsilon/sqrt(3),
tr(V AB)=(8*t-2*epsilon)/3,
tr(V**2 AB)=(8*t-4*epsilon*t*t)/3.
Check all finite eigenvalue-count assignments to the four (bit0,bit3)
input blocks for compatibility with one diagonal D.
"""

from fractions import Fraction
from itertools import product
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
EIGEN=(1,1j,-1,-1j)
TOTAL=(6,3,4,3)
BLOCKS=((0,0),(0,1),(1,0),(1,1))


def signed_moment(rows,mask,power):
    value=0j
    for (b0,b3),row in zip(BLOCKS,rows):
        sign=(-1)**((b0 if mask&1 else 0)+(b3 if mask&2 else 0))
        value+=sign*sum(n*(lam**power) for n,lam in zip(row,EIGEN))
    assert value.real.is_integer() and value.imag.is_integer()
    return int(value.real),int(value.imag)


def main():
    certificate=json.loads((ROOT/'global_lower_bound7_joint_aa_centralizer_certificate_result.json').read_text())
    assert certificate['remaining_labels']==['IIXZ','IIYZ','ZIIX','ZIIY',
                                             'ZIXX','ZIXY','ZIYX','ZIYY']
    rows=[row for row in product(range(5),repeat=4) if sum(row)==4]
    valid_total=0
    marginal_eligible=0
    first_moment_eligible=0
    compatible=[]
    for r00,r01,r10 in product(rows,repeat=3):
        r11=tuple(n-a-b-c for n,a,b,c in zip(TOTAL,r00,r01,r10))
        if any(v<0 for v in r11) or sum(r11)!=4:
            continue
        blocks=(r00,r01,r10,r11)
        valid_total+=1
        a1=signed_moment(blocks,1,1)
        b1=signed_moment(blocks,2,1)
        if a1[0]!=0 or b1[0]!=0 or a1[1]==0 or b1[1]==0:
            continue
        if signed_moment(blocks,1,2)!=(0,0) or signed_moment(blocks,2,2)!=(0,0):
            continue
        m0,m3=a1[1],b1[1]
        # z0=2/(sqrt(3)*m0), z3=2*epsilon/(sqrt(3)*m3).
        if 4>=3*m0*m0 or 4>=3*m3*m3:
            continue
        marginal_eligible+=1
        ab1=signed_moment(blocks,3,1)
        ab2=signed_moment(blocks,3,2)
        if ab1[1]!=0:
            continue
        first_moment_eligible+=1
        R,Q=ab1[0],ab2[0]
        for epsilon in (-1,1):
            t=Fraction(epsilon,4)+Fraction(epsilon*R,2*m0*m3)
            if abs(t)>1:
                continue
            lhs=8*t-4*epsilon*t*t
            rhs=Fraction(4*epsilon*Q,m0*m3)
            if lhs==rhs:
                compatible.append({'blocks':{''.join(map(str,block)):list(row)
                                             for block,row in zip(BLOCKS,blocks)},
                                   'epsilon':epsilon,'t':str(t),
                                   'm0':m0,'m3':m3})
    result={'valid_total_block_assignments':valid_total,
            'assignments_passing_single_axis_marginals':marginal_eligible,
            'assignments_with_real_product_first_moment':first_moment_eligible,
            'compatible_assignments':len(compatible),
            'example_compatible':compatible[0] if compatible else None,
            'all_compatible':compatible}
    (ROOT/'global_lower_bound7_joint_aa_family_trace_result.json').write_text(
        json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='all_compatible'},sort_keys=True))


if __name__=='__main__':main()
