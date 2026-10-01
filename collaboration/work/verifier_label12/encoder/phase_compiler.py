"""Exact rational-pi closed-Gray parity compiler; not a full UV certificate."""
import json,hashlib
from fractions import Fraction
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
def angle(a):
    a=Fraction(a)
    return f'{a.numerator}*pi/{a.denominator}'
def local(name,q,a=None):
    g={'gate':name,'qubit':q}
    if a is not None:g['theta']=angle(a)
    return g
def phase(wires,values,theta):
    """exp(i pi theta projector) up to global exp(-i pi theta/2^m)."""
    m=len(wires);gates=[];words=[]
    for k,target in enumerate(wires):
        gray=[j^(j>>1) for j in range(1<<k)]
        current=0
        for position,mask in enumerate(gray):
            assert mask==current
            subset=[k]+[j for j in range(k) if mask&(1<<j)]
            sign=(-1)**sum(values[j] for j in subset)
            coefficient=Fraction(theta)*sign/(1<<m)
            gates.append(local('rz',target,-2*coefficient));words.append((subset,coefficient))
            nxt=gray[(position+1)%len(gray)]
            delta=mask^nxt
            if delta:
                assert delta&(delta-1)==0
                j=delta.bit_length()-1
                gates.append({'gate':'cx','control':wires[j],'target':target});current=nxt
        assert current==0
    assert sum(g['gate']=='cx' for g in gates)==(1<<m)-2
    # Exact rational phase polynomial on every conditioned-basis word.
    for x in range(1<<m):
        bits=[(x>>j)&1 for j in range(m)]
        accumulated=sum(c*(-1)**sum(bits[j] for j in subset) for subset,c in words)
        wanted=Fraction(theta)*(bits==values)-Fraction(theta)/(1<<m)
        assert accumulated==wanted
    return gates

def mcx(target,controls,values):
    if not controls:return [local('x',target)]
    if len(controls)==1:
        q=controls[0]
        return ([local('x',q)] if values[0]==0 else [])+[{'gate':'cx','control':q,'target':target}]+([local('x',q)] if values[0]==0 else [])
    return [local('h',target)]+phase(controls+[target],values+[1],1)+[local('h',target)]
def mch(target,controls,values):
    return [local('ry',target,Fraction(-1,4))]+phase(controls+[target],values+[1],1)+[local('ry',target,Fraction(1,4))]
def main():
    source=ROOT/'collaboration/work/numerical_label12/encoder/result.json';data=json.loads(source.read_text());gates=[]
    for g in data['native_steps']:gates+=mcx(g['target'],g['controls'],g['values'])
    p_cost=sum(g['gate']=='cx' for g in gates);assert p_cost==25
    gates += [local('h',2)]+phase([2,3],[1,1],Fraction(1,2))+[local('h',3)]
    gates += mch(3,[0,1],[0,0])+phase([0,1,2,3],[0,0,1,1],Fraction(-1,2))+mch(2,[0,1],[0,0])
    gates += mch(3,[0,1,2],[0,0,1])
    assert sum(g['gate']=='cx' for g in gates)==67
    result={'n':4,'gates':gates,'provenance':{'encoder_source':str(source.relative_to(ROOT)),'encoder_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'cost':{'encoder_cx':25,'orbit_mixing_cx':42,'total_cx':67},'scope':'Exact rational-pi construction and focused phase-polynomial proof; full exact UV certificate not executed here.'}}
    (HERE/'orbit_encoder67.json').write_text(json.dumps(result,indent=2)+'\n')
    (HERE/'PHASE_COMPILER_CHECK.json').write_text(json.dumps({'phase_polynomials':'All basis words checked exactly with Fraction for every emitted projector phase','compiled_cx':67,'encoder_source_sha256':result['provenance']['encoder_sha256'],'whole_UV_certificate':False},indent=2)+'\n')
    print('Exported67CX complete construction; exact phase polynomial tests pass')
if __name__=='__main__':main()
