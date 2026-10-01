"""Exact integer basis/ANF verification of a submitted reversible encoder."""
import json,sys,hashlib
from pathlib import Path
EXPECTED=[0,4,7,8,6,2,11,12,5,9,3,13,10,14,15,1]
def apply(g,x):
    target=g['target'];controls=g.get('controls',[])
    assert target in range(4) and all(c['qubit']!=target and c['value'] in (0,1) for c in controls)
    assert len({c['qubit'] for c in controls})==len(controls)
    return x^(1<<(3-target)) if all(((x>>(3-c['qubit']))&1)==c['value'] for c in controls) else x
p=Path(sys.argv[1]);data=json.loads(p.read_text());actual=[]
for x in range(16):
    y=x
    for g in data['gates']:y=apply(g,y)
    actual.append(y)
assert actual==EXPECTED
if 'anf_masks_q0_through_q3' in data:
    for x in range(16):
        value=0
        for q,masks in enumerate(data['anf_masks_q0_through_q3']):
            bit=sum((x&m)==m for m in masks)%2;value|=bit<<(3-q)
        assert value==EXPECTED[x]
costs=[0,1,6,34];cost=sum(costs[len(g.get('controls',[]))] for g in data['gates'])
report={'source':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'exact16basis_mapping':actual,'gates':len(data['gates']),'compiled_cx_upper_bound':cost,'encoder_only_not_diagonalizer':True,'anf_verified': 'anf_masks_q0_through_q3' in data}
Path(sys.argv[2]).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
