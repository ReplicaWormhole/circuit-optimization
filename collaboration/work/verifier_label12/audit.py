"""Read-only source layout, frozen-target and final target-residual audit."""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
import torch
from verify import matrix, cycle

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from delete14_search import u3_to_axis
from search13_adaptive_topology import matrix_for_schedule
from search13_fresh_chain_dynamic_rows import unitary

folder = ROOT / 'collaboration/work/numerical_label12'
cfg = json.loads((folder/'config.json').read_text())
source = ROOT/cfg['source']
data = json.loads(source.read_text())
gates = data['gates']
cursor = 0
locals_ = []
for layer in range(14):
    for wire in range(4):
        g = gates[cursor]; cursor += 1
        assert g['gate'] == 'u3' and g['qubit'] == wire
        locals_.append(g)
    if layer < 13 and cfg['optimizer_schedule'][layer] is not None:
        c,t = cfg['optimizer_schedule'][layer]
        assert gates[cursor] == {'gate':'cx','control':c,'target':t}
        cursor += 1
assert cursor == len(gates)
u,_ = matrix(data)
axis = np.array([u3_to_axis(g) for g in locals_]).reshape(14,4,3)
cx = [torch.eye(16,dtype=torch.complex128) if e is None else matrix_for_schedule([e])[0]
      for e in cfg['optimizer_schedule']]
optimizer_u = unitary(torch.tensor(axis),cx).detach().numpy()
phase = np.vdot(u,optimizer_u); phase /= abs(phase)
roundtrip = float(np.max(abs(optimizer_u-phase*u)))
assert roundtrip < 1e-12
roots=np.array([1,1j,-1,-1j])
b=u@cycle()@u.conj().T
targets=json.loads((folder/'targets.json').read_text())
labels=np.array(targets['base_labels'])
assert np.array_equal(labels,np.argmin(abs(np.diag(b)[:,None]-roots),axis=1))
assert np.bincount(labels,minlength=4).tolist()==[6,3,4,3]
pairs=sorted([(float(abs(b[r,s])**2+abs(b[s,r])**2),r,s)
              for r in range(16) for s in range(r+1,16) if labels[r]!=labels[s]],
             key=lambda p:(-p[0],p[1],p[2]))[:4]
rows=[]
for index,target in enumerate(targets['targets']):
    assert target['pair']==list(pairs[index][1:])
    wanted=labels.copy();r,s=target['pair'];wanted[r],wanted[s]=wanted[s],wanted[r]
    assert wanted.tolist()==target['labels']
    assert abs(target['coupling_score']-pairs[index][0])<1e-12
    path=folder/f'target{index}.json'
    m,count=matrix(json.loads(path.read_text())); assert count==12
    residual=m@cycle()-roots[wanted,None]*m
    rows.append({'path':str(path.relative_to(ROOT)),
                 'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                 'fixed_target_loss':float(np.sum(abs(residual)**2)/16),
                 'fixed_target_max_residual':float(np.max(abs(residual)))})
result={'source_layout_verified':True,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'optimizer_source_max_error_up_to_phase':roundtrip,
        'frozen_target_selection_verified':True,'targets_sha256':hashlib.sha256((folder/'targets.json').read_bytes()).hexdigest(),
        'config_sha256':hashlib.sha256((folder/'config.json').read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256((folder/'search.py').read_bytes()).hexdigest(),
        'final_fixed_target_checks':rows}
Path(__file__).with_name('independent_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
