"""Read-only source roundtrip and final fixed-target residual diagnostics."""
import sys,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT))
import numpy as np
import torch
from check_circuit import right_shift
from delete14_search import u3_to_axis
from search13_adaptive_topology import matrix_for_schedule
from search13_fresh_chain_dynamic_rows import unitary
from search13_fulltail_invariant_gauge_rank import circuit_matrix
HERE=Path(__file__).resolve().parent
cfg=json.loads((HERE/'config.json').read_text());source=json.loads((ROOT/cfg['source']).read_text())
x=np.array([u3_to_axis(g) for g in source['gates'] if g['gate']=='u3']).reshape(14,4,3)
cx=[torch.eye(16,dtype=torch.complex128) if e is None else matrix_for_schedule([e])[0] for e in cfg['optimizer_schedule']]
reconstructed=unitary(torch.tensor(x),cx).detach().numpy();original=circuit_matrix(source['gates'])
overlap=np.vdot(original,reconstructed);phase=overlap/abs(overlap)
rows=[];v=right_shift(4);roots=np.array([1,1j,-1,-1j])
for row in json.loads((HERE/'result.json').read_text())['rows']:
    path=ROOT/row['candidate'];u=circuit_matrix(json.loads(path.read_text())['gates']);d=roots[row['labels']]
    residual=u@v-d[:,None]*u
    rows.append({'candidate':row['candidate'],'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'fixed_target_loss':float(np.sum(abs(residual)**2)/16),'fixed_target_max_entry':float(np.max(abs(residual)))})
result={'source_reconstructed_phase_aligned_max_error':float(np.max(abs(reconstructed-phase*original))),'source_phase':[float(phase.real),float(phase.imag)],'rows':rows,'scope':'Read-only floating diagnostics, no optimization and no exact certificate.'}
assert result['source_reconstructed_phase_aligned_max_error']<1e-12
(HERE/'diagnostics.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
