"""Independent saved coordinates, seeds and restricted chronology audit."""
import json,hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from check_candidates import matrix,embed,PAULI
ROOT=Path(__file__).resolve().parents[4];folder=ROOT/'collaboration/work/numerical_label12/parity12';r=json.loads((folder/'result.json').read_text());rows=[]
for row in r['rows']:
 seed=row['seed'];np.testing.assert_array_equal(np.random.default_rng(seed).normal(0,.3,size=92),row['initial_params'])
 params=np.array(row['params']);assert params.shape==(92,);locals_=params[:84].reshape(7,4,3);f=params[84:].reshape(4,2);u=np.eye(16,dtype=complex);pairs=[(0,2),(1,3),(0,1),(2,3)];cursor=0;data=json.loads((ROOT/row['native']).read_text())
 for layer in range(7):
  for q in range(4):
   g=data['gates'][cursor];cursor+=1;assert g['gate']=='u3' and g['qubit']==q
   axis=locals_[layer,q];local=expm(-.5j*sum(axis[j]*PAULI[name] for j,name in enumerate(['rx','ry','rz'])));op=np.array([[1]],complex)
   for wire in range(4):op=np.kron(op,local if wire==q else np.eye(2))
   u=op@u
  if layer==0 or layer==3:
   for c,t in ([(0,1),(2,1),(3,1)] if layer==0 else [(1,2)]):
    assert data['gates'][cursor]=={'gate':'cx','control':c,'target':t};cursor+=1;u=matrix({'gates':[{'gate':'cx','control':c,'target':t}]})@u
  elif layer in [1,2,4,5]:
   i={1:0,2:1,4:2,5:3}[layer];g=data['gates'][cursor];cursor+=1;assert g['gate']=='xx_yy' and tuple(g['qubits'])==pairs[i];np.testing.assert_array_equal([g['a'],g['b']],f[i]);a,b=f[i];u=embed(expm(1j*a*np.kron(PAULI['rx'],PAULI['rx'])+1j*b*np.kron(PAULI['ry'],PAULI['ry'])),*pairs[i])@u
 assert cursor==len(data['gates']);v=matrix(data);phase=np.vdot(v,u);phase/=abs(phase);error=float(np.max(abs(u-phase*v)));assert error<1e-12;rows.append({'seed':seed,'exact_seed_initialization_verified':True,'restricted_native_chronology_verified':True,'saved_parameter_matrix_phase_error':error})
output={'checks':rows,'helper_sha256':hashlib.sha256((ROOT/'collaboration/work/verifier_label12/parity12/native_helpers.py').read_bytes()).hexdigest(),'config_sha256':hashlib.sha256((folder/'config.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256((folder/'search.py').read_bytes()).hexdigest()};Path(__file__).with_name('parameter_audit.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))
