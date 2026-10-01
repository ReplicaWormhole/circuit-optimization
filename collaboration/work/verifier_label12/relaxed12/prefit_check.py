"""Read-only exact coordinate embedding and numerical matrix agreement."""
import json,sys,hashlib,importlib.util
from pathlib import Path
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[4];folder=ROOT/'collaboration/work/numerical_label12/relaxed12'
approval=Path(__file__).with_name('PREFIT_APPROVAL.json')
if approval.exists() or (folder/'result.json').exists():
    raise RuntimeError('Refusing to recreate prefit timing evidence after approval or execution')
import sqlite3
with sqlite3.connect(f'file:{ROOT / "experiments.sqlite3"}?mode=ro',uri=True) as db:
    if any('numerical_label12/relaxed12' in row[0] for row in db.execute('SELECT config_json FROM attempts')):
        raise RuntimeError('Refusing to create prefit evidence after a relaxed12 reservation')
cfg=json.loads((folder/'config.json').read_text())
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'parity12'));from check_candidates import matrix
spec=importlib.util.spec_from_file_location('relaxeddriver',folder/'search.py');driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
source=ROOT/cfg['source'];assert hashlib.sha256(source.read_bytes()).hexdigest()==cfg['source_sha256'];old=json.loads(source.read_text());rows=[]
for entry in cfg['warmstarts']:
 seed=entry['source_seed'];previous=next(r for r in old['rows'] if r['seed']==seed);before=np.array(previous['params']);after=np.array(entry['params116']);expected=np.zeros(116);expected[:12]=before[:12];expected[36:108]=before[12:84];expected[108:]=before[84:];np.testing.assert_array_equal(after,expected)
 native=driver.serialize(after);u=matrix(native);p=ROOT/previous['native'];assert hashlib.sha256(p.read_bytes()).hexdigest()==previous['native_sha256'];v=matrix(json.loads(p.read_text()));error=float(np.max(abs(u-v)));assert error<1e-12
 compiled=driver.helper.compile_native(native);w=matrix(compiled);assert np.max(abs(u-w))<1e-12
 t=driver.matrix(torch.tensor(after,dtype=torch.float64)).detach().numpy();phase=np.vdot(u,t);phase/=abs(phase);torcherror=float(np.max(abs(t-phase*u)));assert torcherror<1e-12
 rows.append({'seed':seed,'source_native_sha256':previous['native_sha256'],'zero_AB_embedding_verified':True,'source_native_matrix_difference':error,'compiled_native_matrix_difference':float(np.max(abs(u-w))),'torch_native_phase_error':torcherror})
for path,digest in cfg['dependencies'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
report={'before_reserve_fit':True,'source_hash_verified':True,'source':cfg['source'],'source_sha256':cfg['source_sha256'],'script_sha256':hashlib.sha256((folder/'search.py').read_bytes()).hexdigest(),'config_sha256':hashlib.sha256((folder/'config.json').read_bytes()).hexdigest(),'checks':rows}
Path(__file__).with_name('PREFIT_APPROVAL.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
