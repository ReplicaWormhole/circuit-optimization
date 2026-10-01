"""Read-only numerical audit of run398's saved curvature outputs.
No AD/full-Hessian computation, finite differences, optimizer, or additional Hungarian solves.
"""
import hashlib, importlib.util, json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'collaboration/work/verifier98fresh/independent_audit.py'
s=importlib.util.spec_from_file_location('audit398_base',BASE); base=importlib.util.module_from_spec(s); s.loader.exec_module(base)
FIX=ROOT/'collaboration/work/verifier98label/objective_fixtures.py'
s=importlib.util.spec_from_file_location('audit398_fixtures',FIX); fx=importlib.util.module_from_spec(s); s.loader.exec_module(fx)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 run=ROOT/'experiments/runs/398'; cfg=json.loads((run/'config.json').read_text()); respath=run/'result.json'; res=json.loads(respath.read_text())
 if sha(run/'diagnostic.py')!=cfg['code_sha256'] or res.get('config')!=cfg or res.get('run_id')!=398: raise SystemExit('run398 frozen source/config binding failed')
 for rel,h in cfg['dependencies'].items():
  if sha(ROOT/rel)!=h: raise SystemExit(f'dependency changed: {rel}')
 if sha(ROOT/cfg['parent_result'])!=cfg['dependencies'][cfg['parent_result']]: raise SystemExit('run396 parent result hash mismatch')
 parent=json.loads((ROOT/cfg['parent_result']).read_text()); row=next(r for r in parent['rows'] if r['seed']==cfg['seed'])
 x=np.asarray(row['best'],dtype=float); outx=np.asarray(res['point98'],dtype=float)
 if res.get('seed')!=9261703 or not np.array_equal(x,outx) or x.shape!=(98,): raise SystemExit('run398 point differs from frozen run396 seed9261703 best')
 xhash=hashlib.sha256(np.asarray(x,dtype='<f8').tobytes()).hexdigest()
 if res.get('point_sha256_f64le')!=xhash: raise SystemExit('diagnostic point digest mismatch')
 pre=json.loads((ROOT/'collaboration/work/verifier98curvature/PREFIT_POINT_AUDIT.json').read_text())
 if pre['vector_sha256_f64le']!=xhash or pre['result_sha256']!=sha(ROOT/cfg['parent_result']): raise SystemExit('independent point preflight binding mismatch')
 audit396=json.loads((ROOT/'experiments/runs/396/POSTFIT_AUDIT.json').read_text())
 seedrow=next(r for r in audit396['endpoints'] if r['seed']==9261703)
 if res['base_native_sha256']!=seedrow['native_sha256'] or res['base_compiled_sha256']!=seedrow['compiled_sha256']:
  raise SystemExit('run398 base gate list hashes differ from independently audited run396 endpoint')
 npth=run/'base_native.json'; cpth=run/'base_compiled.json'
 if sha(npth)!=res['base_native_sha256'] or sha(cpth)!=res['base_compiled_sha256']: raise SystemExit('run398 base gate-list content hash mismatch')
 u=base.coordinate_matrix(x); un,ng=base.gate_list(npth); uc,cg=base.gate_list(cpth)
 matrix_errors={'coordinate_native':base.phase_error(u,un),'native_compiled':base.phase_error(un,uc)}
 unitary={'native':float(np.max(np.abs(un.conj().T@un-np.eye(16)))), 'compiled':float(np.max(np.abs(uc.conj().T@uc-np.eye(16))))}
 counts={'native_cx':sum(g['gate']=='cx' for g in ng),'native_xx_yy':sum(g['gate']=='xx_yy' for g in ng),'compiled_cx':sum(g['gate']=='cx' for g in cg)}
 if max(matrix_errors.values())>2e-12 or max(unitary.values())>2e-12 or counts!={'native_cx':4,'native_xx_yy':4,'compiled_cx':12}: raise SystemExit('base gate-list matrix/unitarity/count check failed')
 obj=fx.objective_components(u)
 labels=np.asarray([complex(a,b) for a,b in res['labels_real_imag']],dtype=complex)
 expected=np.asarray([complex(a,b) for a,b in obj['labels_by_row_real_imag']],dtype=complex)
 if not np.array_equal(labels,expected) or res['columns']!=obj['assignment_columns']: raise SystemExit('saved base assignment differs from independent preflight assignment')
 if abs(float(res['label_cost'])-obj['raw_assignment_cost_normalized'])>2e-12: raise SystemExit('saved base assigned loss differs')
 if res['base_checker'].get('valid_diagonalizer') is not False: raise SystemExit('base checker flag unexpected')
 # Check the saved constrained alternatives and their costs without solving any new assignment problem.
 v=base.right_shift(); costs=fx.raw_cost_matrix(u,v); slots=fx.LABELS
 alts=res.get('assignment_alternatives',[])
 if len(alts)!=16 or [a.get('banned_row') for a in alts]!=list(range(16)): raise SystemExit('gap record lacks the 16 ordered row bans')
 checked=[]
 for alt in alts:
  i=int(alt['banned_row']); cols=np.asarray(alt['columns'],dtype=int)
  if cols.shape!=(16,) or sorted(cols.tolist())!=list(range(16)): raise SystemExit(f'row ban {i} columns are not a permutation')
  if slots[cols[i]]==labels[i]: raise SystemExit(f'row ban {i} still uses the selected root (duplicate slots included)')
  cost=float(costs[np.arange(16),cols].sum()/16)
  if abs(cost-float(alt['normalized_cost']))>2e-12: raise SystemExit(f'row ban {i} cost mismatch')
  checked.append(cost)
 gap=min(checked)-float(res['label_cost'])
 if res['distinct_root_assignment_gap'] is None or abs(gap-float(res['distinct_root_assignment_gap']))>2e-12: raise SystemExit('recorded gap does not match min of saved alternative costs minus base')
 if len(res['records'])>2: raise SystemExit('more than two derivative records')
 if {r['name'] for r in res['records']}!={'offdiagonal','fixed_label_branch'}: raise SystemExit('diagnostic objective records incomplete or unexpected')
 for rec in res['records']:
  expected_loss=obj['original_offdiagonal_loss'] if rec['name']=='offdiagonal' else obj['label_residual_direct']
  if abs(float(rec['loss'])-expected_loss)>2e-12: raise SystemExit(f"{rec['name']} reported center loss differs from independent preflight")
 hessians=[]; record_reports=[]
 for rec in res['records']:
  grad=np.asarray(rec['gradient'],dtype=float)
  if grad.shape!=(98,) or not np.isfinite(grad).all(): raise SystemExit(f"{rec['name']} gradient malformed")
  if abs(float(np.linalg.norm(grad))-float(rec['gradient_norm']))>2e-10: raise SystemExit(f"{rec['name']} gradient norm mismatch")
  if not rec.get('complete'): raise SystemExit(f"{rec['name']} Hessian incomplete")
  H=np.asarray(rec['hessian'],dtype=float); vals=np.asarray(rec['eigenvalues'],dtype=float); vec=np.asarray(rec['eigenvectors_columns'],dtype=float)
  if H.shape!=(98,98) or vals.shape!=(98,) or vec.shape!=(98,98) or not(np.isfinite(H).all() and np.isfinite(vals).all() and np.isfinite(vec).all()): raise SystemExit(f"{rec['name']} eigensystem shape/nonfinite")
  symerr=float(np.max(np.abs(H-H.T)))
  if abs(symerr-float(rec['hessian_symmetry_max']))>2e-12: raise SystemExit(f"{rec['name']} symmetry diagnostic mismatch")
  Hs=(H+H.T)/2
  ortho=float(np.max(np.abs(vec.T@vec-np.eye(98))))
  residuals=np.linalg.norm(Hs@vec-vec*vals[None,:],axis=0)
  valsref=np.linalg.eigvalsh(Hs)
  evalerr=float(np.max(np.abs(vals-valsref)))
  sign_bad=0
  for j in range(98):
   ix=int(np.argmax(np.abs(vec[:,j])))
   if vec[ix,j]<0: sign_bad+=1
  if ortho>2e-8 or np.max(residuals)>2e-7*max(1.,float(np.linalg.norm(Hs,ord=2))) or evalerr>2e-8*max(1.,float(np.max(np.abs(vals)))) or sign_bad:
   raise SystemExit(f"{rec['name']} independently checked eigensystem failed")
  if abs(float(vals[0])-float(rec['least_eigenvalue']))>2e-10 or bool(vals[0]<cfg['negative_trigger'])!=bool(rec['negative_trigger']): raise SystemExit(f"{rec['name']} least-eigenvalue/trigger mismatch")
  if not np.isfinite(float(rec['least_eigenpair_residual'])) or abs(float(rec['least_eigenpair_residual'])-float(residuals[0]))>2e-7*max(1.,float(np.linalg.norm(Hs,ord=2))): raise SystemExit(f"{rec['name']} least eigenpair residual mismatch")
  hessians.append({'name':rec['name'],'loss':rec['loss'],'gradient_norm':rec['gradient_norm'],
    'hessian_symmetry_max':symerr,'eigenvalue_recompute_max_error':evalerr,
    'eigenvector_column_orthogonality_max_error':ortho,'max_eigenpair_residual':float(np.max(residuals)),
    'least_eigenvalue':float(vals[0]),'least_eigenpair_residual_independent':float(residuals[0]),
    'negative_trigger':bool(vals[0]<cfg['negative_trigger']),'canonical_sign_violations':sign_bad})
  record_reports.append(rec['name'])
 # The saved gap is a row-ban feasible-cost audit; no fresh Hungarian solves are performed here.
 out={'run_id':398,'board_hypothesis_id':96,'diagnostic_sha256':sha(run/'diagnostic.py'),
  'config_sha256':sha(run/'config.json'),'result_sha256':sha(respath),'run396_parent_result_sha256':sha(ROOT/cfg['parent_result']),
  'point_sha256_f64le':xhash,'independent_objective':obj,'gate_list':{'native_sha256':sha(npth),'compiled_sha256':sha(cpth),
   'matrix_errors_up_to_phase':matrix_errors,'unitarity_max_errors':unitary,'counts':counts,'checker_valid':False},
  'assignment_gap':{'base_cost':res['label_cost'],'saved_gap':res['distinct_root_assignment_gap'],
    'checked_16_row_ban_feasible_costs':checked,'computed_gap_from_saved_costs':gap,
    'independent_hungarian_reoptimization_per_ban':'not performed; outside frozen 17-assignment budget'},
  'eigensystems':hessians,'termination':res['termination'],'elapsed_seconds':res['elapsed_seconds'],
  'scope':'Read-only output integrity/eigenpair audit; no AD/full-Hessian calculation, optimizer, finite differences, or new assignment solves'}
 path=ROOT/'collaboration/work/verifier98curvature/POSTFIT398_AUDIT.json'; path.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
if __name__=='__main__':main()
