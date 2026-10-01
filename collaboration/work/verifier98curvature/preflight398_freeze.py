"""Hash/source-only approval for repaired run398. No matrix evaluations or AD/full-Hessian computations."""
import ast, hashlib, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
RUN=ROOT/'experiments/runs/398'
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 cp=RUN/'config.json'; sp=RUN/'diagnostic.py'; cfg=json.loads(cp.read_text())
 if sha(sp)!=cfg['code_sha256'] or cfg['code_sha256']!='5e67087e57f9b7d3a8915ac41339db7877b22df958a1c60e67896456cac9b086': raise SystemExit('run398 script hash mismatch')
 if sha(cp)!='a5e6ecbf165e548321c5b178edb0c4d2941cfa72dfbce2179fd89b95acaeae13': raise SystemExit('run398 config hash mismatch')
 for rel,h in cfg['dependencies'].items():
  if sha(ROOT/rel)!=h: raise SystemExit(f'dependency changed: {rel}')
 parent=ROOT/cfg['parent_result']; data=json.loads(parent.read_text())
 if sha(parent)!=cfg['dependencies'][cfg['parent_result']]: raise SystemExit('run396 parent hash mismatch')
 row=next(r for r in data['rows'] if r['seed']==cfg['seed'])
 old=json.loads((ROOT/'collaboration/work/verifier98curvature/PREFIT_POINT_AUDIT.json').read_text())
 p_hash=hashlib.sha256(__import__('numpy').asarray(row['best'],dtype='<f8').tobytes()).hexdigest()
 if cfg['seed']!=9261703 or data['best_index']!=data['rows'].index(row) or p_hash!=old['vector_sha256_f64le'] or old['result_sha256']!=sha(parent): raise SystemExit('run398 point does not match independently preflighted run396 global best')
 source=sp.read_text(); ast.parse(source)
 required=("row = next(r for r in parent['rows'] if r['seed'] == cfg['seed'])",
  "labels, columns, label_cost = assignment_module.assignment(u, v_np)",
  "constrained[i, slots == labels[i]] = np.inf",
  "for i in range(16):","gap = min(r['normalized_cost'] for r in alternatives) - label_cost",
  "gradient_tensor = torch.autograd.grad(value, x, create_graph=True)[0]",
  "for j in range(98):","torch.autograd.grad(gradient_tensor[j], x, retain_graph=True)[0]",
  "except BudgetStop:","termination = 'internal_time_limit'",
  "rec['partial_hessian_rows'].append(hrow.tolist())",
  "(HERE / 'result.json').write_text(json.dumps(out, indent=2) + '\\n')")
 for frag in required:
  if frag not in source: raise SystemExit(f'repaired implementation missing expected expression {frag}')
 if 'minimize(' in source or 'scipy.optimize import minimize' in source: raise SystemExit('unexpected optimizer')
 if source.count('torch.autograd.grad(value, x, create_graph=True)')!=1 or source.count('torch.autograd.grad(gradient_tensor[j], x, retain_graph=True)')!=1:
  raise SystemExit('derivative calls exceed frozen first/second-derivative loops')
 # Catch must enclose both assignment and rowwise Hessian loops so the shared handler can serialize partial data.
 catch=source.index('except BudgetStop:')
 trypos=source.index('        try:\n            for i in range(16):')
 looppos=source.index('            for j in range(98):')
 if not (trypos<looppos<catch): raise SystemExit('assignment/Hessian row loops are not within common timeout handler')
 preloop=source[source.index('            for j in range(98):'):source.index('        except BudgetStop:')]
 if 'budget()' not in preloop or 'partial_hessian_rows' not in preloop: raise SystemExit('per-row deadline/partial storage missing')
 if source.index("native = family.serialize(point)") < catch: raise SystemExit('base serialization occurs before timeout catch')
 if cfg['max_hessians']!=2 or cfg['max_first_gradients']!=2 or cfg['max_second_derivative_rows']!=196 or cfg['max_assignments']!=17 or cfg['internal_seconds']!=110 or cfg['external_seconds']!=120 or cfg['threads']!=1:
  raise SystemExit('repaired diagnostic budget differs from frozen assignment')
 if cfg['objective_names']!=['offdiagonal','fixed_label_branch']: raise SystemExit('objective branch names differ')
 if cfg['assignment_gap']!='minimum of16 Hungarian optima banning selectedroot in each row minus basecost; normalized by16; detects distinctroot alternatives, ignoring duplicate-slot permutations': raise SystemExit('assignment-gap definition changed')
 if (RUN/'result.json').exists(): raise SystemExit('run398 already has output; preflight must precede execution')
 out={'run_id':398,'board_hypothesis_id':98,'diagnostic_sha256':sha(sp),'config_sha256':sha(cp),
  'parent_run396_result_sha256':sha(parent),'point_preflight_sha256':sha(ROOT/'collaboration/work/verifier98curvature/PREFIT_POINT_AUDIT.json'),
  'point_f64le_sha256':p_hash,'seed':cfg['seed'],'source_hashes_and_dependencies_match':True,
  'branch_definition_review':{'same_98_vector_and_chart':True,'fixed_D_from_frozen_assignment':True,
   'all_same_root_slots_banned_per_row':True,'base_plus_16_assignment_alternatives':True,
   'offdiagonal_and_fixed_label_branch':True,'hessian_rows_are_manual_second_derivatives':True,
   'eigenvectors_stored_as_columns_and_sign_canonicalized':"vec[:, j] *= -1" in source},
  'timeout_repair_review':{'assignment_and_hessian_loops_inside_shared_catch':True,'budget_checked_before_each_hessian_row':True,
   'partial_rows_retained':True,'base_lists_serialized_after_catch':True},
  'budget':{'max_hessians':2,'max_first_gradients':2,'max_second_derivative_rows':196,
   'max_assignments':17,'internal_seconds':110,'external_seconds':120,'threads':1},
  'scope':'hash/point/source preflight only; no objective, matrix, assignment gap, gradient, or Hessian evaluation'}
 out['preflight_script_sha256']=sha(Path(__file__))
 dst=ROOT/'collaboration/work/verifier98curvature/PREFIT398_APPROVAL.json'; dst.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
