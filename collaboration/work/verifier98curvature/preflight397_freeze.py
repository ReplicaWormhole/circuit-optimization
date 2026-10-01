"""Hash and source preflight for run397; deliberately no objective/derivative calls."""
import ast, hashlib, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
RUN=ROOT/'experiments/runs/397'
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    cp=RUN/'config.json'; sp=RUN/'diagnostic.py'; cfg=json.loads(cp.read_text())
    if sha(sp)!=cfg['code_sha256'] or cfg['code_sha256']!='2e518feeaa4a57bf6a3497a17696643049bf8123750f81628606478ac9f0ccda': raise SystemExit('diagnostic code hash mismatch')
    if sha(cp)!='f70c638de7ac8fecd9ffc23750ae4b5b1ceecd92236855ecc3f31052f9933d1a': raise SystemExit('diagnostic config hash mismatch')
    for rel,expected in cfg['dependencies'].items():
        if sha(ROOT/rel)!=expected: raise SystemExit(f'dependency mismatch: {rel}')
    parent_path=ROOT/cfg['parent_result']; parent=json.loads(parent_path.read_text())
    if sha(parent_path)!=cfg['dependencies'][cfg['parent_result']]: raise SystemExit('run396 parent result hash mismatch')
    if cfg['seed']!=9261703 or cfg['initialization']!='exact run396 seed9261703 best98 vector; no perturbation/RNG/fits': raise SystemExit('diagnostic point selector differs')
    row=next(r for r in parent['rows'] if r['seed']==cfg['seed'])
    if parent['best_index']!=parent['rows'].index(row) or row['seed']!=9261703: raise SystemExit('seed9261703 is no longer run396 best')
    pre=json.loads((ROOT/'collaboration/work/verifier98curvature/PREFIT_POINT_AUDIT.json').read_text())
    if pre['result_sha256']!=sha(parent_path) or pre['vector_sha256_f64le']!=hashlib.sha256(__import__('numpy').asarray(row['best'],dtype='<f8').tobytes()).hexdigest(): raise SystemExit('point preflight does not bind this run396 best vector')
    source=sp.read_text(); ast.parse(source)
    required=("row = next(r for r in parent['rows'] if r['seed'] == cfg['seed'])",
              "constrained[i, slots == labels[i]] = np.inf",
              "for i in range(16):",
              "gap = min(r['normalized_cost'] for r in alternatives) - label_cost",
              "labels, columns, label_cost = assignment_module.assignment(u, v_np)",
              "d = torch.tensor(labels, dtype=torch.complex128)",
              "residual = m @ v - d[:, None] * m",
              "for name, fun in [('offdiagonal', off), ('fixed_label_branch', branch)]:",
              "gradient = torch.autograd.grad(value, x)[0].detach().numpy()",
              "torch.autograd.functional.hessian(fun, x, vectorize=False)")
    for frag in required:
        if frag not in source: raise SystemExit(f'missing frozen branch/input expression: {frag}')
    if 'scipy.optimize import minimize' in source or 'minimize(' in source: raise SystemExit('unexpected fit/optimizer')
    if source.count('torch.autograd.grad(value, x)')!=1 or source.count('torch.autograd.functional.hessian(')!=1:
        raise SystemExit('AD calls differ from at-most-two-per-objective loop')
    if 'torch.set_num_threads(1)' not in source or 'torch.set_num_interop_threads(1)' not in source or 'threadpool_limits(limits=1)' not in source:
        raise SystemExit('one-thread controls missing')
    if cfg['internal_seconds']!=110 or cfg['external_seconds']!=120 or cfg['threads']!=1 or cfg['max_hessians']!=2 or cfg['max_gradients']!=2 or cfg['max_assignments']!=17:
        raise SystemExit('budget differs from frozen diagnostic assignment')
    if cfg['assignment_gap']!='minimum of16 Hungarian optima banning selectedroot in each row minus basecost; normalized by16; detects distinctroot alternatives, ignoring duplicate-slot permutations': raise SystemExit('assignment gap definition mismatch')
    if cfg['objective_names']!=['offdiagonal','fixed_label_branch']: raise SystemExit('objective branch list mismatch')
    if (RUN/'result.json').exists(): raise SystemExit('run397 result already exists; preflight should precede execution')
    out={'run_id':397,'board_hypothesis_id':93,'diagnostic_sha256':sha(sp),'config_sha256':sha(cp),
         'parent_run396_result_sha256':sha(parent_path),'run396_point_audit_sha256':sha(ROOT/'collaboration/work/verifier98curvature/PREFIT_POINT_AUDIT.json'),
         'point_f64le_sha256':pre['vector_sha256_f64le'],'seed':cfg['seed'],
         'code_hashes_match':True,'all_dependencies_match':True,
         'static_branch_review':{'objective_names':cfg['objective_names'],'point_source':'exact run396 seed9261703 best; no RNG or perturbation',
           'same_frozen_98_coordinate_chart':True,'fixed_label_d_from_pinned_assignment':True,
           'assignment_gap_bans_selected_root_across_all_duplicate_slots_per_row':True,
           'base_plus_16_row_bans':True,'max_explicit_gradients':2,'max_hessians':2,
           'column_eigenvectors_canonicalized_by_largest_magnitude_entry':"for j in range(98)" in source and "vec[:, j] *= -1" in source,
           'negative_trigger':cfg['negative_trigger']},
         'budget':{'internal_seconds':cfg['internal_seconds'],'external_seconds':cfg['external_seconds'],'threads':cfg['threads'],
                   'max_assignments':cfg['max_assignments'],'max_gradients':cfg['max_gradients'],'max_hessians':cfg['max_hessians']},
         'scope':'static script/config/source binding only; no objective evaluations, AD/full-Hessian computations, or gap calculations'}
    outpath=ROOT/'collaboration/work/verifier98curvature/PREFIT397_APPROVAL.json'
    outpath.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
