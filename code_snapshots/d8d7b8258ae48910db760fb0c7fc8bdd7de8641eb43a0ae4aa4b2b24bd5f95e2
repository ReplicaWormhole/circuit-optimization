"""Deterministic bounded beam for one exact reversible orbit encoder."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import json,time,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
def main():
    cfg=json.loads((HERE/'config.json').read_text());assert not (HERE/'result.json').exists()
    assert hashlib.sha256((ROOT/cfg['source']).read_bytes()).hexdigest()==cfg['source_sha256']
    target=tuple(cfg['permutation']);catalog=cfg['catalog'];identity=tuple(range(16));begun=time.monotonic();expansions=0
    def metrics(state):return sum(a!=b for a,b in zip(state,target)),sum((a^b).bit_count() for a,b in zip(state,target))
    def rank(node):
        state,cx,steps,local=node;m,h=metrics(state)
        return (h+cx,m,h,cx,len(steps),local,state,steps)
    beam=[(identity,0,(),0)];best=beam[0];found=None;trace=[];seen={identity:(0,0)};stop='depthcap'
    for depth in range(1,cfg['depth']+1):
        nxt={}
        for state,cx,steps,local in beam:
            for index,g in enumerate(catalog):
                if expansions>=cfg['expansion_cap'] or time.monotonic()-begun>cfg['wall_seconds']:stop='expansioncap' if expansions>=cfg['expansion_cap'] else 'wallcap';break
                expansions+=1;new=tuple(g['permutation'][x] for x in state);cost=cx+g['cx_upper'];seq=steps+(index,);locals_=local+g['local_x'];node=(new,cost,seq,locals_)
                if seen.get(new,(10**9,10**9))<=(cost,len(seq)):continue
                seen[new]=(cost,len(seq))
                if new==target:
                    if found is None or (cost,len(seq),locals_,seq)<(found[1],len(found[2]),found[3],found[2]):found=node
                old=nxt.get(new)
                if old is None or rank(node)<rank(old):nxt[new]=node
                if (metrics(new),cost,len(seq),locals_,seq)<(metrics(best[0]),best[1],len(best[2]),best[3],best[2]):best=node
            if stop!='depthcap':break
        if not nxt:break
        beam=sorted(nxt.values(),key=rank)[:cfg['beam_width']];trace.append({'depth':depth,'expansions':expansions,'best_mismatches':metrics(best[0])[0],'best_hamming':metrics(best[0])[1],'beam':len(beam)})
        if found is not None:best=found;stop='exactfound';break
        if stop!='depthcap':break
    state,cx,indices,locals_=best
    result={'config':cfg,'elapsed_seconds':time.monotonic()-begun,'expansions':expansions,'stop':stop,'trace':trace,'exact_encoder':state==target,'native_steps':[catalog[k] for k in indices],'native_step_count':len(indices),'compiled_cx_upper':cx,'local_x_conjugations':locals_,'actual_permutation':state,'mismatched_rows':[{'input':i,'actual':a,'target':b} for i,(a,b) in enumerate(zip(state,target)) if a!=b],'scope':'Finite heuristic beam. Native reversible encoder only; no circuit optimality or diagonalizer certification.'}
    (HERE/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['config','native_steps']}))
if __name__=='__main__':main()
