"""Independent 100-digit simulation of radical ansatz; not exact proof."""
import sys,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
import mpmath as mp
import sympy as sp
HERE=Path(__file__).resolve().parent
SOURCE=ROOT/'collaboration/work/algebraic13_new/canonical_algebraic_ansatz.json'
def main():
 out=HERE/'radical13_high_precision_result.json';assert not out.exists();mp.mp.dps=100;c=json.loads(SOURCE.read_text());u=[[mp.mpc(int(i==j)) for j in range(16)] for i in range(16)];halferror=mp.mpf(0);count=0
 for g in c['gates']:
  if g['gate']=='cx':
   count+=1;mask=1<<(3-g['target']);control=1<<(3-g['control']);u=[u[i^mask] if i&control else u[i] for i in range(16)];continue
  cs=mp.mpf(str(sp.sympify(g['cos_half_exact']).evalf(110)));sn=mp.mpf(str(sp.sympify(g['sin_half_exact']).evalf(110)));halferror=max(halferror,abs(cs*cs+sn*sn-1));axis=g['gate'][1]
  if axis=='x':m=((cs,-mp.j*sn),(-mp.j*sn,cs))
  elif axis=='y':m=((cs,-sn),(sn,cs))
  elif axis=='z':m=((cs-mp.j*sn,0),(0,cs+mp.j*sn))
  else:raise ValueError(g)
  mask=1<<(3-g['qubit'])
  for i in range(16):
   if i&mask:continue
   a=u[i];b=u[i|mask];u[i]=[m[0][0]*x+m[0][1]*y for x,y in zip(a,b)];u[i|mask]=[m[1][0]*x+m[1][1]*y for x,y in zip(a,b)]
 assert count==13
 shift=lambda x:((x&1)<<3)|(x>>1);roots=[mp.mpc(1),mp.j,mp.mpc(-1),-mp.j];diagonal=[sum(u[i][shift(j)]*mp.conj(u[i][j]) for j in range(16)) for i in range(16)];labels=[min(range(4),key=lambda k:abs(z-roots[k])) for z in diagonal];error=max(abs(u[i][shift(j)]-roots[labels[i]]*u[i][j]) for i in range(16) for j in range(16));unity=max(abs(sum(u[i][k]*mp.conj(u[j][k]) for k in range(16))-int(i==j)) for i in range(16) for j in range(16))
 result={'source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'decimal_digits':100,'cnot_count':count,'labels':labels,'uv_du_max_error':str(error),'unitarity_max_error':str(unity),'halfangle_normalization_max_error':str(halferror),'scope':'Highprecision numerical evidence, not exact certificate'};out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
