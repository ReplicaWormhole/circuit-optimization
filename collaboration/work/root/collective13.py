"""Collective input SU2 symmetry; numerical simplification only."""
import sys,json,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT))
import numpy as np
from check_circuit import evaluate,one_qubit_matrix,right_shift
from search13_fulltail_invariant_gauge_rank import circuit_matrix
import canonical13
SOURCE=HERE/'labelswap_svd_refine_candidate.json'
OUT=HERE/'collective';
def main():
 OUT.mkdir(exist_ok=True);report=OUT/'collective_report.json';assert not report.exists()
 c=json.loads(SOURCE.read_text());first=c['gates'][:4];assert all(g['gate']=='u3' and g['qubit']==i for i,g in enumerate(first))
 ms=[]
 for g in first:
  m=one_qubit_matrix(g);ms.append(m/np.sqrt(np.linalg.det(m)))
 common=ms[0].conj().T;gates=[]
 for wire,m in enumerate(ms):
  normalized=m@common;a,b,g=canonical13.euler(normalized);reconstructed=canonical13.rot('z',a)@canonical13.rot('y',b)@canonical13.rot('z',g);assert np.max(abs(reconstructed-normalized))<1e-12
  gates.extend([{'gate':'rz','qubit':wire,'theta':float(g)},{'gate':'ry','qubit':wire,'theta':float(b)},{'gate':'rz','qubit':wire,'theta':float(a)}])
 gates.extend(c['gates'][4:]);transformed={'n':4,'gates':gates};p=OUT/'collective_candidate.json';p.write_text(json.dumps(transformed,indent=2)+'\n');q=evaluate(transformed);assert q['cnot_count']==13 and q['valid_diagonalizer']
 m4=np.array([[1]],complex)
 for _ in range(4):m4=np.kron(m4,common)
 v=right_shift(4);commutator=float(np.max(abs(m4@v-v@m4)));assert commutator<1e-12
 u=circuit_matrix(c['gates']);w=circuit_matrix(gates);reference=u@m4;phase=np.trace(reference.conj().T@w)/16;equivalence=float(np.max(abs(w-phase*reference)));assert equivalence<1e-12
 canonical13.SOURCE=p;canonical13.HERE=OUT;canonical13.main()
 result={'source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'commutator_error':commutator,'equivalence_up_global_phase_error':equivalence,'common_matrix':[[[float(z.real),float(z.imag)] for z in row] for row in common],'checker':q,'scope':'Collective input rotation symmetry, not exact angle recognition;13CNOT remains numerical'};report.write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
