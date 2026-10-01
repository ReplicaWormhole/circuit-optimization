"""Independent explicit-column matrix comparison of native bridge gate lists."""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

def matrix(data):
    assert data["n"] == 4
    result = np.eye(16, dtype=complex)
    count = 0
    for g in data["gates"]:
        if g["gate"] == "cx":
            count += 1
            op = np.zeros((16,16), dtype=complex)
            for x in range(16):
                bits = [(x >> (3-q)) & 1 for q in range(4)]
                bits[g["target"]] ^= bits[g["control"]]
                y = sum(b << (3-q) for q,b in enumerate(bits))
                op[y,x] = 1
        else:
            assert g["gate"] == "u3"
            t,p,l = [g[k] for k in ("theta", "phi", "lam")]
            c,s = np.cos(t/2), np.sin(t/2)
            local = np.array([[c,-np.exp(1j*l)*s],[np.exp(1j*p)*s,np.exp(1j*(p+l))*c]])
            op = np.array([[1]],dtype=complex)
            for q in range(4): op=np.kron(op,local if q==g["qubit"] else np.eye(2))
        result=op@result
    return result,count

if __name__ == "__main__":
    a,b=map(Path,sys.argv[1:3])
    ma,ca=matrix(json.loads(a.read_text())); mb,cb=matrix(json.loads(b.read_text()))
    overlap=np.vdot(mb,ma)
    phase=overlap/abs(overlap) if abs(overlap)>1e-15 else 1
    output={"bridge":str(a),"collapsed":str(b),"bridge_sha256":hashlib.sha256(a.read_bytes()).hexdigest(),"collapsed_sha256":hashlib.sha256(b.read_bytes()).hexdigest(),"cx_counts":[ca,cb],"phase":[float(phase.real),float(phase.imag)],"maximum_difference_up_to_global_phase":float(np.max(np.abs(ma-phase*mb))),"bridge_unitarity_error":float(np.max(np.abs(ma.conj().T@ma-np.eye(16)))),"collapsed_unitarity_error":float(np.max(np.abs(mb.conj().T@mb-np.eye(16))))}
    print(json.dumps(output,indent=2))
    if len(sys.argv)>3: Path(sys.argv[3]).write_text(json.dumps(output,indent=2)+"\n")
    assert ca==13 and cb==12
    assert output["maximum_difference_up_to_global_phase"]<1e-12
