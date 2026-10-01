import json
from pathlib import Path
from check_circuit import evaluate

g=[]
def rz(q,a): g.append({'gate':'rz','qubit':q,'theta':a})
def rx(q,a): g.append({'gate':'rx','qubit':q,'theta':a})
def cx(c,t): g.append({'gate':'cx','control':c,'target':t})
def h(q): g.append({'gate':'h','qubit':q})
def cz(a,b): h(b); cx(a,b); h(b)
def h_lift(a,b):
    # Eq. (4.4), P_b(pi/2) then L†, C, Rx_a(pi/4) Rz_b(pi/4), C, L, P_b(pi/2)
    rz(b,'pi/2')
    rx(a,'-pi/2'); rx(b,'-pi/2')
    cx(a,b)
    rx(a,'pi/4'); rz(b,'pi/4')
    cx(a,b)
    rx(a,'pi/2'); rx(b,'pi/2')
    rz(b,'pi/2')
# B4, Appendix A
for j,a in [(1,'pi/8'),(2,'pi/4'),(3,'3*pi/8')]: rz(j,a)
cx(2,0); cx(3,0); rz(0,'pi/8')
cx(1,0); rz(0,'-3*pi/4'); cx(2,0); rz(0,'pi/4')
cx(2,0); cx(3,0); rz(0,'3*pi/8')
cx(1,0); cx(2,0)
# F4, Appendix B, decimation in frequency
cz(2,1)
h_lift(0,2); h_lift(1,3)
cz(2,1)
rz(3,'pi/2')
h_lift(0,1); h_lift(2,3)
result=evaluate({'n':4,'gates':g})
print({k:v for k,v in result.items() if k != 'diagonal_eigenvalues'})
with (Path(__file__).parent / 'baseline_18.json').open('w') as f: json.dump({'n':4,'gates':g},f,indent=2)
