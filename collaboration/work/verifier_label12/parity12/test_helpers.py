"""Focused kernel derivative, nonadjacent embedding and two-CX identity tests."""
import sys,unittest
from pathlib import Path
import numpy as np
import torch
sys.path.insert(0,str(Path(__file__).resolve().parents[4]))
from native_gate_check import circuit_matrix,embedded_two_qubit,native_matrix
from native_helpers import f_matrix,native_gate,compile_native
class Tests(unittest.TestCase):
 def test_matrix_compilation(self):
  for c,t in [(0,2),(1,3),(0,1),(2,3),(3,0)]:
   for a,b in [(0.,0.),(.27,-.43),(-.38,.12)]:
    g=native_gate(c,t,a,b);native={'n':4,'gates':[g]};compiled=compile_native(native)
    self.assertEqual(sum(x['gate']=='cx' for x in compiled['gates']),2)
    ta,tb=torch.tensor(a,dtype=torch.float64),torch.tensor(b,dtype=torch.float64)
    np.testing.assert_allclose(f_matrix(ta,tb,c,t).detach().numpy(),embedded_two_qubit(native_matrix(g),c,t,4),atol=1e-14)
    np.testing.assert_allclose(circuit_matrix(native),circuit_matrix(compiled),atol=1e-14)
 def test_gradients(self):
  a=torch.tensor(.27,dtype=torch.float64,requires_grad=True);b=torch.tensor(-.43,dtype=torch.float64,requires_grad=True)
  self.assertTrue(torch.autograd.gradcheck(lambda x,y:torch.view_as_real(f_matrix(x,y,0,2)),(a,b),eps=1e-6,atol=1e-6))
if __name__=='__main__':unittest.main()
