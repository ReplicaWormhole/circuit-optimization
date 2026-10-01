"""Independent prefix composition checks for the saved animation matrices."""
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from check_circuit import one_qubit_matrix


def main():
    output = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / 'video'
    arrays = np.load(output / 'matrices.npz')
    units, shifted = arrays['U'], arrays['conjugated_shift']
    candidate = json.loads((ROOT / 'experiments/runs/411/candidate.json').read_text())
    assert units.shape == shifted.shape == (31, 16, 16)
    u = np.eye(16, dtype=complex)
    np.testing.assert_allclose(units[0], u, atol=1e-12, rtol=0)
    # Build the right-shift permutation from bit strings, not the shared helper.
    v = np.zeros((16, 16), dtype=complex)
    for b in range(16):
        bits = f'{b:04b}'
        v[int(bits[-1] + bits[:-1], 2), b] = 1
    for k, gate in enumerate(candidate['gates'], 1):
        if gate['gate'] == 'cx':
            updated = np.zeros_like(u)
            for b in range(16):
                bits = list(f'{b:04b}')
                if bits[gate['control']] == '1':
                    bits[gate['target']] = str(1 - int(bits[gate['target']]))
                updated[int(''.join(bits), 2)] = u[b]
            u = updated
        else:
            tensor = u.reshape(2, 2, 2, 2, 16)
            tensor = np.tensordot(one_qubit_matrix(gate), tensor, axes=(1, gate['qubit']))
            u = np.moveaxis(tensor, 0, gate['qubit']).reshape(16, 16)
        np.testing.assert_allclose(units[k], u, atol=1e-12, rtol=0)
    for k, u in enumerate(units):
        np.testing.assert_allclose(u @ u.conj().T, np.eye(16), atol=1e-12, rtol=0)
        np.testing.assert_allclose(shifted[k], u @ v @ u.conj().T, atol=1e-12, rtol=0)
    acceptance = json.loads((ROOT / 'collaboration/EXACT11_ACCEPTANCE.json').read_text())
    expected = np.array([1, 1j, -1, -1j])[acceptance['output_labels']]
    np.testing.assert_allclose(shifted[-1], np.diag(expected), atol=1e-12, rtol=0)
    assert len(list((output / 'frames').glob('step_*.png'))) == 31
    print('PASS: all prefixes, all unitarities, all conjugations, accepted final diagonal, 31 frames.')
    print('Final diagonal error:', np.max(np.abs(shifted[-1] - np.diag(expected))))


if __name__ == '__main__':
    main()
