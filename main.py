import math
import sys

from qiskit import QuantumCircuit

from n_toffoli_decomp_utils import *

# This main function implements our compilation scheme Qulin and counts the gates.
n = 8 # for n in range(8, 34, 1):
qr = QuantumRegister(n)
qc = QuantumCircuit(qr)

qc.append(HGate(), [qr[-1]])
ccz_qulin_combined(qc, qr)
qc.append(HGate(), [qr[-1]])

qc = transpile(qc, basis_gates=['u1', 'u2', 'u3', 'cx'], optimization_level=0)  # IBM

cx_cnt = 0
total_cnt = 0

for gate in qc.data:
    if 'cx' == gate.operation.name:
        cx_cnt += 1
    total_cnt += 1

print(cx_cnt)
print(total_cnt)
print(qc.depth())

'''state = Statevector.from_label('+' * n)
state = state.evolve(qc).data
state -= state[0]
state = np.round(state, 5)

print(n, np.count_nonzero(state))'''


