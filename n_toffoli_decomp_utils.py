from qiskit.circuit.library.standard_gates import *
import math
import numpy as np
import sys
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.quantum_info import Statevector
from qiskit import transpile
import time
from qiskit.circuit.library import *


# This function implements our compilation scheme Qulin.
def ccz_qulin_combined(qc, operated_qubits):
    n = len(operated_qubits)
    theta = 400001 * math.pi

    # This implements SmallQulin for up to 22 qubits.
    if n <= 22:
        ccz_smallqulin(qc=qc, operated_qubits=operated_qubits, theta=theta)

    # This implements Qulin for 23 or more qubits.
    else:
        ccz_qulin(qc=qc, operated_qubits=operated_qubits, theta=theta)


# This function implements our compilation scheme Qulin for 23 or more qubits.
def ccz_qulin(qc, operated_qubits, theta):
    n = len(operated_qubits)

    # This implements C^{n-1}(V) in Fig.7.
    if n % 2 == 0:

        # This implements U^{n-1}_{+1} in Fig.7.
        large_increment_gate(qc=qc, operated_qubits=operated_qubits[:-1],
                             ancilla_qubits=[operated_qubits[-1]], flag_add=True)

        # This implements U^{n-1}_{V_{n-1}}(-1) in Fig.7.
        theta /= 2
        idx_theta = theta
        for qubit in operated_qubits[1:-1][::-1]:
            c_u1(qc=qc, theta=- idx_theta, operated_qubits=[qubit, operated_qubits[-1]])
            idx_theta /= 2

        # This implements U^{n-1}_{-1} in Fig.7.
        large_increment_gate(qc=qc, operated_qubits=operated_qubits[:-1],
                             ancilla_qubits=[operated_qubits[-1]], flag_add=False)

        # This implements U^{n-1}_{V_{n-1}}(1) in Fig.7.
        idx_theta = theta
        for qubit in operated_qubits[1:-1][::-1]:
            c_u1(qc=qc, theta=idx_theta, operated_qubits=[qubit, operated_qubits[-1]])
            idx_theta /= 2

        # This implements V_{n-1}(1) in Fig.7.
        c_u1(qc=qc, theta=idx_theta * 2, operated_qubits=[operated_qubits[0], operated_qubits[-1]])

    # This implements C^{n-1}(V) in Fig.8.
    else:

        # This implements U^{n-2}_{+1} in Fig.8.
        large_increment_gate(qc=qc, operated_qubits=operated_qubits[:-2],
                             ancilla_qubits=operated_qubits[-2:], flag_add=True)

        # This implements C(U^{n-2}_{V_{n-2}}(-1)) in Fig.8.
        theta /= 2
        idx_theta = theta
        for qubit in operated_qubits[1:-2][::-1]:
            cc_u1(qc=qc, theta=- idx_theta, operated_qubits=[qubit] + operated_qubits[-2:])
            idx_theta /= 2

        # This implements U^{n-2}_{-1} in Fig.8.
        large_increment_gate(qc=qc, operated_qubits=operated_qubits[:-2],
                             ancilla_qubits=operated_qubits[-2:], flag_add=False)

        # This implements C(U^{n-2}_{V_{n-2}}(1)) in Fig.8.
        idx_theta = theta
        for qubit in operated_qubits[1:-2][::-1]:
            cc_u1(qc=qc, theta=idx_theta, operated_qubits=[qubit] + operated_qubits[-2:])
            idx_theta /= 2

        # This implements C(V_{n-2}(1)) in Fig.8
        cc_u1(qc=qc, theta=idx_theta * 2, operated_qubits=[operated_qubits[0]] + operated_qubits[-2:])


# This function implements U^{n-1}_{+1} and U^{n-1}_{-1} in Fig.5 and Fig.6.
def large_increment_gate(qc, operated_qubits, ancilla_qubits, flag_add=True):
    first_half_qubits = operated_qubits[:len(operated_qubits) // 2 + 1]
    second_half_qubits = operated_qubits[len(operated_qubits) // 2 + 1:]

    # This implements U^{n-1}_{-1} in Fig.5.
    if not flag_add:
        for qubit in operated_qubits:
            qc.append(XGate(), [qubit])

    # This implements U^{k_2}_{+1} in Fig.6 by Eq.(6).
    large_plus_1_gate_w_enough_ancilla(qc=qc, operated_qubits=[ancilla_qubits[0]] + second_half_qubits,
                                       ancilla_qubits=first_half_qubits + ancilla_qubits[1:])

    qc.append(XGate(), [ancilla_qubits[0]])

    for qubit in second_half_qubits:
        qc.append(CXGate(), [ancilla_qubits[0], qubit])

    # This implements C^{k_1}(X) in Fig.6 by Eq.(2).
    large_toffoli_w_enough_ancilla(qc=qc, operated_qubits=first_half_qubits + [ancilla_qubits[0]],
                                   ancilla_qubits=second_half_qubits)

    # This implements U^{k_2}_{+1} in Fig.6 by Eq.(6).
    large_plus_1_gate_w_enough_ancilla(qc=qc, operated_qubits=[ancilla_qubits[0]] + second_half_qubits,
                                       ancilla_qubits=first_half_qubits + ancilla_qubits[1:])

    qc.append(XGate(), [ancilla_qubits[0]])

    # This implements C^{k_1}(X) in Fig.6 by Eq.(2).
    large_toffoli_w_enough_ancilla(qc=qc, operated_qubits=first_half_qubits + [ancilla_qubits[0]],
                                   ancilla_qubits=second_half_qubits)

    for qubit in second_half_qubits:
        qc.append(CXGate(), [ancilla_qubits[0], qubit])

    # This implements U^{k_1}_{+1} in Fig.6 by Eq.(6).
    large_plus_1_gate_w_enough_ancilla(qc=qc, operated_qubits=first_half_qubits,
                                       ancilla_qubits=second_half_qubits + ancilla_qubits)

    # This implements U^{n-1}_{-1} in Fig.5.
    if not flag_add:
        for qubit in operated_qubits:
            qc.append(XGate(), [qubit])


# This function implements U^{n/2}_{+1} by Eq.(6) in Gn scheme as shown in Fig.18.
def large_plus_1_gate_w_enough_ancilla(qc, operated_qubits, ancilla_qubits):
    qc.append(XGate(), [ancilla_qubits[0]])
    for qubit in operated_qubits:
        qc.append(XGate().control(1), [ancilla_qubits[0], qubit])
    qc.append(XGate(), [ancilla_qubits[0]])

    # This implements U^{n}_{z+y+x} in Fig.19 and Fig.23.
    for idx in range(len(operated_qubits) - 1):
        ux_gate(qc=qc, operated_qubits=[ancilla_qubits[0], ancilla_qubits[idx + 1], operated_qubits[idx]])
    qc.append(CXGate(), [ancilla_qubits[0], operated_qubits[-1]])
    for idx in range(len(operated_qubits) - 2, -1, -1):
        uz_gate(qc=qc, operated_qubits=[ancilla_qubits[0], ancilla_qubits[idx + 1], operated_qubits[idx]])

    for idx in range(len(operated_qubits) - 1):
        qc.append(XGate(), [ancilla_qubits[idx + 1]])

    # This implements U^{n}_{z+y+x} in Fig.19 and Fig.23.
    for idx in range(len(operated_qubits) - 1):
        ux_gate(qc=qc, operated_qubits=[ancilla_qubits[0], ancilla_qubits[idx + 1], operated_qubits[idx]])
    qc.append(CXGate(), [ancilla_qubits[0], operated_qubits[-1]])
    for idx in range(len(operated_qubits) - 2, -1, -1):
        uz_gate(qc=qc, operated_qubits=[ancilla_qubits[0], ancilla_qubits[idx + 1], operated_qubits[idx]])

    for idx in range(len(operated_qubits) - 1):
        qc.append(XGate(), [ancilla_qubits[idx + 1]])

    qc.append(XGate(), [operated_qubits[-1]])

    qc.append(XGate(), [ancilla_qubits[0]])
    for qubit in operated_qubits:
        qc.append(XGate().control(1), [ancilla_qubits[0], qubit])
    qc.append(XGate(), [ancilla_qubits[0]])


# This function implements U^{3}_{x} in Fig.22.
def ux_gate(qc, operated_qubits):
    qc.append(CXGate(), [operated_qubits[0], operated_qubits[2]])
    qc.append(CXGate(), [operated_qubits[0], operated_qubits[1]])
    qc.append(CCXGate(), [operated_qubits[1], operated_qubits[2], operated_qubits[0]])

# This function implements U^{3}_{z} in Fig.22.
def uz_gate(qc, operated_qubits):
    qc.append(CCXGate(), [operated_qubits[1], operated_qubits[2], operated_qubits[0]])
    qc.append(CXGate(), [operated_qubits[0], operated_qubits[1]])
    qc.append(CXGate(), [operated_qubits[1], operated_qubits[2]])


# This function implements C^{n/2}(X) by Eq.(2) in Iten scheme.
def large_toffoli_w_enough_ancilla(qc, operated_qubits, ancilla_qubits):
    if len(operated_qubits) <= 3:
        qc.append(XGate().control(len(operated_qubits) - 1), operated_qubits)

    else:
        toffoli_qubits = operated_qubits[:]

        for idx in range(len(operated_qubits) - 3):
            toffoli_qubits.insert(2 * (idx + 1), ancilla_qubits[idx])

        qc.append(XGate().control(2), toffoli_qubits[len(toffoli_qubits) - 3:])

        for idx in list(range(len(toffoli_qubits) - 3, 1, -2))[1:]:
            qc.append(RYGate(- math.pi / 4), [toffoli_qubits[idx + 2]])
            qc.append(CXGate(), [toffoli_qubits[idx + 1], toffoli_qubits[idx + 2]])
            qc.append(RYGate(- math.pi / 4), [toffoli_qubits[idx + 2]])
            qc.append(CXGate(), [toffoli_qubits[idx + 0], toffoli_qubits[idx + 2]])

        qc.append(RYGate(- math.pi / 4), [toffoli_qubits[2]])
        qc.append(CXGate(), [toffoli_qubits[1], toffoli_qubits[2]])
        qc.append(RYGate(- math.pi / 4), [toffoli_qubits[2]])
        qc.append(CXGate(), [toffoli_qubits[0], toffoli_qubits[2]])
        qc.append(RYGate(math.pi / 4), [toffoli_qubits[2]])
        qc.append(CXGate(), [toffoli_qubits[1], toffoli_qubits[2]])
        qc.append(RYGate(math.pi / 4), [toffoli_qubits[2]])

        for idx in list(range(0, len(toffoli_qubits) - 4, 2))[1:]:
            qc.append(CXGate(), [toffoli_qubits[idx + 0], toffoli_qubits[idx + 2]])
            qc.append(RYGate(math.pi / 4), [toffoli_qubits[idx + 2]])
            qc.append(CXGate(), [toffoli_qubits[idx + 1], toffoli_qubits[idx + 2]])
            qc.append(RYGate(math.pi / 4), [toffoli_qubits[idx + 2]])

        qc.append(XGate().control(2), toffoli_qubits[len(toffoli_qubits) - 3:])

        for idx in list(range(len(toffoli_qubits) - 3, 1, -2))[1:]:
            qc.append(RYGate(- math.pi / 4), [toffoli_qubits[idx + 2]])
            qc.append(CXGate(), [toffoli_qubits[idx + 1], toffoli_qubits[idx + 2]])
            qc.append(RYGate(- math.pi / 4), [toffoli_qubits[idx + 2]])
            qc.append(CXGate(), [toffoli_qubits[idx + 0], toffoli_qubits[idx + 2]])

        qc.append(RYGate(- math.pi / 4), [toffoli_qubits[2]])
        qc.append(CXGate(), [toffoli_qubits[1], toffoli_qubits[2]])
        qc.append(RYGate(- math.pi / 4), [toffoli_qubits[2]])
        qc.append(CXGate(), [toffoli_qubits[0], toffoli_qubits[2]])
        qc.append(RYGate(math.pi / 4), [toffoli_qubits[2]])
        qc.append(CXGate(), [toffoli_qubits[1], toffoli_qubits[2]])
        qc.append(RYGate(math.pi / 4), [toffoli_qubits[2]])

        for idx in list(range(0, len(toffoli_qubits) - 4, 2))[1:]:
            qc.append(CXGate(), [toffoli_qubits[idx + 0], toffoli_qubits[idx + 2]])
            qc.append(RYGate(math.pi / 4), [toffoli_qubits[idx + 2]])
            qc.append(CXGate(), [toffoli_qubits[idx + 1], toffoli_qubits[idx + 2]])
            qc.append(RYGate(math.pi / 4), [toffoli_qubits[idx + 2]])


# This function implements C(U^{1}_{diag}) by Eq.(3) in Shende and Markov scheme.
def c_u1(qc, theta, operated_qubits):
    qc.append(U1Gate(theta / 2), [operated_qubits[0]])

    qc.append(U1Gate(theta / 2), [operated_qubits[1]])

    qc.append(CXGate(), operated_qubits)

    qc.append(U1Gate(- theta / 2), [operated_qubits[1]])

    qc.append(CXGate(), operated_qubits)


# This function implements C(U^{2}_{diag}) by Eq.(4) in Shende and Markov scheme.
def cc_u1(qc, theta, operated_qubits):
    qc.append(CXGate(), [operated_qubits[0], operated_qubits[2]])

    qc.append(U1Gate(- theta / 4), [operated_qubits[2]])

    qc.append(CXGate(), [operated_qubits[1], operated_qubits[2]])

    qc.append(U1Gate(theta / 4), [operated_qubits[2]])

    qc.append(CXGate(), [operated_qubits[0], operated_qubits[2]])

    qc.append(U1Gate(- theta / 4), [operated_qubits[2]])

    qc.append(CXGate(), [operated_qubits[1], operated_qubits[2]])

    qc.append(U1Gate(theta / 4), [operated_qubits[2]])

    c_u1(qc=qc, theta=theta / 2, operated_qubits=operated_qubits[:-1])


# This function implements our compilation scheme SmallQulin for up to 22 qubits in Fig.8.
def ccz_smallqulin(qc, operated_qubits, theta):

    # This implements U^{n-2}_{+1} in SmallQulin.
    large_increment_gate_smallqulin(qc=qc, operated_qubits=operated_qubits[:-2],
                                    ancilla_qubits=operated_qubits[-2:], flag_add=True)

    theta /= 2
    idx_theta = theta
    for qubit in operated_qubits[1:-2][::-1]:
        cc_u1(qc=qc, theta=- idx_theta, operated_qubits=[qubit] + operated_qubits[-2:])
        idx_theta /= 2

    # This implements U^{n-2}_{-1} in SmallQulin.
    large_increment_gate_smallqulin(qc=qc, operated_qubits=operated_qubits[:-2],
                                    ancilla_qubits=operated_qubits[-2:], flag_add=False)

    idx_theta = theta
    for qubit in operated_qubits[1:-2][::-1]:
        cc_u1(qc=qc, theta=idx_theta, operated_qubits=[qubit] + operated_qubits[-2:])
        idx_theta /= 2

    cc_u1(qc=qc, theta=idx_theta * 2, operated_qubits=[operated_qubits[0]] + operated_qubits[-2:])


# This function implements U^{n-2}_{+1} and U^{n-2}_{-1} in SmallQulin in Fig.6.
def large_increment_gate_smallqulin(qc, operated_qubits, ancilla_qubits, flag_add=True):
    first_half_qubits = operated_qubits[:len(operated_qubits) // 2 + 1]
    second_half_qubits = operated_qubits[len(operated_qubits) // 2 + 1:]

    if not flag_add:
        for qubit in operated_qubits:
            qc.append(XGate(), [qubit])

    large_increment_gate_smallqulin_w_ancilla(qc=qc, operated_qubits=[ancilla_qubits[0]] + second_half_qubits,
                                           ancilla_qubits=first_half_qubits + ancilla_qubits[1:])

    qc.append(XGate(), [ancilla_qubits[0]])

    for qubit in second_half_qubits:
        qc.append(CXGate(), [ancilla_qubits[0], qubit])

    # This implements C^{k_1}(R_z) if 8 <= n <= 19
    if len(first_half_qubits) < 10:
        qc.append(HGate(), [ancilla_qubits[0]])
        large_rz_8_CX(qc=qc, operated_qubits=first_half_qubits + [ancilla_qubits[0]], theta=math.pi / 2)
        qc.append(HGate(), [ancilla_qubits[0]])
    else:
        large_toffoli_w_enough_ancilla(qc=qc, operated_qubits=first_half_qubits + [ancilla_qubits[0]],
                                       ancilla_qubits=second_half_qubits + ancilla_qubits[1:])

    large_increment_gate_smallqulin_w_ancilla(qc=qc, operated_qubits=[ancilla_qubits[0]] + second_half_qubits,
                                           ancilla_qubits=first_half_qubits + ancilla_qubits[1:])

    qc.append(XGate(), [ancilla_qubits[0]])

    # This implements C^{k_1}(R_z) if 8 <= n <= 19
    if len(first_half_qubits) < 10:
        qc.append(HGate(), [ancilla_qubits[0]])
        large_rz_8_CX(qc=qc, operated_qubits=first_half_qubits + [ancilla_qubits[0]], theta=- math.pi / 2)
        qc.append(HGate(), [ancilla_qubits[0]])
    else:
        large_toffoli_w_enough_ancilla(qc=qc, operated_qubits=first_half_qubits + [ancilla_qubits[0]],
                                       ancilla_qubits=second_half_qubits + ancilla_qubits[1:])

    for qubit in second_half_qubits:
        qc.append(CXGate(), [ancilla_qubits[0], qubit])

    large_increment_gate_smallqulin_w_ancilla(qc=qc, operated_qubits=first_half_qubits,
                                           ancilla_qubits=second_half_qubits + ancilla_qubits)

    if not flag_add:
        for qubit in operated_qubits:
            qc.append(XGate(), [qubit])


# This function implements U^{k}_{+1} in SmallQulin scheme.
def large_increment_gate_smallqulin_w_ancilla(qc, operated_qubits, ancilla_qubits):
    n = len(operated_qubits)

    # This implements U^{k}_{+1} in Gn scheme if n >= 11.
    if n >= 11:
        large_plus_1_gate_w_enough_ancilla(qc=qc, operated_qubits=operated_qubits, ancilla_qubits=ancilla_qubits)

    # This implements U^{k}_{+1} in SmallQulin scheme if n < 11.
    else:
        for idx in range(len(operated_qubits), 4, - 1):
            large_toffoli_w_enough_ancilla(qc=qc, operated_qubits=operated_qubits[:idx], ancilla_qubits=ancilla_qubits)

        if n >= 4:
            qc.append(HGate(), [operated_qubits[3]])
            ccc_u1(qc=qc, theta=math.pi, operated_qubits=operated_qubits[:4])
            qc.append(HGate(), [operated_qubits[3]])

        if n >= 3:
            qc.append(HGate(), [operated_qubits[2]])
            cc_u1(qc=qc, theta=math.pi, operated_qubits=operated_qubits[:3])
            qc.append(HGate(), [operated_qubits[2]])

        qc.append(CXGate(), operated_qubits[:2])
        qc.append(XGate(), [operated_qubits[0]])


# This function implements C^3(U^1_{diag}) needed in C^3(X) in SmallQulin scheme.
def ccc_u1(qc, theta, operated_qubits):
    qc.append(U1Gate(theta / 8), [operated_qubits[-1]])
    qc.append(CXGate(), [operated_qubits[1], operated_qubits[3]])

    qc.append(U1Gate(-theta / 8), [operated_qubits[-1]])
    qc.append(CXGate(), [operated_qubits[0], operated_qubits[3]])

    qc.append(U1Gate(theta / 8), [operated_qubits[-1]])
    qc.append(CXGate(), [operated_qubits[1], operated_qubits[3]])

    qc.append(U1Gate(-theta / 8), [operated_qubits[-1]])
    qc.append(CXGate(), [operated_qubits[2], operated_qubits[3]])

    qc.append(U1Gate(theta / 8), [operated_qubits[-1]])
    qc.append(CXGate(), [operated_qubits[1], operated_qubits[3]])

    qc.append(U1Gate(-theta / 8), [operated_qubits[-1]])
    qc.append(CXGate(), [operated_qubits[0], operated_qubits[3]])

    qc.append(U1Gate(theta / 8), [operated_qubits[-1]])
    qc.append(CXGate(), [operated_qubits[1], operated_qubits[3]])

    qc.append(U1Gate(-theta / 8), [operated_qubits[-1]])
    qc.append(CXGate(), [operated_qubits[2], operated_qubits[3]])

    cc_u1(qc=qc, theta=theta / 2, operated_qubits=operated_qubits[:-1])


# This function implements C^{k_1}(R_z) in SmallQulin scheme.
def large_rz_8_CX(qc, operated_qubits, theta):

    # This function implements C^{2}(R_z) in Shende and Markov scheme.
    if len(operated_qubits) == 3:
        qc.append(CXGate(), [operated_qubits[0], operated_qubits[2]])
        qc.append(U1Gate(- theta / 2), [operated_qubits[2]])
        qc.append(CXGate(), [operated_qubits[1], operated_qubits[2]])
        qc.append(U1Gate(theta / 2), [operated_qubits[2]])
        qc.append(CXGate(), [operated_qubits[0], operated_qubits[2]])
        qc.append(U1Gate(- theta / 2), [operated_qubits[2]])
        qc.append(CXGate(), [operated_qubits[1], operated_qubits[2]])
        qc.append(U1Gate(theta / 2), [operated_qubits[2]])

    # This function implements C^{n-1}(R_z) in SmallQulin scheme when n >= 4.
    else:
        control_n = len(operated_qubits) - 1
        qubit_3rd = control_n // 3
        qubit_2nd = (control_n - qubit_3rd) // 2
        qubit_1st = control_n - qubit_3rd - qubit_2nd

        qubit_list_1st = operated_qubits[:qubit_1st] + [operated_qubits[-1]]
        qubit_list_2nd = operated_qubits[qubit_1st: qubit_1st + qubit_2nd] + [operated_qubits[-1]]
        qubit_list_3rd = operated_qubits[qubit_1st + qubit_2nd: -1] + [operated_qubits[-1]]

        qc.append(U1Gate(theta / 4), [operated_qubits[-1]])
        large_toffoli_8_CX(qc=qc, operated_qubits=qubit_list_3rd, theta=math.pi / 2)

        qc.append(U1Gate(-theta / 4), [operated_qubits[-1]])
        large_toffoli_8_CX(qc=qc, operated_qubits=qubit_list_2nd, theta=math.pi / 2)

        qc.append(U1Gate(theta / 4), [operated_qubits[-1]])
        large_toffoli_8_CX(qc=qc, operated_qubits=qubit_list_3rd, theta=math.pi / 2)

        qc.append(U1Gate(-theta / 4), [operated_qubits[-1]])
        large_toffoli_8_CX(qc=qc, operated_qubits=qubit_list_1st, theta=math.pi / 2)

        qc.append(U1Gate(theta / 4), [operated_qubits[-1]])
        large_toffoli_8_CX(qc=qc, operated_qubits=qubit_list_3rd, theta=- math.pi / 2)

        qc.append(U1Gate(-theta / 4), [operated_qubits[-1]])
        large_toffoli_8_CX(qc=qc, operated_qubits=qubit_list_2nd, theta=- math.pi / 2)

        qc.append(U1Gate(theta / 4), [operated_qubits[-1]])
        large_toffoli_8_CX(qc=qc, operated_qubits=qubit_list_3rd, theta=- math.pi / 2)

        qc.append(U1Gate(-theta / 4), [operated_qubits[-1]])
        large_toffoli_8_CX(qc=qc, operated_qubits=qubit_list_1st, theta=- math.pi / 2)


#This function iteratively implements C^{n-1}(R_z) in SmallQulin scheme.
def large_toffoli_8_CX(qc, operated_qubits, theta):
    if len(operated_qubits) > 2:
        qc.append(HGate(), [operated_qubits[-1]])
        large_rz_8_CX(qc=qc, operated_qubits=operated_qubits, theta=theta)
        qc.append(HGate(), [operated_qubits[-1]])
    else:
        qc.append(CXGate(), operated_qubits)

























