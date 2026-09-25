# single-qubit channels; each returns a Kraus list
import math
import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)


def amplitude_damping(gamma):
    # amplitude damping, relaxation toward |0>
    return [np.array([[1, 0], [0, math.sqrt(1 - gamma)]], complex),
            np.array([[0, math.sqrt(gamma)], [0, 0]], complex)]


def amplitude_pumping(gamma):
    # relaxation toward |1>
    return [np.array([[math.sqrt(1 - gamma), 0], [0, 1]], complex),
            np.array([[0, 0], [math.sqrt(gamma), 0]], complex)]


def generalized_amplitude_damping(gamma, p=0.5):
    # finite-temperature amplitude damping; excited-state population 1-p
    a = [np.array([[1, 0], [0, math.sqrt(1 - gamma)]], complex), np.array([[0, math.sqrt(gamma)], [0, 0]], complex)]
    b = [np.array([[math.sqrt(1 - gamma), 0], [0, 1]], complex), np.array([[0, 0], [math.sqrt(gamma), 0]], complex)]
    return [math.sqrt(p) * a[0], math.sqrt(p) * a[1], math.sqrt(1 - p) * b[0], math.sqrt(1 - p) * b[1]]


def channel_from_t1t2(T1, T2, gate_time):
    # build the amplitude-plus-phase-damping Kraus set from T1, T2, and gate time
    g1 = 1.0 - math.exp(-gate_time / T1)
    inv_tphi = max(0.0, 1.0 / T2 - 1.0 / (2.0 * T1))
    lam = 1.0 - math.exp(-2.0 * gate_time * inv_tphi)
    A = [np.array([[1, 0], [0, math.sqrt(1 - g1)]], complex), np.array([[0, math.sqrt(g1)], [0, 0]], complex)]
    P = [np.array([[1, 0], [0, math.sqrt(1 - lam)]], complex), np.array([[0, 0], [0, math.sqrt(lam)]], complex)]
    return [p @ a for p in P for a in A]


def per_qubit_t1t2(t1s, t2s, gate_time):
    """One relaxation channel per qubit, from lists of T1 and T2."""
    return [channel_from_t1t2(a, b, gate_time) for a, b in zip(t1s, t2s)]


def per_qubit_damping(gammas):
    """One amplitude-damping channel per qubit."""
    return [amplitude_damping(g) for g in gammas]


def coherent_diagonal(theta):
    # coherent Z-rotation diag(1, e^{i theta}) on every qubit; single-Kraus unitary channel
    return [np.array([[1, 0], [0, np.exp(1j * theta)]], complex)]


def depolarizing(p):
    # Pauli channel
    s = math.sqrt(p / 3.0)
    return [math.sqrt(1 - p) * I2, s * X, s * Y, s * Z]


def dephasing(p):
    # Z dephasing
    return [math.sqrt(1 - p) * I2, math.sqrt(p) * Z]


LIBRARY = {
    "amplitude_damping": amplitude_damping,
    "generalized_amplitude_damping": generalized_amplitude_damping,
    "coherent_diagonal": coherent_diagonal,
    "depolarizing": depolarizing,
    "dephasing": dephasing,
}
