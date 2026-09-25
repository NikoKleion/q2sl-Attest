# stabilizer Renyi entropy M2 (Leone, Oliviero and Hamma 2022)
import itertools

import numpy as np

_I = np.eye(2, dtype=complex)
_X = np.array([[0, 1], [1, 0]], complex)
_Y = np.array([[0, -1j], [1j, 0]], complex)
_Z = np.array([[1, 0], [0, -1]], complex)
_H = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)
_T = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], complex)
_P = [_I, _X, _Y, _Z]


def _kron(ms):
    m = np.array([[1]], complex)
    for x in ms:
        m = np.kron(m, x)
    return m


def sre_m2(psi, n):
    """M2 = -log2(2^-n sum_P <P>^4) over the 4^n Pauli strings; zero on stabilizer states."""
    ops = (_kron([_P[i] for i in c]) for c in itertools.product(range(4), repeat=n))
    xi = np.array([abs(np.vdot(psi, O @ psi)) ** 2 for O in ops]) / (2 ** n)
    return float(-np.log2((2 ** n) * np.sum(xi ** 2)))


def reproduce():
    magic = _T @ (_H @ np.array([1, 0], complex))
    m2 = sre_m2(magic, 1)
    return ("magic (stabilizer Renyi entropy of T|+>)", round(m2, 4), f"{round(-np.log2(0.75), 4)} (analytic)",
            abs(m2 + np.log2(0.75)) < 1e-3)


if __name__ == "__main__":
    print(reproduce())
