# syndrome distribution from stabilizer expectations, without syndrome projectors
import math

import numpy as np

from .core import I2, X, Y, Z, per_qubit_channels

_MAT = {"I": I2, "X": X, "Y": Y, "Z": Z}
# single-qubit Pauli products with phase: (a, b) -> (phase, c) with a b = phase c
_MUL = {("X", "Y"): (1j, "Z"), ("Y", "X"): (-1j, "Z"), ("Y", "Z"): (1j, "X"), ("Z", "Y"): (-1j, "X"),
        ("Z", "X"): (1j, "Y"), ("X", "Z"): (-1j, "Y")}


def pauli_product(a, b):
    """(phase, string) with a b = phase * string."""
    phase, out = 1.0 + 0j, []
    for x, y in zip(a, b):
        if x == "I":
            out.append(y)
        elif y == "I" or x == y:
            out.append("I" if x == y else x)
        else:
            ph, c = _MUL[(x, y)]
            phase *= ph
            out.append(c)
    return phase, "".join(out)


def group_elements(generators, n):
    """(phase, string) of S^a for a = 0 .. 2^m - 1, bit j of a selecting generator j."""
    m = len(generators)
    out = []
    for a in range(2 ** m):
        phase, s = 1.0 + 0j, "I" * n
        for j in range(m):
            if (a >> j) & 1:
                ph, s = pauli_product(s, generators[j])
                phase *= ph
        out.append((phase, s))
    return out


def apply_pauli_string(psi, s):
    """P |psi> for a Pauli string P, on a state vector with qubit 0 the most significant axis."""
    n = len(s)
    T = psi.reshape([2] * n)
    for q, ch in enumerate(s):
        if ch != "I":
            T = np.moveaxis(np.tensordot(_MAT[ch], T, axes=([1], [q])), 0, q)
    return T.reshape(-1)


def logical_state_vector(code, t=0.0, phi=0.0):
    """cos(t/2) |0_L> + e^{i phi} sin(t/2) |1_L>, from the stabilizer and logical strings."""
    n = code.n
    v = np.zeros(2 ** n, complex)
    v[0] = 1.0
    fixed = list(code.stab_strings) + list(getattr(code, "zl_strs", [code.zl_str]))
    for g in fixed:
        v = (v + apply_pauli_string(v, g)) / 2
    if np.linalg.norm(v) < 1e-9:
        v = np.ones(2 ** n, complex) / math.sqrt(2 ** n)
        for g in fixed:
            v = (v + apply_pauli_string(v, g)) / 2
    v0 = v / np.linalg.norm(v)
    v1 = apply_pauli_string(v0, code.xl_str)
    v1 = v1 / np.linalg.norm(v1)
    out = math.cos(t / 2) * v0 + np.exp(1j * phi) * math.sin(t / 2) * v1
    return out / np.linalg.norm(out)


def adjoint_channel(kraus):
    """E^dagger(P) = sum_k K^dagger P K for P in I, X, Y, Z."""
    return {p: sum(K.conj().T @ _MAT[p] @ K for K in kraus) for p in "IXYZ"}


def expectation(psi, s, adjoint):
    """<psi| (x)_i E_i^dagger(P_i) |psi> for the Pauli string s, with one adjoint per qubit."""
    n = len(s)
    T = psi.reshape([2] * n)
    for q, ch in enumerate(s):
        if ch != "I":
            a = adjoint[q] if isinstance(adjoint, list) else adjoint
            T = np.moveaxis(np.tensordot(a[ch], T, axes=([1], [q])), 0, q)
    return complex(np.vdot(psi, T.reshape(-1)))


def stabilizer_expectations(code, kraus, t=0.0, phi=0.0, psi=None):
    """<S^a>_E(rho) for every element of the stabilizer group, indexed by a."""
    psi = logical_state_vector(code, t, phi) if psi is None else psi
    adj = [adjoint_channel(k) for k in per_qubit_channels(kraus, code.n)]
    elems = group_elements(code.stab_strings, code.n)
    return np.array([phase * expectation(psi, s, adj) for phase, s in elems])


def walsh_hadamard(v):
    """W[s] = sum_a (-1)^(s.a) v[a], in place order."""
    v = np.array(v, dtype=complex)
    h = 1
    while h < len(v):
        for i in range(0, len(v), 2 * h):
            a, b = v[i:i + h].copy(), v[i + h:i + 2 * h].copy()
            v[i:i + h], v[i + h:i + 2 * h] = a + b, a - b
        h *= 2
    return v


def readout_factors(m, q):
    """(1 - 2q_j) per generator in each group element, the effect of flipping syndrome bit j."""
    qs = np.full(m, float(q)) if np.ndim(q) == 0 else np.asarray(q, float)
    assert len(qs) == m, f"expected {m} readout rates, got {len(qs)}"
    f = np.ones(2 ** m)
    for a in range(2 ** m):
        for j in range(m):
            if (a >> j) & 1:
                f[a] *= 1.0 - 2.0 * qs[j]
    return f


def apply_readout(d, q):
    """A syndrome distribution after each bit flips on its own, with probability q."""
    m = int(round(math.log2(len(d))))
    e = walsh_hadamard(np.asarray(d, float)) * readout_factors(m, q)
    out = np.real(walsh_hadamard(e)) / len(e)
    out = np.clip(out, 0.0, None)
    return out / out.sum()


def syndrome_dist(code, kraus, t=0.0, phi=0.0, psi=None, readout=None):
    """P(s) over the 2^m syndromes, bit j of s for generator j, matching Code.syndrome_dist order."""
    e = stabilizer_expectations(code, kraus, t, phi, psi)
    if readout is not None:
        e = e * readout_factors(len(code.stab_strings), readout)
    d = np.real(walsh_hadamard(e)) / len(e)
    d = np.clip(d, 0.0, None)
    return d / d.sum()


def population_leak(code, kraus, readout=None):
    d0 = syndrome_dist(code, kraus, 0.0, 0.0, readout=readout)
    d1 = syndrome_dist(code, kraus, math.pi, 0.0, readout=readout)
    return 0.5 * float(np.abs(d0 - d1).sum()), d0, d1


def rotated_surface_code(d=3):
    """Rotated surface code checks on a d x d grid of data qubits, index d r + c, as (Hx, Hz)."""
    assert d >= 3 and d % 2 == 1, f"odd d of 3 or more, got {d}"
    q = lambda r, c: d * r + c
    zs, xs = [], []
    for r in range(d - 1):
        for c in range(d - 1):
            plaq = (q(r, c), q(r, c + 1), q(r + 1, c), q(r + 1, c + 1))
            (zs if (r + c) % 2 == 0 else xs).append(plaq)
    for r in range(0, d - 1, 2):
        zs.append((q(r, d - 1), q(r + 1, d - 1)))
    for r in range(1, d - 1, 2):
        zs.append((q(r, 0), q(r + 1, 0)))
    for c in range(0, d - 1, 2):
        xs.append((q(0, c), q(0, c + 1)))
    for c in range(1, d - 1, 2):
        xs.append((q(d - 1, c), q(d - 1, c + 1)))
    Hx = np.zeros((len(xs), d * d), np.uint8)
    Hz = np.zeros((len(zs), d * d), np.uint8)
    for i, group in enumerate(xs):
        Hx[i, list(group)] = 1
    for i, group in enumerate(zs):
        Hz[i, list(group)] = 1
    return Hx, Hz


def surface_code_3():
    """The rotated [[9,1,3]] surface code as a strings-only code."""
    from .css import css_strings
    return css_strings(*rotated_surface_code(3), name="rotated surface [[9,1,3]]")
