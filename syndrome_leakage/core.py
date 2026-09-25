# Code and the exact simulation engine
import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)
_P = {"I": I2, "X": X, "Y": Y, "Z": Z}


def op(s):
    # tensor a Pauli string into a 2^n x 2^n operator
    M = np.array([[1]], complex)
    for ch in s:
        M = np.kron(M, _P[ch])
    return M


def tvd(p, q):
    return 0.5 * float(np.abs(np.asarray(p) - np.asarray(q)).sum())


def _pmul_char(a, b):
    if a == "I":
        return b
    if b == "I":
        return a
    if a == b:
        return "I"
    return {("X", "Y"): "Z", ("Y", "X"): "Z", ("Y", "Z"): "X",
            ("Z", "Y"): "X", ("Z", "X"): "Y", ("X", "Z"): "Y"}[(a, b)]


def pmul(s1, s2):
    # product of two Pauli strings, phase ignored
    return "".join(_pmul_char(a, b) for a, b in zip(s1, s2))


def per_qubit_channels(kraus, n):
    """One Kraus list per qubit, from a single channel or from a list of n channels."""
    first = kraus[0]
    if isinstance(first, (list, tuple)) or np.ndim(first) == 3:
        assert len(kraus) == n, f"expected {n} channels, got {len(kraus)}"
        return [list(k) for k in kraus]
    return [list(kraus)] * n


class Code:
    def __init__(self, name, n, stabilizers, zl, xl, logical_pairs=None, active=0):
        self.name = name
        self.n = n
        self.dim = 2 ** n
        self.stab_strings = list(stabilizers)
        self.zl_str = zl
        self.xl_str = xl
        self.logical_pairs = list(logical_pairs) if logical_pairs else [(zl, xl)]
        self.zl_strs = [z for z, _x in self.logical_pairs]
        self.active = active
        self.GS = [op(s) for s in stabilizers]
        self.ZL = op(zl)
        self.XL = op(xl)
        self.V0, self.V1 = self._codewords()
        self._verify()
        self.PROJ = self._projectors()

    def _codewords(self):
        Pc = np.eye(self.dim, dtype=complex)
        for g in self.GS:
            Pc = Pc @ ((np.eye(self.dim) + g) / 2)
        P0 = Pc
        for z in self.zl_strs:
            P0 = P0 @ ((np.eye(self.dim) + op(z)) / 2)
        ref = np.zeros(self.dim, complex); ref[0] = 1.0
        v0 = P0 @ ref
        if np.linalg.norm(v0) < 1e-9:
            v0 = P0 @ (np.ones(self.dim, complex) / np.sqrt(self.dim))
        v0 = v0 / np.linalg.norm(v0)
        v1 = self.XL @ v0; v1 = v1 / np.linalg.norm(v1)
        return v0, v1

    def _verify(self):
        for g in self.GS:
            assert np.allclose(g @ self.V0, self.V0) and np.allclose(g @ self.V1, self.V1), f"{self.name}: not codewords"
        assert np.allclose(self.ZL @ self.V0, self.V0) and np.allclose(self.ZL @ self.V1, -self.V1), f"{self.name}: Z_L wrong"
        assert abs(np.vdot(self.V0, self.V1)) < 1e-9, f"{self.name}: codewords not orthogonal"

    def _projectors(self):
        projs = {}
        for s in range(2 ** len(self.GS)):
            bits = tuple((s >> j) & 1 for j in range(len(self.GS)))
            P = np.eye(self.dim, dtype=complex)
            for j, g in enumerate(self.GS):
                P = P @ ((np.eye(self.dim) + (-1) ** bits[j] * g) / 2)
            projs[bits] = P
        return projs

    def verify_projectors(self):
        # checks: projectors resolve the identity, are idempotent, commute with both logicals
        Isum = np.zeros((self.dim, self.dim), complex)
        for P in self.PROJ.values():
            Isum += P
            assert np.allclose(P @ P, P), f"{self.name}: projector not idempotent"
            assert np.allclose(P @ self.ZL, self.ZL @ P), f"{self.name}: [Pi_s, Z_L] != 0"
            assert np.allclose(P @ self.XL, self.XL @ P), f"{self.name}: [Pi_s, X_L] != 0 (projector carries a logical)"
        assert np.allclose(Isum, np.eye(self.dim)), f"{self.name}: projectors do not sum to I"
        return True

    @property
    def k(self):
        return len(self.logical_pairs)

    def with_logical(self, i):
        """The same code with logical qubit i as the one the measures read."""
        z, x = self.logical_pairs[i]
        return Code(f"{self.name} [logical {i}]", self.n, self.stab_strings, z, x, self.logical_pairs, i)

    def logical_state(self, t, phi):
        v = np.cos(t / 2) * self.V0 + np.exp(1j * phi) * np.sin(t / 2) * self.V1
        v = v / np.linalg.norm(v)
        return np.outer(v, v.conj())

    def _apply_1q(self, rho, kraus, q):
        n = self.n
        T = rho.reshape([2] * n + [2] * n)
        out = np.zeros_like(T)
        for K in kraus:
            A = np.moveaxis(np.tensordot(K, T, axes=([1], [q])), 0, q)
            A = np.moveaxis(np.tensordot(K.conj(), A, axes=([1], [n + q])), 0, n + q)
            out += A
        return out.reshape(self.dim, self.dim)

    def apply(self, rho, single_qubit_kraus):
        # one single-qubit channel on every physical qubit, or one channel per qubit
        chans = per_qubit_channels(single_qubit_kraus, self.n)
        for K in chans:
            assert np.allclose(sum(k.conj().T @ k for k in K), np.eye(2)), "channel not trace-preserving"
        r = rho
        for q in range(self.n):
            r = self._apply_1q(r, chans[q], q)
        return r

    def syndrome_dist(self, rho):
        d = np.array([float(np.real(np.trace(P @ rho))) for P in self.PROJ.values()])
        d = np.clip(d, 0.0, None)
        return d / d.sum()

    def stabilizer_group(self):
        # all Pauli strings in <stabilizers, Z_L>, phase ignored
        gens = self.stab_strings + [self.zl_str]
        group = {"I" * self.n}
        for g in gens:
            group |= {pmul(x, g) for x in group}
        changed = True
        while changed:
            changed = False
            for a in list(group):
                for g in gens:
                    p = pmul(a, g)
                    if p not in group:
                        group.add(p); changed = True
        return group

    def weight_profile(self, v):
        probs = np.abs(v) ** 2
        wt = np.array([bin(x).count("1") for x in range(self.dim)])
        ni = np.zeros(self.n)
        for x in range(self.dim):
            if probs[x] > 1e-15:
                for qb in range(self.n):
                    if (x >> (self.n - 1 - qb)) & 1:
                        ni[qb] += probs[x]
        return float((probs * wt).sum()), ni
