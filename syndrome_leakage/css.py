# Code from CSS parity-check matrices over GF(2)
import numpy as np
from .core import Code


def _rref(M):
    # reduced row echelon form over GF(2); returns (R, pivot_columns)
    R = (np.asarray(M) % 2).astype(np.uint8).copy()
    rows, cols = R.shape
    piv, r = [], 0
    for c in range(cols):
        pr = next((i for i in range(r, rows) if R[i, c]), None)
        if pr is None:
            continue
        R[[r, pr]] = R[[pr, r]]
        for i in range(rows):
            if i != r and R[i, c]:
                R[i] ^= R[r]
        piv.append(c); r += 1
        if r == rows:
            break
    return R, piv


def _rank(M):
    return len(_rref(M)[1])


def _kernel(M):
    # basis of the right null space of M over GF(2)
    M = (np.asarray(M) % 2).astype(np.uint8)
    if M.size == 0:
        return np.eye(M.shape[1], dtype=np.uint8) if M.shape[1] else np.zeros((0, 0), np.uint8)
    R, piv = _rref(M)
    cols = M.shape[1]
    free = [c for c in range(cols) if c not in piv]
    basis = []
    for f in free:
        v = np.zeros(cols, np.uint8); v[f] = 1
        for i, p in enumerate(piv):
            if R[i, f]:
                v[p] = 1
        basis.append(v)
    return np.array(basis, np.uint8) if basis else np.zeros((0, cols), np.uint8)


def _in_rowspace(v, basis):
    # is v in the GF(2) row space of basis?
    if len(basis) == 0:
        return not np.any(v % 2)
    stack = np.vstack([np.asarray(basis) % 2, (np.asarray(v) % 2)[None, :]])
    return _rank(stack) == _rank(basis)


def _logical_reps(Hx, Hz):
    # logical X reps = ker(Hz) not in rowspace(Hx); logical Z reps = ker(Hx) not in rowspace(Hz)
    kx = [v for v in _kernel(Hz) if not _in_rowspace(v, Hx)]
    kz = [v for v in _kernel(Hx) if not _in_rowspace(v, Hz)]
    return kx, kz


def _symplectic_pairs(kx, kz, k):
    # pair the logical reps so that <z_i, x_j> = delta_ij, by Gram-Schmidt over GF(2)
    zs = [np.asarray(v, np.uint8) % 2 for v in kz]
    xs = [np.asarray(v, np.uint8) % 2 for v in kx]
    pairs = []
    while zs and len(pairs) < k:
        z = zs.pop(0)
        if not np.any(z):
            continue
        j = next((i for i, x in enumerate(xs) if int(np.dot(z, x)) % 2 == 1), None)
        if j is None:
            continue
        x = xs.pop(j)
        pairs.append((z, x))
        zs = [(v ^ z if int(np.dot(v, x)) % 2 else v) for v in zs]
        xs = [(v ^ x if int(np.dot(v, z)) % 2 else v) for v in xs]
    assert len(pairs) == k, f"found {len(pairs)} logical pairs for k={k}"
    for i, (zi, xi) in enumerate(pairs):
        for j, (zj, xj) in enumerate(pairs):
            want = 1 if i == j else 0
            assert int(np.dot(zi, xj)) % 2 == want, "logical pairs are not symplectic"
    return pairs


def _gf2_inverse(M):
    # inverse of a square 0/1 matrix over GF(2), by Gauss-Jordan on [M | I]
    M = (np.asarray(M) % 2).astype(np.uint8)
    k = M.shape[0]
    A = np.hstack([M, np.eye(k, dtype=np.uint8)])
    for c in range(k):
        pr = next((r for r in range(c, k) if A[r, c]), None)
        assert pr is not None, "the supplied logical operators do not pair: lz lx^T is singular over GF(2)"
        A[[c, pr]] = A[[pr, c]]
        for r in range(k):
            if r != c and A[r, c]:
                A[r] ^= A[c]
    return A[:, k:]


def _independent_reps(reps, stab, k):
    # k representatives that stay independent once the stabilizer rows are divided out
    chosen, stack = [], [np.asarray(r, np.uint8) % 2 for r in stab if np.any(r)]
    base = _rank(np.array(stack)) if stack else 0
    for v in reps:
        v = np.asarray(v, np.uint8) % 2
        if _rank(np.array(stack + [v])) > base + len(chosen):
            chosen.append(v)
            stack.append(v)
        if len(chosen) == k:
            break
    assert len(chosen) == k, f"found {len(chosen)} independent logical representatives for k={k}"
    return chosen


def _pairs_from_given(Hx, Hz, lx, lz, k):
    # keep the supplied Z logicals as they are, and re-pair the X logicals inside their own span
    lz = (np.atleast_2d(np.asarray(lz)) % 2).astype(np.uint8)
    assert lz.shape[0] == k, f"expected {k} Z logicals, got {lz.shape[0]}"
    assert not np.any((Hx @ lz.T) % 2), "a supplied Z logical anticommutes with an X check"
    assert _rank(np.vstack([Hz, lz])) == _rank(Hz) + k, "the supplied Z logicals are not independent"
    if lx is None:
        kx, _kz = _logical_reps(Hx, Hz)
        lx = np.array(_independent_reps(kx, Hx, k), np.uint8)
    else:
        lx = (np.atleast_2d(np.asarray(lx)) % 2).astype(np.uint8)
        assert lx.shape[0] == k, f"expected {k} X logicals, got {lx.shape[0]}"
        assert not np.any((Hz @ lx.T) % 2), "a supplied X logical anticommutes with a Z check"
        assert _rank(np.vstack([Hx, lx])) == _rank(Hx) + k, "the supplied X logicals are not independent"
    M = (lz @ lx.T) % 2
    paired = (_gf2_inverse(M).T @ lx) % 2
    assert np.array_equal((lz @ paired.T) % 2, np.eye(k, dtype=np.uint8))
    return [(lz[i], paired[i].astype(np.uint8)) for i in range(k)]


def _to_str(vec, sym):
    return "".join(sym if b else "I" for b in np.asarray(vec) % 2)


def css_from_matrices(Hx, Hz, name="CSS", lx=None, lz=None):
    # build a one-logical-qubit Code from CSS check matrices; requires Hx Hz^T = 0 and k=1
    Hx = (np.asarray(Hx) % 2).astype(np.uint8)
    Hz = (np.asarray(Hz) % 2).astype(np.uint8)
    assert Hx.shape[1] == Hz.shape[1], "Hx and Hz must have the same number of qubits (columns)"
    n = Hx.shape[1]
    assert not np.any((Hx @ Hz.T) % 2), "CSS condition Hx Hz^T = 0 violated (stabilizers do not commute)"
    k = n - _rank(Hx) - _rank(Hz)
    assert k >= 1, f"no logical qubit in this code; got k={k}"
    kx, kz = _logical_reps(Hx, Hz)
    assert kx and kz, "no logical operators found"
    found = _pairs_from_given(Hx, Hz, lx, lz, k) if lz is not None else _symplectic_pairs(kx, kz, k)
    pairs = [(_to_str(z, "Z"), _to_str(x, "X")) for z, x in found]
    stabs = [_to_str(r, "X") for r in Hx if np.any(r)] + [_to_str(r, "Z") for r in Hz if np.any(r)]
    return Code(name, n, stabs, pairs[0][0], pairs[0][1], pairs)


class _StringCode:
    # strings-only code for the analytic path, no density matrix
    def __init__(self, name, n, stab_strings, zl_str, xl_str, logical_pairs=None, active=0):
        self.name = name
        self.n = n
        self.stab_strings = list(stab_strings)
        self.zl_str = zl_str
        self.xl_str = xl_str
        self.logical_pairs = list(logical_pairs) if logical_pairs else [(zl_str, xl_str)]
        self.zl_strs = [z for z, _x in self.logical_pairs]
        self.active = active

    @property
    def k(self):
        return len(self.logical_pairs)

    def with_logical(self, i):
        """The same code with logical qubit i as the one the measures read."""
        z, x = self.logical_pairs[i]
        return _StringCode(f"{self.name} [logical {i}]", self.n, self.stab_strings, z, x,
                           self.logical_pairs, i)


def _min_coset_weight(rep, basis):
    # minimum Hamming weight of rep + span(basis) over GF(2); basis is small (2^rank enumeration)
    rep = np.asarray(rep) % 2
    rows = [np.asarray(b) % 2 for b in basis]
    m = len(rows)
    best = int(np.count_nonzero(rep))
    for mask in range(1, 1 << m):
        v = rep.copy()
        for i in range(m):
            if (mask >> i) & 1:
                v = (v + rows[i]) % 2
        w = int(np.count_nonzero(v))
        if w < best:
            best = w
    return best


def css_matrices(stab_strings):
    """(Hx, Hz) from CSS stabilizer strings."""
    hx, hz = [], []
    for s in stab_strings:
        if set(s) <= {"I", "X"}:
            hx.append([int(ch == "X") for ch in s])
        elif set(s) <= {"I", "Z"}:
            hz.append([int(ch == "Z") for ch in s])
        else:
            raise ValueError(f"not a CSS generator: {s}")
    n = len(stab_strings[0])
    return np.array(hx, np.uint8).reshape(-1, n), np.array(hz, np.uint8).reshape(-1, n)


def code_distance(Hx, Hz):
    # min weight of a nontrivial logical operator (X or Z coset), for verifying a constructed code
    Hx = (np.asarray(Hx) % 2).astype(np.uint8)
    Hz = (np.asarray(Hz) % 2).astype(np.uint8)
    kx, kz = _logical_reps(Hx, Hz)
    dz = min(_min_coset_weight(z, Hz) for z in kz)
    dx = min(_min_coset_weight(x, Hx) for x in kx)
    return min(dx, dz)


def css_strings(Hx, Hz, name="CSS", lx=None, lz=None):
    # same construction and GF(2) verification as css_from_matrices, returning strings only (no Code).
    Hx = (np.asarray(Hx) % 2).astype(np.uint8)
    Hz = (np.asarray(Hz) % 2).astype(np.uint8)
    assert Hx.shape[1] == Hz.shape[1], "Hx and Hz must have the same number of qubits (columns)"
    n = Hx.shape[1]
    assert not np.any((Hx @ Hz.T) % 2), "CSS condition Hx Hz^T = 0 violated"
    k = n - _rank(Hx) - _rank(Hz)
    assert k >= 1, f"no logical qubit in this code; got k={k}"
    kx, kz = _logical_reps(Hx, Hz)
    assert kx and kz, "no logical operators found"
    found = _pairs_from_given(Hx, Hz, lx, lz, k) if lz is not None else _symplectic_pairs(kx, kz, k)
    pairs = [(_to_str(z, "Z"), _to_str(x, "X")) for z, x in found]
    stabs = [_to_str(r, "X") for r in Hx if np.any(r)] + [_to_str(r, "Z") for r in Hz if np.any(r)]
    return _StringCode(name, n, stabs, pairs[0][0], pairs[0][1], pairs)


def hamming_css(r=3):
    # [[2^r-1, 1, 3]] quantum Hamming code from the classical Hamming check matrix, Hx=Hz=H
    ncol = 2 ** r - 1
    H = np.array([[(col >> (r - 1 - bit)) & 1 for bit in range(r)] for col in range(1, ncol + 1)], np.uint8).T
    return css_from_matrices(H, H, name=f"Hamming-CSS[[{ncol},1,3]] (r={r})")
