# the L2 syndrome distance of a CSS code by tensor network contraction
import math
import random

import numpy as np

_X = np.array([[0, 1], [1, 0]], complex)
_Z = np.array([[1, 0], [0, -1]], complex)


def _xz(x, z):
    """X^x Z^z on one qubit."""
    return np.linalg.matrix_power(_X, x) @ np.linalg.matrix_power(_Z, z)


def local_table(kraus):
    """t[u, w, x, z] = Tr(X^u Z^w E^dagger(X^x Z^z)) for one qubit's channel."""
    t = np.zeros((2, 2, 2, 2), complex)
    for x in (0, 1):
        for z in (0, 1):
            A = sum(K.conj().T @ _xz(x, z) @ K for K in kraus)
            for u in (0, 1):
                for w in (0, 1):
                    t[u, w, x, z] = np.trace(_xz(u, w) @ A)
    return t


def _structure(code):
    """Generator strings, logical Z strings and the active logical, checking the code is CSS."""
    stabs = list(code.stab_strings)
    for s in stabs:
        if "Y" in s or ("X" in s and "Z" in s):
            raise ValueError(f"{code.name} is not CSS; generator {s} mixes X and Z")
    zl = list(getattr(code, "zl_strs", [code.zl_str]))
    active = getattr(code, "active", 0)
    if active >= len(zl) or zl[active] != code.zl_str:
        active = zl.index(code.zl_str)
    for z in zl:
        if set(z) - {"I", "Z"}:
            raise ValueError(f"{code.name} has a logical Z that is not Z type: {z}")
    return stabs, zl, active


def _parity(out, variables, flip):
    """T[out, v1..vr] = 1 when out = flip + v1 + ... + vr mod 2."""
    r = len(variables)
    bits = (np.arange(2 ** (r + 1))[:, None] >> np.arange(r, -1, -1)) & 1
    ok = (bits.sum(axis=1) + flip) % 2 == 0
    return ok.astype(float).reshape([2] * (r + 1))


def _locks_x(table):
    """True when the channel never changes a Pauli's X part: t[u, w, x, z] = 0 whenever u != x."""
    scale = max(np.abs(table).max(), 1.0)
    return all(np.abs(table[u, :, x, :]).max() <= 1e-14 * scale for u in (0, 1) for x in (0, 1) if u != x)


def _network(code, kraus):
    """(inputs, arrays, scale): the contraction's value times scale is ||d0 - d1||_2^2."""
    from .css import _rank
    stabs, zl, active = _structure(code)
    n, M, k = code.n, len(stabs), len(zl)
    per_qubit = isinstance(kraus[0], (list, tuple))
    tables = [local_table(kraus[q]) for q in range(n)] if per_qubit else [local_table(kraus)] * n
    merge = all(_locks_x(t) for t in tables)
    others = [l for l in range(k) if l != active]
    xrows = [s for s in stabs if "X" in s]
    MX = len(xrows)
    rX = _rank(np.array([[1 if c == "X" else 0 for c in s] for s in xrows], np.uint8)) if xrows else 0

    inputs, arrays = [], []
    for q in range(n):
        F = np.einsum("uwxz,vyxz->xzuwvy", tables[q], tables[q])
        if merge:
            F = np.stack([F[x, :, x, :, x, :] for x in (0, 1)])          # F[x, z, w, w2]
        jx = [j for j, s in enumerate(stabs) if s[q] == "X"]
        jz = [j for j, s in enumerate(stabs) if s[q] == "Z"]
        lz = [l for l in others if zl[l][q] == "Z"]
        flip = 1 if zl[active][q] == "Z" else 0
        bits = [("x", [("a", j) for j in jx], 0), ("z", [("a", j) for j in jz], 0)]
        if not merge:
            bits.append(("u", [("b", j) for j in jx], 0))
        bits.append(("w", [("b", j) for j in jz] + [("bl", l) for l in lz], flip))
        if not merge:
            bits.append(("u2", [("c", j) for j in jx], 0))
        bits.append(("w2", [("c", j) for j in jz] + [("cl", l) for l in lz], flip))
        legs, index = [], []
        for name, variables, fl in bits:
            if not variables:
                index.append(fl)                              # a constant bit: slice the qubit tensor
            elif len(variables) == 1 and fl == 0:
                index.append(slice(None))
                legs.append(variables[0])                     # the bit is that one variable
            else:
                leg = (name, q)
                index.append(slice(None))
                legs.append(leg)
                inputs.append((leg, *variables))
                arrays.append(_parity(leg, variables, fl))
        inputs.append(tuple(legs))
        arrays.append(F[tuple(index)])
    scale = 2.0 ** (-M) * 4.0 ** (1 - M - k)
    if merge:
        scale *= 4.0 ** (MX - rX)
    return inputs, arrays, scale


_PLANS = {}


def _search():
    """The default plan search: 64 random greedy trials, each refined by subtree reconfiguration, with their parameters
    drawn from a fixed seed."""
    import cotengra as ctg
    return ctg.HyperOptimizer(methods=["random-greedy"], optlib="random", seed=0, max_repeats=64,
                              minimize="combo", reconf_opts={"subtree_size": 4, "maxiter": 100},
                              parallel=False, progbar=False)


def _plan(inputs, arrays, optimize):
    """A contraction tree, kept per network so a sweep over the channel strength searches once."""
    import cotengra as ctg
    key = (tuple(inputs), None if optimize is None else id(optimize))
    if key in _PLANS:
        return _PLANS[key]
    size = {ix: 2 for term in inputs for ix in term}
    state = random.getstate()
    random.seed(0)                # cotengra draws from the global generator wherever no seed reaches it
    try:
        tree = ctg.array_contract_tree(inputs, (), size, optimize=_search() if optimize is None else optimize)
    finally:
        random.setstate(state)
    _PLANS[key] = tree
    return tree


def contraction_width(code, kraus, optimize=None):
    """log2 of the largest intermediate tensor the contraction would build."""
    inputs, arrays, _scale = _network(code, kraus)
    return float(_plan(inputs, arrays, optimize).contraction_width())


def l2_leak(code, kraus, max_width=26, optimize=None):
    """||d0 - d1||_2 between the syndrome distributions of |0_L> and |1_L>, by contraction."""
    inputs, arrays, scale = _network(code, kraus)
    tree = _plan(inputs, arrays, optimize)
    width = float(tree.contraction_width())
    if width > max_width:
        raise MemoryError(f"{code.name}: the contraction needs an intermediate of 2^{width:.1f} entries, "
                          f"past the limit of 2^{max_width}")
    total = tree.contract(arrays)
    norm2 = float(np.real(total)) * scale
    l2 = math.sqrt(max(norm2, 0.0))
    return {"l2": l2, "lower": l2 / 2, "width": width}
