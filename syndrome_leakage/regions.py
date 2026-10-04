# the leak of a part of the record: chosen generators, the smallest sets that leak, one qubit made ideal
import itertools
import math

import numpy as np

from . import zchecks
from .core import per_qubit_channels, tvd
from .expectations import (adjoint_channel, expectation, group_elements, logical_state_vector, population_leak,
                           readout_factors, walsh_hadamard)
from .hardware import marginal_dist, select_checks

IDEAL = [np.eye(2, dtype=complex)]


def region_dist(code, kraus, checks, t=0.0, phi=0.0, readout=None):
    """P over the outcomes of the chosen generators, bit j for the j-th of them, from the products of those generators
    alone. checks is "all", "z", "x" or a list of generator positions; readout is one rate or one per generator."""
    idx = select_checks(code, checks)
    psi = logical_state_vector(code, t, phi)
    adj = [adjoint_channel(k) for k in per_qubit_channels(kraus, code.n)]
    elems = group_elements([code.stab_strings[j] for j in idx], code.n)
    e = np.array([phase * expectation(psi, s, adj) for phase, s in elems])
    if readout is not None:
        q = np.broadcast_to(np.asarray(readout, float), (len(code.stab_strings),))[idx]
        e = e * readout_factors(len(idx), q)
    d = np.clip(np.real(walsh_hadamard(e)) / len(e), 0.0, None)
    return d / d.sum()


def region_leak(code, kraus, checks, readout=None):
    """(leak, d0, d1) for an observer of the chosen generators. It is at most the leak of the whole record, so a
    value of zero does not show that the code is protected."""
    d0 = region_dist(code, kraus, checks, 0.0, 0.0, readout)
    d1 = region_dist(code, kraus, checks, math.pi, 0.0, readout)
    return tvd(d0, d1), d0, d1


def leaking_sets(code, kraus, checks="all", readout=None, tol=1e-12, max_size=None):
    """The least number of generators whose joint record leaks, and every set of that size with its leak, largest
    first: {"size", "sets": [(generator positions, leak)], "whole"}; size None where no set up to max_size leaks."""
    idx = select_checks(code, checks)
    whole, d0, d1 = population_leak(code, kraus, readout)
    for size in range(1, (max_size or len(idx)) + 1):
        found = []
        for c in itertools.combinations(idx, size):
            v = tvd(marginal_dist(d0, c), marginal_dist(d1, c))
            if v > tol:
                found.append((c, v))
        if found:
            return {"size": size, "sets": sorted(found, key=lambda cv: -cv[1]), "whole": whole}
    return {"size": None, "sets": [], "whole": whole}


def qubit_drops(code, kraus, checks="all", readout=None):
    """The leak of the chosen generators, the leak with each qubit's channel made the identity in turn, and the
    difference: {"leak", "without", "drop"}. The drops are not shares: they do not add up to the leak, and at strong
    damping a drop can be negative, the leak rising when that qubit is noiseless."""
    chans = per_qubit_channels(kraus, code.n)
    whole = region_leak(code, chans, checks, readout)[0]
    without = [region_leak(code, chans[:q] + [IDEAL] + chans[q + 1:], checks, readout)[0] for q in range(code.n)]
    return {"leak": whole, "without": without, "drop": [whole - w for w in without]}


def _z_rows(code, checks):
    # the chosen Z generators as rows of Hz, their positions among the code's generators, the qubits they touch
    P = zchecks.css_parts(code)
    idx = select_checks(code, checks)
    pos = {g: i for i, g in enumerate(P["zpos"])}
    if any(g not in pos for g in idx):
        raise ValueError(f"{code.name}: generators {[g for g in idx if g not in pos]} are not of Z type")
    Hz = P["Hz"][[pos[g] for g in idx]].reshape(len(idx), code.n)
    return P, idx, Hz, np.flatnonzero(Hz.any(axis=0))


def z_region_dist(code, p_reset, bit, checks="z", max_terms=2 ** 26):
    """P over the outcomes of the chosen Z generators of one logical state under the reset process, for a CSS code of
    any size. The cost is the number of codewords on the qubits those generators touch times 2^(generators)."""
    P, idx, Hz, U = _z_rows(code, checks)
    r = len(idx)
    basis = zchecks._row_basis(P["Hx"][:, U])
    if 2 ** (len(basis) + r) > max_terms:
        raise ValueError(f"{code.name}: 2^{len(basis)} codewords on {len(U)} qubits times 2^{r} outcomes is too many")
    X = zchecks._row_space(basis) ^ (bit * P["xl"][U])
    T = zchecks._mod2_matmul((np.arange(2 ** r)[:, None] >> np.arange(r)) & 1, Hz[:, U]).astype(np.float64)
    f = 1.0 - 2.0 * np.broadcast_to(np.asarray(p_reset, float), (code.n,))[U]
    d = np.clip(np.real(walsh_hadamard(zchecks._mean_character(T, X, f))) / 2 ** r, 0.0, None)
    return d / d.sum()


def z_region_leak(code, p_reset, checks="z", max_terms=2 ** 26):
    """(leak, d0, d1) for an observer of the chosen Z generators under the reset process: a lower bound on the leak
    of the whole record, at a cost set by the region."""
    d0 = z_region_dist(code, p_reset, 0, checks, max_terms)
    d1 = z_region_dist(code, p_reset, 1, checks, max_terms)
    return tvd(d0, d1), d0, d1


def z_qubit_drops(code, p_reset, checks="z", max_terms=2 ** 26):
    """qubit_drops for the chosen Z generators under the reset process. A qubit no chosen generator touches has no
    effect on their record, so its drop is zero."""
    U = _z_rows(code, checks)[3]
    p = np.array(np.broadcast_to(np.asarray(p_reset, float), (code.n,)))
    whole = z_region_leak(code, p, checks, max_terms)[0]
    without = np.full(code.n, whole)
    for q in U:
        ideal = p.copy()
        ideal[q] = 0.0
        without[q] = z_region_leak(code, ideal, checks, max_terms)[0]
    return {"leak": whole, "without": without.tolist(), "drop": (whole - without).tolist()}


def _walsh(v):
    # expectations.walsh_hadamard on a real vector, without a loop over the blocks
    v, h = np.array(v, float), 1
    while h < len(v):
        w = v.reshape(-1, 2, h)
        w[:, 0, :], w[:, 1, :] = w[:, 0, :] + w[:, 1, :], w[:, 0, :] - w[:, 1, :]
        h *= 2
    return v


def _css_characteristic(code, p_reset, bit, checks, coherence, max_terms):
    # E[a] over the products of the chosen generators; bit None gives E of 0_L minus E of 1_L, summed term by term
    P = zchecks.css_parts(code)
    idx, n = select_checks(code, checks), code.n
    r = len(idx)
    kind = ["X" if "X" in code.stab_strings[g] else "Z" for g in idx]
    xpos, zpos = [j for j, k in enumerate(kind) if k == "X"], [j for j, k in enumerate(kind) if k == "Z"]
    rows = lambda pos, ch: np.array([[c == ch for c in code.stab_strings[idx[j]]] for j in pos],
                                    np.uint8).reshape(len(pos), n)
    X, Z = rows(xpos, "X"), rows(zpos, "Z")
    p = np.array(np.broadcast_to(np.asarray(p_reset, float), (n,)))
    if p.max() >= 1.0:
        raise ValueError("css_region_dist needs every decay probability below 1")
    c = np.sqrt(1.0 - p) if coherence is None else np.array(np.broadcast_to(np.asarray(coherence, float), (n,)))
    U = np.flatnonzero(Z.any(axis=0))
    if len(U) > 63:
        raise ValueError(f"{code.name}: the chosen Z generators touch {len(U)} qubits, more than 63")
    pairs = list(getattr(code, "logical_pairs", [(code.zl_str, code.xl_str)]))
    others = [[ch == "Z" for ch in z] for z, _x in pairs if z != code.zl_str]
    G = np.array(list(P["Hz"]) + others + [P["zl"]], np.uint8).reshape(-1, n)      # the active logical last
    outside = np.setdiff1d(np.arange(n), U)
    null = zchecks._null_space(G[:, outside].T) if len(outside) else np.eye(len(G), dtype=np.uint8)
    if len(null) > 26 or r > 26:
        raise ValueError(f"{code.name}: 2^{len(null)} Z elements or 2^{r} outcomes is too many")
    C = zchecks._row_space(null)
    V = zchecks._mod2_matmul(C, G[:, U]).astype(np.uint64)
    shift = np.arange(len(U), dtype=np.uint64)
    mask = lambda M: (M.astype(np.uint64) << shift).sum(axis=1, dtype=np.uint64)
    vmask = mask(V)
    weight = np.prod(np.where(V == 1, p[U] / (1.0 - p[U]), 1.0), axis=1)
    weight = weight * (2.0 * C[:, -1] if bit is None else 1.0 - 2.0 * bit * C[:, -1])
    bits = lambda k: (np.arange(2 ** k)[:, None] >> np.arange(k)) & 1
    B = zchecks._mod2_matmul(bits(len(xpos)), X)
    A = zchecks._mod2_matmul(bits(len(zpos)), Z)[:, U]
    bmask, amask = mask(B[:, U]), mask(A)
    place = lambda pos: (bits(len(pos)) << np.array(pos, np.int64)).sum(axis=1) if pos else np.zeros(1, np.int64)
    xindex, zindex = place(xpos), place(zpos)
    clear = [np.flatnonzero((vmask & bm) == 0) for bm in bmask]
    if len(amask) * sum(len(k) for k in clear) > max_terms:
        raise ValueError(f"{code.name}: {len(amask) * sum(len(k) for k in clear)} terms, more than max_terms")
    damp = np.prod(np.where(B == 1, c, 1.0), axis=1)
    keep = np.log(1.0 - p[U])
    e = np.zeros(2 ** r)
    for bx, (bm, v) in enumerate(zip(bmask, clear)):
        wmask = amask & ~bm
        w = ((wmask[:, None] >> shift) & np.uint64(1)).astype(float)
        step = max(1, (1 << 22) // max(len(v), 1))
        total = np.concatenate([((vmask[v][None, :] & ~wmask[i:i + step, None]) == 0) @ weight[v]
                                for i in range(0, len(wmask), step)])
        e[xindex[bx] + zindex] = damp[bx] * np.exp(w @ keep) * total
    return e


def css_region_dist(code, p_reset, bit, checks="all", coherence=None, max_terms=2 ** 26):
    """P over the outcomes of the chosen generators of a CSS code, X and Z type together, for one logical state under
    relaxation, at any code size. Each qubit has a decay probability p_reset below 1 and a factor coherence by which X
    and Y shrink, sqrt(1 - p) for amplitude damping when not given. For each product of the chosen X generators the
    cost is 2^(Z generators) times the number of Z stabilizers and Z logicals that lie on the qubits the chosen Z
    generators touch and off that product's support."""
    e = _css_characteristic(code, p_reset, bit, checks, coherence, max_terms)
    d = np.clip(_walsh(e) / len(e), 0.0, None)
    return d / d.sum()


def css_region_leak(code, p_reset, checks="all", coherence=None, max_terms=2 ** 26):
    """(leak, d0, d1) for an observer of the chosen generators of a CSS code under relaxation, X and Z type together:
    a lower bound on the leak of the whole record, at a cost set by the region. The leak is the transform of the
    difference of the two characteristic functions, taken term by term, so it keeps its digits at weak damping."""
    d0 = css_region_dist(code, p_reset, 0, checks, coherence, max_terms)
    d1 = css_region_dist(code, p_reset, 1, checks, coherence, max_terms)
    return css_region_tv(code, p_reset, checks, coherence, max_terms), d0, d1


def css_region_tv(code, p_reset, checks="all", coherence=None, max_terms=2 ** 26):
    """The leak of css_region_leak alone, without the two records."""
    e = _css_characteristic(code, p_reset, None, checks, coherence, max_terms)
    return 0.5 * float(np.abs(_walsh(e)).sum()) / len(e)


def logical_cover(code, max_terms=2 ** 22):
    """The leak order and the qubits that lie on a Z logical of that weight inside the support of a product of Z
    generators, by enumeration: (order, one flag per qubit), or (None, all False) for a protected code. A qubit off
    the cover has a negative drop at weak damping."""
    P = zchecks.css_parts(code)
    n, Hz = code.n, P["Hz"]
    if 2 ** (n + len(Hz)) > max_terms * 2 ** 10:
        raise ValueError(f"{code.name}: {n} qubits and {len(Hz)} Z generators are too many to enumerate")
    z = ((np.arange(2 ** n)[:, None] >> np.arange(n)) & 1).astype(np.uint8)
    ok = ~zchecks._mod2_matmul(z, P["Hx"].T).any(axis=1) & (zchecks._mod2_matmul(z, P["xl"][:, None])[:, 0] == 1)
    for other in P["others"]:
        ok &= zchecks._mod2_matmul(z, other[:, None])[:, 0] == 0
    logicals = z[ok]
    supports = zchecks._row_space(zchecks._row_basis(Hz))
    best, cover = None, np.zeros(n, bool)
    for A in supports:
        inside = logicals[~(logicals & (1 - A)).any(axis=1)]
        for w, v in zip(inside.sum(axis=1), inside):
            if best is None or w < best:
                best, cover = int(w), v.astype(bool).copy()
            elif w == best:
                cover |= v.astype(bool)
    return best, cover


def _logical_rows(P):
    rows = list(P["Hx"]) + [P["xl"]] + list(P["others"])
    return rows, [0.0] * len(P["Hx"]) + [1.0] + [0.0] * len(P["others"])


def region_order(code, checks="z", time_limit=None):
    """The least weight of a Z logical of the active qubit inside the support of a product of the chosen Z generators,
    by integer programming: {"order", "logical", "optimal", "dual_bound"}. order None: their record does not depend
    on the state under relaxation."""
    P, _idx, Hz, _U = _z_rows(code, checks)
    rows, rhs = _logical_rows(P)
    r = zchecks._milp_min_weight(code.n, rows, rhs, extra=Hz, time_limit=time_limit)
    if r["infeasible"]:
        return {"order": None, "logical": None, "optimal": True, "dual_bound": None}
    return {"order": r["weight"], "logical": r["a"], "optimal": r["optimal"] and r["a"] is not None,
            "dual_bound": r["dual_bound"]}


def smallest_leaking_set(code, time_limit=None):
    """The fewest Z generators whose joint record depends on the state under relaxation, by integer programming:
    {"size", "checks", "logical", "optimal", "dual_bound"}. size None with optimal True: the Z record is independent
    of the state; with optimal False: the solver stopped at time_limit before it found a set."""
    P = zchecks.css_parts(code)
    rows, rhs = _logical_rows(P)
    r = zchecks._milp_min_weight(code.n, rows, rhs, extra=P["Hz"], time_limit=time_limit, count="checks")
    if r["infeasible"]:
        return {"size": None, "checks": [], "logical": None, "optimal": True, "dual_bound": None}
    if r["y"] is None:
        return {"size": None, "checks": [], "logical": None, "optimal": False, "dual_bound": r["dual_bound"]}
    checks = [P["zpos"][j] for j in np.flatnonzero(r["y"])]
    return {"size": len(checks), "checks": checks, "logical": r["a"], "optimal": r["optimal"],
            "dual_bound": r["dual_bound"]}
