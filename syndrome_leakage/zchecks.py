# the Z checks of a CSS code under relaxation, as a classical process
import numpy as np

from . import estimate


def css_parts(code):
    """Hx and Hz (nonzero rows), the positions of the Z generators among the code's generators, and the X and Z
    supports of the active logical and the X supports of the other logical qubits."""
    stabs = list(code.stab_strings)
    for s in stabs:
        if "Y" in s or ("X" in s and "Z" in s):
            raise ValueError(f"{code.name} is not CSS; generator {s} mixes X and Z")
    vec = lambda s, c: np.array([1 if ch == c else 0 for ch in s], np.uint8)
    xrows = [vec(s, "X") for s in stabs if "X" in s]
    zpos = [j for j, s in enumerate(stabs) if "Z" in s]
    pairs = list(getattr(code, "logical_pairs", [(code.zl_str, code.xl_str)]))
    for z, x in pairs:
        if set(z) - {"I", "Z"} or set(x) - {"I", "X"}:
            raise ValueError(f"{code.name}: a logical is not of pure Z or X type")
    active = [z for z, _x in pairs].index(code.zl_str)
    n = code.n
    return {"Hx": np.array(xrows, np.uint8).reshape(-1, n), "Hz": np.array([vec(stabs[j], "Z") for j in zpos],
            np.uint8).reshape(-1, n), "zpos": zpos, "xl": vec(pairs[active][1], "X"),
            "zl": vec(pairs[active][0], "Z"),
            "others": [vec(x, "X") for i, (_z, x) in enumerate(pairs) if i != active]}


def _mod2_matmul(a, b):
    # exact for sums below 2^24
    return (a.astype(np.float32) @ b.astype(np.float32)).astype(np.int64) % 2


def _row_basis(M):
    """A reduced basis of the GF(2) row space of M."""
    B = np.atleast_2d(np.array(M, np.uint8) % 2)
    r = 0
    for c in range(B.shape[1]):
        piv = next((i for i in range(r, len(B)) if B[i, c]), None)
        if piv is None:
            continue
        B[[r, piv]] = B[[piv, r]]
        for i in range(len(B)):
            if i != r and B[i, c]:
                B[i] ^= B[r]
        r += 1
    return B[:r]


def _null_space(M):
    """A basis of {c : M c = 0} over GF(2), one vector per row."""
    M = np.atleast_2d(np.array(M, np.uint8) % 2)
    rows, cols = M.shape
    pivots, r = [], 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        pivots.append(c)
        r += 1
    out = []
    for f in (c for c in range(cols) if c not in pivots):
        v = np.zeros(cols, np.uint8)
        v[f] = 1
        for i, c in enumerate(pivots):
            v[c] = M[i, f]
        out.append(v)
    return np.array(out, np.uint8).reshape(len(out), cols)


def _row_space(basis):
    """Every GF(2) combination of the rows of a basis, one per row."""
    r = len(basis)
    return _mod2_matmul((np.arange(2 ** r)[:, None] >> np.arange(r)) & 1, basis).astype(np.uint8)


def _mean_character(T, X, f, chunk=1 << 22):
    """For each row a of T, the mean over the codewords x in X of the product over q in x of f_q where a has q."""
    logf, neg, zero = T * np.log(np.abs(np.where(f == 0, 1.0, f))), T * (f < 0), T * (f == 0)
    total, step = np.zeros(len(T)), max(1, chunk // len(T))
    for i in range(0, len(X), step):
        x = X[i:i + step].T.astype(np.float64)
        chi = np.exp(logf @ x)
        chi *= 1.0 - 2.0 * ((neg @ x) % 2)
        chi[(zero @ x) > 0] = 0.0
        total += chi.sum(axis=1)
    return total / len(X)


def sample_syndromes(code, p_reset, shots, seed=0, chunk=65536, errors=False):
    """(s0, s1): Z syndromes of |0_L> and |1_L>, one row per shot, one column per Z generator in the code's order."""
    P = css_parts(code)
    Hx, Hz, xl, n = P["Hx"], P["Hz"], P["xl"], code.n
    p = np.broadcast_to(np.asarray(p_reset, float), (n,))
    out, errs = [], []
    for bit in (0, 1):
        rng = np.random.default_rng([*np.atleast_1d(seed), bit])
        parts, eparts = [], []
        for start in range(0, shots, chunk):
            m = min(chunk, shots - start)
            coef = rng.integers(0, 2, (m, len(Hx)), dtype=np.uint8)
            x = _mod2_matmul(coef, Hx).astype(np.uint8) ^ (bit * xl)
            E = x & (rng.random((m, n)) < p)
            parts.append(_mod2_matmul(E, Hz.T).astype(np.uint8))
            if errors:
                eparts.append(E.astype(np.uint8))
        out.append(np.vstack(parts))
        errs.append(np.vstack(eparts) if errors else None)
    return (out[0], out[1], errs[0], errs[1]) if errors else (out[0], out[1])


def exact_distribution(code, p_reset, bit, statistic=None, max_terms=2 ** 22):
    """{value: probability} of a statistic of the Z syndrome, by enumerating every stabilizer codeword and every reset
    set."""
    P = css_parts(code)
    Hx, Hz, xl, n = P["Hx"], P["Hz"], P["xl"], code.n
    span = np.unique(_mod2_matmul(((np.arange(2 ** len(Hx))[:, None] >> np.arange(len(Hx))) & 1), Hx), axis=0)
    if len(span) * 2 ** n > max_terms:
        raise ValueError(f"{code.name}: {len(span)} codewords times 2^{n} reset sets is too many to enumerate")
    R = ((np.arange(2 ** n)[:, None] >> np.arange(n)) & 1).astype(np.uint8)
    p = np.broadcast_to(np.asarray(p_reset, float), (n,))
    wR = np.prod(np.where(R == 1, p, 1 - p), axis=1) / len(span)
    out = {}
    for c in span:
        x = (c.astype(np.uint8) ^ (bit * xl)).astype(np.uint8)
        s = _mod2_matmul(R & x, Hz.T).astype(np.uint8)
        vals = s if statistic is None else statistic(s)
        keys = [tuple(v) for v in vals] if statistic is None else [int(v) for v in vals]
        for key, w in zip(keys, wR):
            out[key] = out.get(key, 0.0) + w
    return out


def exact_z_distribution(code, p_reset, bit, max_terms=2 ** 26):
    """P(s) over the Z syndromes of one logical state under the reset process, bit j of s for Z generator j, from the
    characteristic function E[(-1)^(a.s)] = mean over codewords x of the product over q in x of (1 - 2 p_q [a.Hz_q])."""
    from .expectations import walsh_hadamard
    P = css_parts(code)
    Hx, Hz, xl, n = P["Hx"], P["Hz"], P["xl"], code.n
    mz = len(Hz)
    span = np.unique(_mod2_matmul(((np.arange(2 ** len(Hx))[:, None] >> np.arange(len(Hx))) & 1), Hx), axis=0)
    if len(span) * 2 ** mz > max_terms:
        raise ValueError(f"{code.name}: {len(span)} codewords times 2^{mz} syndromes is too many")
    X = span.astype(np.uint8) ^ (bit * xl)
    T = _mod2_matmul((np.arange(2 ** mz)[:, None] >> np.arange(mz)) & 1, Hz).astype(np.float64)
    f = 1.0 - 2.0 * np.broadcast_to(np.asarray(p_reset, float), (n,))
    d = np.clip(np.real(walsh_hadamard(_mean_character(T, X, f))) / 2 ** mz, 0.0, None)
    return d / d.sum()


def least_weight_x(code):
    """For every Z syndrome, the least weight of an X error that has it and one such error as a qubit mask, by a
    breadth-first search over the syndromes."""
    Hz = css_parts(code)["Hz"]
    col = [int(sum(int(Hz[j, q]) << j for j in range(len(Hz)))) for q in range(code.n)]
    weight, error = np.full(2 ** len(Hz), -1, np.int64), np.zeros(2 ** len(Hz), np.int64)
    weight[0], frontier = 0, [0]
    while frontier:
        nxt = []
        for s in frontier:
            for q in range(code.n):
                t = s ^ col[q]
                if weight[t] < 0:
                    weight[t], error[t] = weight[s] + 1, error[s] | (1 << q)
                    nxt.append(t)
        frontier = nxt
    return weight, error


def bposd(H, prior, max_iter=48, osd_order=7):
    """A syndrome -> error decoder from ldpc's BP+OSD (Roffe et al.), under ldpc 2 or ldpc 0.1."""
    import ldpc
    from .css import _rank
    order = max(0, min(osd_order, H.shape[1] - _rank(H)))
    prior = [float(v) for v in np.broadcast_to(np.asarray(prior, float), (H.shape[1],))]
    if hasattr(ldpc, "BpOsdDecoder"):
        dec = ldpc.BpOsdDecoder(H.astype(np.uint8), error_channel=prior, max_iter=max_iter,
                                bp_method="minimum_sum", osd_method="osd_cs", osd_order=order)
    else:
        dec = ldpc.bposd_decoder(H.astype(int), channel_probs=prior, max_iter=max_iter, bp_method="ms",
                                 osd_method="osd_cs", osd_order=order)
    return lambda s: np.asarray(dec.decode(np.asarray(s, np.uint8)), np.uint8)


def matching(H, prior):
    """A batch syndrome -> error decoder from PyMatching (Higgott and Gidney), for checks where every qubit touches at
    most two of them."""
    import pymatching
    p = np.broadcast_to(np.asarray(prior, float), (H.shape[1],))
    m = pymatching.Matching.from_check_matrix(H, weights=np.log((1 - p) / p))
    return m


def decoded_parity(code, syndromes, p_reset, decoder="bposd"):
    """Per shot, the parity of the decoded X error against the active Z logical."""
    P = css_parts(code)
    Hz, zl = P["Hz"], P["zl"].astype(np.int64)
    s = np.asarray(syndromes, np.uint8)
    uniq, inverse = np.unique(s, axis=0, return_inverse=True)
    prior = np.asarray(p_reset, float) / 2
    if decoder == "matching":
        E = matching(Hz, prior).decode_batch(uniq)
    else:
        dec = bposd(Hz, prior) if decoder == "bposd" else decoder
        E = np.array([dec(row) for row in uniq], np.uint8).reshape(len(uniq), code.n)
    return (np.asarray(E, np.int64) @ zl % 2)[np.ravel(inverse)]


def statistic_values(code, s, kind="parity", p_reset=None, decoder="bposd"):
    """One integer per shot: "parity" the decoded parity, "z_weight" the number of Z checks that fire, "parity_weight"
    the two together."""
    if kind == "z_weight":
        return s.sum(axis=1).astype(np.int64)
    par = decoded_parity(code, s, p_reset, decoder)
    if kind == "parity":
        return par
    if kind == "parity_weight":
        return 2 * s.sum(axis=1).astype(np.int64) + par
    raise ValueError(f"unknown statistic {kind!r}")


def sample_counts(code, p_reset, shots, kind="z_weight", decoder="bposd", seed=0, chunk=250000):
    """(c0, c1): counts of one statistic over the shots of |0_L> and |1_L>, drawn a chunk at a time so that no more
    than one chunk of syndromes is held."""
    size = 2 * len(css_parts(code)["Hz"]) + 2
    counts = [np.zeros(size, np.int64), np.zeros(size, np.int64)]
    for j, start in enumerate(range(0, shots, chunk)):
        s0, s1 = sample_syndromes(code, p_reset, min(chunk, shots - start), seed=(seed, j) if j else seed)
        for c, s in zip(counts, (s0, s1)):
            c += np.bincount(statistic_values(code, s, kind, p_reset, decoder), minlength=size)
    last = max(int(np.flatnonzero(counts[0] + counts[1]).max()) + 1, 2)
    return counts[0][:last], counts[1][:last]


def leak_from_z(code, p_reset, shots, kind="z_weight", decoder="bposd", seed=0, boots=2000, splits=20,
                null_reps=40):
    """The leak through one statistic of sampled Z syndromes: plug-in distance with its permutation floor and p-value,
    and the distance an observer achieves on held-out shots."""
    c0, c1 = sample_counts(code, p_reset, shots, kind, decoder, seed)
    out = {"shots": shots, "statistic": kind}
    out.update(estimate.permutation_test(c0, c1, boots, seed))
    out.update(estimate.achieved_distance(c0, c1, splits, seed + 1, null_reps=null_reps))
    return out


def _gf2_solvable(A, b):
    """True when A y = b has a solution over GF(2)."""
    M = np.concatenate([A % 2, (b % 2)[:, None]], axis=1).astype(np.uint8)
    rows, cols = M.shape
    r = 0
    for c in range(cols - 1):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
    return not any(M[i, -1] and not M[i, :-1].any() for i in range(rows))


def _milp_min_weight(n, parity_rows, rhs, extra=None, time_limit=None, count="weight"):
    """min |a| over binary a with parity_rows a = rhs mod 2, by HiGHS through scipy.optimize.milp. With extra, a lies
    inside the support of a product of the rows of extra, and count="checks" minimises the rows in that product."""
    try:
        from scipy.optimize import Bounds, LinearConstraint, milp
    except ImportError as e:
        raise ImportError("leak_order and z_distance need scipy 1.9 or later: pip install scipy") from e
    A = np.asarray(parity_rows, float)
    m = len(A)
    # variables: a (n), k (m) with A a - 2 k = rhs
    nv = n + m
    rows, lo, hi = [], [], []
    for i in range(m):
        row = np.zeros(nv)
        row[:n] = A[i]
        row[n + i] = -2.0
        rows.append(row)
        lo.append(rhs[i])
        hi.append(rhs[i])
    ub = [1.0] * n + [np.floor(A[i].sum() / 2) + 1 for i in range(m)]
    if extra is not None:
        Hz = np.asarray(extra, float)
        mz = len(Hz)
        base = nv
        nv = base + mz + 2 * n          # y (mz), h (n), t (n): Hz^T y - 2 t - h = 0, a - h <= 0
        rows = [np.concatenate([r, np.zeros(mz + 2 * n)]) for r in rows]
        for q in range(n):
            row = np.zeros(nv)
            row[base:base + mz] = Hz[:, q]
            row[base + mz + q] = -1.0
            row[base + mz + n + q] = -2.0
            rows.append(row)
            lo.append(0.0)
            hi.append(0.0)
            row = np.zeros(nv)
            row[q] = 1.0
            row[base + mz + q] = -1.0
            rows.append(row)
            lo.append(-np.inf)
            hi.append(0.0)
        ub += [1.0] * mz + [1.0] * n + [np.floor(Hz[:, q].sum() / 2) + 1 for q in range(n)]
    c = np.zeros(nv)
    if count == "checks":
        c[base:base + mz] = 1.0
    else:
        c[:n] = 1.0
    options = {"disp": False}
    if time_limit is not None:
        options["time_limit"] = float(time_limit)
    res = milp(c, constraints=LinearConstraint(np.array(rows), lo, hi), integrality=np.ones(nv),
               bounds=Bounds(np.zeros(nv), np.array(ub)), options=options)
    a = None if res.x is None else np.rint(res.x[:n]).astype(np.uint8)
    bound = getattr(res, "mip_dual_bound", None)
    y = None if res.x is None or extra is None else np.rint(res.x[base:base + mz]).astype(np.uint8)
    return {"status": int(res.status), "weight": None if a is None else int(a.sum()), "a": a, "y": y,
            "dual_bound": None if bound is None or not np.isfinite(bound) else float(bound),
            "optimal": res.status == 0, "infeasible": res.status == 2}


def leak_order(code, time_limit=None):
    """The leak order under relaxation: the least weight of a Z logical of the active qubit lying inside the support of
    some Z stabilizer, or None where no such logical exists and the code is protected."""
    P = css_parts(code)
    n = code.n
    rows = list(P["Hx"]) + [P["xl"]] + list(P["others"])
    rhs = [0.0] * len(P["Hx"]) + [1.0] + [0.0] * len(P["others"])
    first = _milp_min_weight(n, rows, rhs, time_limit=time_limit)
    if first["a"] is not None and first["optimal"]:
        support = np.flatnonzero(first["a"])
        if _gf2_solvable(P["Hz"][:, support].T.astype(np.uint8), np.ones(len(support), np.uint8)):
            return {"order": first["weight"], "a": first["a"], "optimal": True, "dual_bound": first["weight"]}
    full = _milp_min_weight(n, rows, rhs, extra=P["Hz"], time_limit=time_limit)
    if full["infeasible"]:
        return {"order": None, "a": None, "optimal": True, "dual_bound": None}
    return {"order": full["weight"], "a": full["a"], "optimal": full["optimal"], "dual_bound": full["dual_bound"]}


def z_distance(code, time_limit=None):
    """The least weight of a nontrivial Z logical, over all logical qubits, by integer programming: {"distance",
    "optimal", "dual_bound"}."""
    P = css_parts(code)
    pairs = list(getattr(code, "logical_pairs", [(code.zl_str, code.xl_str)]))
    xs = [np.array([1 if ch == "X" else 0 for ch in x], np.uint8) for _z, x in pairs]
    best, optimal, bound = None, True, None
    for j in range(len(xs)):
        rows = list(P["Hx"]) + [xs[j]]
        r = _milp_min_weight(code.n, rows, [0.0] * len(P["Hx"]) + [1.0], time_limit=time_limit)
        if r["weight"] is not None and (best is None or r["weight"] < best):
            best = r["weight"]
        optimal = optimal and r["optimal"]
        if r["dual_bound"] is not None:
            bound = r["dual_bound"] if bound is None else min(bound, r["dual_bound"])
    return {"distance": best, "optimal": optimal, "dual_bound": bound}
