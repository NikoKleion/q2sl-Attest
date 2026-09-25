# syndrome records by stabilizer simulation under T1 and T2 relaxation, for codes of any size
import math

import numpy as np

from . import estimate


def relaxation_mixture(t, T1, T2):
    """(p_reset, p_z) for relaxation over time t."""
    t, T1, T2 = (np.asarray(v, float) for v in (t, T1, T2))
    if np.any(T2 > T1 * (1 + 1e-12)):
        raise ValueError("T2 > T1: relaxation is not a mixture of reset and Z on that qubit")
    p_reset = 1.0 - np.exp(-t / T1)
    p_z = (1.0 - np.exp(-t * np.maximum(1.0 / T2 - 1.0 / T1, 0.0))) / 2.0
    if p_reset.ndim == 0:
        return float(p_reset), float(p_z)
    return p_reset, p_z


def mixture_channel(p_reset, p_z=0.0):
    """Kraus operators of reset with probability p_reset and Z with probability p_z, for the exact engines."""
    I2 = np.eye(2, dtype=complex)
    Z = np.diag([1.0, -1.0]).astype(complex)
    to0 = [np.array([[1, 0], [0, 0]], complex), np.array([[0, 1], [0, 0]], complex)]
    out = []
    for w, U in ((1.0 - p_z, I2), (p_z, Z)):
        out.append(math.sqrt(w * (1.0 - p_reset)) * U)
        out += [math.sqrt(w * p_reset) * U @ k for k in to0]
    return out


def _per_qubit(value, n, name):
    v = np.asarray(value, float)
    if v.ndim == 0:
        return np.full(n, float(v))
    if v.shape != (n,):
        raise ValueError(f"{name}: expected one value or {n}, got shape {v.shape}")
    return v


def _logical_state(code, bit):
    """A Stim simulator holding |0_L> or |1_L>: +1 on every generator and every logical Z, except the active logical Z,
    whose sign is (-1)^bit."""
    import stim
    zl = list(getattr(code, "zl_strs", [code.zl_str]))
    active = zl.index(code.zl_str)
    stabs = [stim.PauliString(s) for s in code.stab_strings]
    stabs += [stim.PauliString(("-" if (l == active and bit) else "+") + z) for l, z in enumerate(zl)]
    sim = stim.TableauSimulator()
    sim.set_state_from_stabilizers(stabs, allow_redundant=True)
    return sim


def _run(code, bit, reset, z, coins):
    """Records of one logical state: reset[i, q] and z[i, q] mark the operations on shot i, and coins[i] supplies every
    outcome the simulator leaves random, first for the resets and then the generators."""
    import stim
    base = _logical_state(code, bit)
    gens = [stim.PauliString(s) for s in code.stab_strings]
    n, m = code.n, len(gens)
    rec = np.zeros((reset.shape[0], m), np.uint8)
    for i in range(reset.shape[0]):
        sim = base.copy()
        for q in np.flatnonzero(reset[i]):
            q = int(q)
            v = sim.peek_z(q)
            one = bool(coins[i, q]) if v == 0 else v < 0
            if v == 0:
                sim.postselect_z(q, desired_value=one)
            if one:
                sim.x(q)
        zq = [int(q) for q in np.flatnonzero(z[i])]
        if zq:
            sim.z(*zq)
        for j, g in enumerate(gens):
            v = sim.peek_observable_expectation(g)
            if v == 0:
                b = bool(coins[i, n + j])
                sim.postselect_observable(g, desired_value=b)
            else:
                b = v < 0
            rec[i, j] = b
    return rec


def sample_records(code, p_reset, p_z=0.0, shots=10000, seed=0, readout=None):
    """(r0, r1): syndrome records of |0_L> and |1_L>, one row per shot and one column per generator, a 1 where the
    generator reads -1."""
    n, m = code.n, len(code.stab_strings)
    pr, pz = _per_qubit(p_reset, n, "p_reset"), _per_qubit(p_z, n, "p_z")
    out = []
    for bit in (0, 1):
        rng = np.random.default_rng([seed, bit])
        reset = rng.random((shots, n)) < pr
        z = rng.random((shots, n)) < pz
        coins = rng.random((shots, n + m)) < 0.5
        rec = _run(code, bit, reset, z, coins)
        if readout is not None:
            rec ^= (rng.random((shots, m)) < _per_qubit(readout, m, "readout")).astype(np.uint8)
        out.append(rec)
    return out[0], out[1]


def statistic(records, code, kind="z_weight"):
    """One integer per shot."""
    r = np.asarray(records)
    if callable(kind):
        return np.asarray(kind(r), np.int64)
    if kind == "syndrome":
        if r.shape[1] > 24:
            raise ValueError(f"{r.shape[1]} generators is too many outcomes for the full syndrome")
        return (r.astype(np.int64) << np.arange(r.shape[1])).sum(axis=1)
    select = {"weight": lambda s: True,
              "z_weight": lambda s: set(s) <= {"I", "Z"},
              "x_weight": lambda s: set(s) <= {"I", "X"}}[kind]
    mask = np.array([select(s) for s in code.stab_strings])
    return r[:, mask].sum(axis=1).astype(np.int64)


def leak_from_records(r0, r1, code, kind="z_weight", boots=2000, splits=20, seed=0, null_reps=40):
    """The leak read through one statistic of the records: plug-in distance, permutation floor and p-value, and the
    distance an observer achieves on held-out shots."""
    s0, s1 = statistic(r0, code, kind), statistic(r1, code, kind)
    size = int(max(s0.max(), s1.max())) + 1
    c0, c1 = np.bincount(s0, minlength=size), np.bincount(s1, minlength=size)
    out = {"shots": (len(s0), len(s1)), "statistic": kind if isinstance(kind, str) else "custom"}
    out.update(estimate.permutation_test(c0, c1, boots, seed))
    out.update(estimate.achieved_distance(c0, c1, splits, seed + 1, null_reps=null_reps))
    return out


def sample_floor(p, shots, trials=200, seed=0):
    """Mean distance between a sample of `shots` draws and the distribution it came from."""
    rng = np.random.default_rng(seed)
    p = np.asarray(p, float)
    return float(np.mean([0.5 * np.abs(rng.multinomial(shots, p) / shots - p).sum() for _ in range(trials)]))


def fit_p_value(counts, p, trials=2000, seed=0):
    """Probability that a sample of the same size drawn from p lies at least as far from p as counts do."""
    rng = np.random.default_rng(seed)
    p = np.asarray(p, float)
    c = np.asarray(counts, float)
    shots = int(c.sum())
    seen = 0.5 * np.abs(c / shots - p).sum()
    draws = 0.5 * np.abs(rng.multinomial(shots, p, size=trials) / shots - p).sum(axis=1)
    return float((draws >= seen).mean())
