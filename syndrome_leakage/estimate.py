# the leak from measured records: plug-in distance, floor, permutation test, achieved distance
import numpy as np

from .core import tvd


def outcome_index(key):
    """The integer outcome behind a result key: an int, a hex string, or bits with register spaces."""
    if isinstance(key, (int, np.integer)):
        return int(key)
    text = str(key).strip().replace(" ", "").replace("_", "")
    if text.lower().startswith("0x"):
        return int(text, 16)
    if text.lower().startswith("0b"):
        return int(text, 2)
    return int(text, 2)


def dist_from_counts(counts, n_outcomes):
    """Counts -> (distribution, shots); keys may be bit strings, hex strings or ints."""
    if not hasattr(counts, "items"):
        return dist_from_memory(counts, n_outcomes)
    d = np.zeros(n_outcomes, float)
    for key, n in counts.items():
        d[outcome_index(key)] += n
    total = d.sum()
    return d / max(total, 1), int(total)


def dist_from_memory(memory, n_outcomes):
    """Per-shot results, one key per shot -> (distribution, shots)."""
    d = np.zeros(n_outcomes, float)
    for key in memory:
        d[outcome_index(key)] += 1
    total = d.sum()
    return d / max(total, 1), int(total)


def null_floor(p, shots, trials=400, seed=0, transform=None):
    """Mean distance between two samples of one distribution, the reading at zero leak."""
    rng = np.random.default_rng(seed)
    f = transform or (lambda x: x)
    return float(np.mean([tvd(f(rng.multinomial(shots, p) / shots), f(rng.multinomial(shots, p) / shots))
                          for _ in range(trials)]))


def shots_for_leak(p0, p1, factor=2.0, grid=(1000, 2000, 4000, 8000, 16000, 32000, 64000),
                   trials=200, seed=0, transform=None):
    """Smallest shot count in the grid whose null floor is under the leak by the given factor."""
    f = transform or (lambda x: x)
    leak = tvd(f(p0), f(p1))
    for n in grid:
        if leak >= factor * null_floor(p0, n, trials, seed, transform):
            return n
    return None


def permutation_test(c0, c1, boots=2000, seed=0, transform=None):
    """p-value for one distance under random reassignment of the pooled shots."""
    rng = np.random.default_rng(seed)
    f = transform or (lambda x: x)
    c0, c1 = np.asarray(c0, np.int64), np.asarray(c1, np.int64)
    n0, pooled = int(c0.sum()), c0 + c1
    raw = tvd(f(c0 / max(n0, 1)), f(c1 / max(int(c1.sum()), 1)))
    vals = []
    for _ in range(boots):
        a = rng.multivariate_hypergeometric(pooled, n0)
        b = pooled - a
        vals.append(tvd(f(a / max(a.sum(), 1)), f(b / max(b.sum(), 1))))
    vals = np.array(vals)
    return {"tvd": raw, "null": float(vals.mean()), "null_95": float(np.percentile(vals, 95)),
            "p_value": float((vals >= raw).mean())}


def achieved_distance(c0, c1, splits=20, seed=0, transform=None, prior=0.5, null_reps=40):
    """Distance an eavesdropper achieves: fit the likelihood ratio on half of each record, score the rest."""
    rng = np.random.default_rng(seed)
    f = transform or (lambda x: x)
    c0, c1 = np.asarray(c0, np.int64), np.asarray(c1, np.int64)
    errs = []
    for _ in range(splits):
        for fit, test in ((0, 1), (1, 0)):
            a0 = rng.multivariate_hypergeometric(c0, int(c0.sum()) // 2)
            a1 = rng.multivariate_hypergeometric(c1, int(c1.sum()) // 2)
            halves = ((a0, c0 - a0), (a1, c1 - a1))
            q0 = f(halves[0][fit] + prior) / f(halves[0][fit] + prior).sum()
            q1 = f(halves[1][fit] + prior) / f(halves[1][fit] + prior).sum()
            t0, t1 = f(halves[0][test]), f(halves[1][test])
            call = np.where(q1 > q0, 1.0, np.where(q1 < q0, 0.0, 0.5))
            e0 = float((t0 * call).sum() / max(t0.sum(), 1))
            e1 = float((t1 * (1 - call)).sum() / max(t1.sum(), 1))
            errs.append(0.5 * (e0 + e1))
    errs = np.array(errs)
    n_test = (int(c0.sum()) + int(c1.sum())) // 2
    err = float(errs.mean())
    se = float(np.sqrt(max(err * (1 - err), 1e-12) / n_test))
    out = {"attack_error": err, "attack_se": se, "distance": max(1.0 - 2.0 * err, 0.0),
           "distance_ci": (max(1.0 - 2.0 * (err + 1.96 * se), 0.0), max(1.0 - 2.0 * (err - 1.96 * se), 0.0)),
           "split_sd": float(errs.std())}
    if null_reps:
        pooled, n0 = c0 + c1, int(c0.sum())
        nulls = []
        for _ in range(null_reps):
            a = rng.multivariate_hypergeometric(pooled, n0)
            nulls.append(achieved_distance(a, pooled - a, max(splits // 4, 1), int(rng.integers(1 << 30)),
                                           transform, prior, null_reps=0)["distance"])
        out["distance_null"] = float(np.mean(nulls))
        out["distance_corrected"] = max(out["distance"] - out["distance_null"], 0.0)
    return out


def leak_from_counts(counts0, counts1, n_outcomes, transform=None, boots=2000, splits=20, seed=0,
                     null_reps=40):
    """Leak between two measured records: plug-in distance, floor, permutation p, achieved distance."""
    p0, n0 = dist_from_counts(counts0, n_outcomes)
    p1, n1 = dist_from_counts(counts1, n_outcomes)
    return leak_from_dists(p0, n0, p1, n1, transform, boots, splits, seed, null_reps)


def leak_from_dists(p0, n0, p1, n1, transform=None, boots=2000, splits=20, seed=0, null_reps=40):
    """leak_from_counts on distributions and shot counts already in hand."""
    c0 = np.rint(np.asarray(p0) * n0).astype(np.int64)
    c1 = np.rint(np.asarray(p1) * n1).astype(np.int64)
    out = {"shots": (int(c0.sum()), int(c1.sum()))}
    out.update(permutation_test(c0, c1, boots, seed, transform))
    out.update(achieved_distance(c0, c1, splits, seed + 1, transform, null_reps=null_reps))
    return out


def report(rows, title="leak from measured records"):
    """One line per row of (label, result of leak_from_counts)."""
    L = [title, ""]
    head = (f"  {'record':>12} {'shots':>14} {'distance':>9} {'floor':>7} {'p':>7} "
            f"{'achieved':>9} {'corrected':>10} {'95% interval':>18}")
    L.append(head)
    L.append("  " + "-" * (len(head) - 2))
    for label, r in rows:
        L.append(f"  {label:>12} {str(r['shots']):>14} {r['tvd']:9.4f} {r['null']:7.4f} "
                 f"{r['p_value']:7.4f} {r['distance']:9.4f} {r['distance_corrected']:10.4f} "
                 f"{'(%.4f, %.4f)' % r['distance_ci']:>18}")
    L.append("")
    L.append("  distance is the plug-in reading and sits above the floor of its shot count even at zero leak")
    L.append("  achieved: the likelihood ratio fitted on half of each record and scored on the other half")
    return "\n".join(L)
