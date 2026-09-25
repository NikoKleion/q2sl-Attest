# syndrome_leakage.eavesdrop: error of an observer reading only the syndrome record
import math

import numpy as np

from .core import tvd
from . import channels as _ch
from .analyze import population_leak, eavesdropper_error


def chernoff_exponent(d0, d1, grid=201, tol=1e-12):
    # probabilities at or below tol are zero; the exact engine leaves residue near 1e-16
    d0 = np.asarray(d0, float); d1 = np.asarray(d1, float)
    m = (d0 > tol) & (d1 > tol)
    if not m.any():
        return math.inf, 0.5
    p = d0[m]; q = d1[m]
    ss = np.linspace(0.0, 1.0, grid)
    vals = np.array([np.sum(p ** s * q ** (1.0 - s)) for s in ss])
    i = int(np.argmin(vals))
    return float(-math.log(max(vals[i], 1e-300))), float(ss[i])


def bhattacharyya_coefficient(d0, d1):
    return float(np.sum(np.sqrt(np.asarray(d0, float) * np.asarray(d1, float))))


def ml_attack(d0, d1, rounds, trials=4000, seed=0):
    d0 = np.asarray(d0, float); d1 = np.asarray(d1, float)
    rng = np.random.default_rng(seed)
    k = len(d0)
    with np.errstate(divide="ignore"):
        llr = np.log(np.where(d1 > 0, d1, 1e-300)) - np.log(np.where(d0 > 0, d0, 1e-300))
    half = trials // 2
    s0 = rng.choice(k, size=(half, rounds), p=d0 / d0.sum())
    s1 = rng.choice(k, size=(trials - half, rounds), p=d1 / d1.sum())
    t0 = llr[s0].sum(axis=1)
    t1 = llr[s1].sum(axis=1)
    err0 = float(np.mean(t0 > 0) + 0.5 * np.mean(t0 == 0))
    err1 = float(np.mean(t1 < 0) + 0.5 * np.mean(t1 == 0))
    return 0.5 * (err0 + err1)


LEAK_FLOOR = 1e-12


def rounds_for_error(d0, d1, target=0.01):
    if tvd(d0, d1) < LEAK_FLOOR:
        return None
    C, _s = chernoff_exponent(d0, d1)
    if not math.isfinite(C) or C <= 0:
        return None
    return int(math.ceil(math.log(1.0 / target) / C))


def attack_point(code, gamma=0.2, channel="amplitude_damping", rounds=(1, 10, 50, 200),
                 trials=4000, seed=0):
    kraus = _ch.LIBRARY[channel](gamma)
    leak, d0, d1 = population_leak(code, kraus)
    C, s_star = chernoff_exponent(d0, d1)
    rho = bhattacharyya_coefficient(d0, d1)
    rows = []
    for r in rounds:
        rows.append({"rounds": int(r),
                     "attack_error": ml_attack(d0, d1, r, trials=trials, seed=seed + r),
                     "bhattacharyya_bound": eavesdropper_error(d0, d1, r),
                     "chernoff_error": 0.5 * math.exp(-C * r)})
    return {"code": code.name, "gamma": float(gamma), "channel": channel, "tvd": float(leak),
            "chernoff": C, "s_star": s_star, "bhattacharyya_coefficient": rho,
            "rounds_for_1pct": rounds_for_error(d0, d1, 0.01), "rows": rows}


def report(code=None, gamma=0.2, trials=4000, seed=0):
    from .codes import STANDARD
    code = code if code is not None else STANDARD["hamming_7"]()
    res = attack_point(code, gamma=gamma, trials=trials, seed=seed)
    L = ["Likelihood-ratio test on the syndrome record",
         f"  code {res['code']}, {res['channel']}(gamma={res['gamma']}), syndrome TVD per shot {res['tvd']:.3e}",
         f"  Chernoff exponent {res['chernoff']:.4f} nats/shot at s* = {res['s_star']:.2f}; "
         f"Bhattacharyya coefficient {res['bhattacharyya_coefficient']:.6f}",
         "",
         "  rounds   ML attack error   Chernoff rate   Bhattacharyya bound",
         "  " + "-" * 62]
    for r in res["rows"]:
        L.append(f"  {r['rounds']:6d}   {r['attack_error']:15.4f}   {r['chernoff_error']:13.4f}   "
                 f"{r['bhattacharyya_bound']:19.4f}")
    need = res["rounds_for_1pct"]
    reach = (f"  reaching 1 percent error needs about {need} syndrome rounds at this gamma."
             if need is not None else
             " no round count reaches 1 percent error: the two syndrome distributions are equal.")
    L += ["", reach]
    return "\n".join(L)


def state_grid(n_theta=7, n_phi=6):
    import itertools
    thetas = [math.pi * i / (n_theta - 1) for i in range(n_theta)]
    phis = [2.0 * math.pi * j / n_phi for j in range(n_phi)]
    out = []
    for t in thetas:
        if t == 0.0 or abs(t - math.pi) < 1e-12:
            out.append((t, 0.0))
        else:
            out.extend((t, ph) for ph in phis)
    return out


def worst_case_pair(code, gamma=0.2, channel="amplitude_damping", grid=None):
    kraus = _ch.LIBRARY[channel](gamma)
    states = grid if grid is not None else state_grid()
    dists = [(t, ph, code.syndrome_dist(code.apply(code.logical_state(t, ph), kraus))) for t, ph in states]
    best = None
    for i in range(len(dists)):
        for j in range(i + 1, len(dists)):
            t0, p0, d0 = dists[i]
            t1, p1, d1 = dists[j]
            C, _s = chernoff_exponent(d0, d1)
            if best is None or C > best["chernoff"]:
                best = {"chernoff": C, "tvd": float(tvd(d0, d1)),
                        "state_a": (t0, p0), "state_b": (t1, p1),
                        "rounds_for_1pct": rounds_for_error(d0, d1, 0.01)}
    axis = population_leak(code, kraus)
    C_axis, _s = chernoff_exponent(axis[1], axis[2])
    best["axis_chernoff"] = C_axis
    best["axis_tvd"] = float(axis[0])
    best["axis_is_worst"] = bool(C_axis >= best["chernoff"] - 1e-9)
    return best


def worst_case_report(gamma=0.2):
    from .codes import STANDARD
    L = ["Largest syndrome distance over pairs of logical states",
         f"  amplitude_damping(gamma={gamma}), searched over a Bloch-sphere grid of logical states",
         "",
         "  code                worst TVD   |0>/|1> TVD   worst pair (theta/pi, phi/pi)   axis optimal",
         "  " + "-" * 88]
    for name in ("repetition", "code_4_1_2", "hamming_7", "steane", "five_qubit"):
        w = worst_case_pair(STANDARD[name](), gamma=gamma)
        a = w["state_a"]; b = w["state_b"]
        pair = f"({a[0]/math.pi:.2f},{a[1]/math.pi:.2f}) vs ({b[0]/math.pi:.2f},{b[1]/math.pi:.2f})"
        L.append(f"  {name:<18} {w['tvd']:9.3e}   {w['axis_tvd']:11.3e}   {pair:<30}   "
                 f"{'yes' if w['axis_is_worst'] else 'no'}")
    return chr(10).join(L)


def recovery_table(code):
    from .core import op
    keys = list(code.PROJ.keys())
    table = {k: None for k in keys}
    ident = "I" * code.n
    table[tuple(0 for _ in code.GS)] = np.eye(code.dim, dtype=complex)
    for q in range(code.n):
        for pauli in ("X", "Y", "Z"):
            estr = ident[:q] + pauli + ident[q + 1:]
            E = op(estr)
            bits = []
            for g in code.GS:
                comm = np.allclose(E @ g, g @ E)
                bits.append(0 if comm else 1)
            k = tuple(bits)
            if k in table and table[k] is None:
                table[k] = E
    for k in keys:
        if table[k] is None:
            table[k] = np.eye(code.dim, dtype=complex)
    return [table[k] for k in keys]


def _measurement_channel(code, rho, recov=None):
    projs = list(code.PROJ.values())
    out = np.zeros_like(rho)
    for i, P in enumerate(projs):
        branch = P @ rho @ P.conj().T
        if recov is not None:
            R = recov[i]
            branch = R @ branch @ R.conj().T
        out += branch
    return out


def repeated_extraction(code, gamma=0.2, channel="amplitude_damping", rounds=8, seed=0, correct=False,
                        trajectories=None, kraus_override=None):
    # kraus_override takes a channel outside the name library, e.g. channel_from_t1t2
    kraus = kraus_override if kraus_override is not None else _ch.LIBRARY[channel](gamma)
    recov = recovery_table(code) if correct else None
    states = {}
    for label, theta in ((0, 0.0), (1, math.pi)):
        states[label] = code.logical_state(theta, 0.0)
    per_round = []
    for r in range(rounds):
        ds = {}
        for label in (0, 1):
            rho = code.apply(states[label], kraus)
            ds[label] = code.syndrome_dist(rho)
            states[label] = _measurement_channel(code, rho, recov=recov)
            tr = float(np.real(np.trace(states[label])))
            if tr > 1e-15:
                states[label] = states[label] / tr
        C, _s = chernoff_exponent(ds[0], ds[1])
        per_round.append({"round": r + 1, "tvd": float(tvd(ds[0], ds[1])), "chernoff": C})
    return {"code": code.name, "gamma": float(gamma), "corrected": bool(correct), "rounds": per_round}


def repeated_extraction_report(code=None, gamma=0.2, rounds=8, correct=True, seed=0, trajectories=None):
    from .codes import STANDARD
    code = code if code is not None else STANDARD["hamming_7"]()
    res = repeated_extraction(code, gamma=gamma, rounds=rounds, seed=seed, correct=correct)
    first = res["rounds"][0]["tvd"] if res["rounds"] else 0.0
    L = ["Repeated syndrome extraction on one held logical state",
         f"  code {res['code']}, amplitude_damping(gamma={res['gamma']}), one logical state held across rounds,",
         f"  each round: noise, syndrome measurement, {'recovery' if correct else 'no recovery'}. "
         f"Ensemble evolved exactly, no sampling.",
         "",
         "  round   syndrome TVD   ratio to round 1",
         "  " + "-" * 44]
    for r in res["rounds"]:
        ratio = r["tvd"] / first if first > 0 else 0.0
        L.append(f"  {r['round']:5d}   {r['tvd']:12.3e}   {ratio:16.3f}")
    return chr(10).join(L)


def pauli_control(code=None, p=0.1, rounds=200, trials=4000, seed=0):
    from .codes import STANDARD
    code = code if code is not None else STANDARD["hamming_7"]()
    kraus = _ch.depolarizing(p)
    _leak, d0, d1 = population_leak(code, kraus)
    return ml_attack(d0, d1, rounds, trials=trials, seed=seed), float(tvd(d0, d1))


def main():
    print(report())


if __name__ == "__main__":
    main()
