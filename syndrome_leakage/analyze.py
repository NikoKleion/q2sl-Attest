# leak measures in exact and analytic modes
import math
import itertools
import numpy as np

from .core import tvd
from . import channels as _ch


def _group(strings, n):
    from .core import pmul
    group = {"I" * n}
    frontier = list(strings)
    changed = True
    while changed:
        changed = False
        for a in list(group):
            for g in frontier:
                p = pmul(a, g)
                if p not in group:
                    group.add(p); changed = True
    return group


_TOO_LARGE = object()


def _span(rows):
    # every element of the GF(2) span of these bit masks
    out = [0]
    for r in rows:
        out += [x ^ r for x in out]
    return out


def _css_leak(code, V, max_weight, max_group):
    # the same search on a CSS code, over the Z-side span
    import numpy as np
    from .css import css_matrices
    try:
        _Hx, Hz = css_matrices(code.stab_strings)
    except ValueError:
        return None
    rows = [int("".join("1" if b else "0" for b in r[::-1]), 2) for r in Hz if np.any(r)]
    zl = int("".join("1" if c == "Z" else "0" for c in code.zl_str[::-1]), 2)
    if 2 ** (len(rows) + 1) > max_group:
        return _TOO_LARGE
    stab = _span(rows)
    supports, drop = set(stab), set()
    for s in supports:
        for t in supports:
            if s != t and s & ~t == 0:
                drop.add(s)
                break
    maximal = sorted(supports - drop, key=lambda x: -bin(x).count("1"))
    vmask = 0
    for i in V:
        vmask |= 1 << i
    cap = max_weight if max_weight is not None else code.n
    best = None
    for a in _span(rows + [zl]):
        w = bin(a).count("1")
        if a == 0 or w > cap or (best is not None and w >= best):
            continue
        if bin(a & vmask).count("1") % 2 == 1 and any(a & ~s == 0 for s in maximal):
            best = w
    return (best is not None), best


def analytic_leak(code, max_weight=None, max_group=2 ** 20):
    # amplitude-damping leak order: least |A| with Z_A in <S, Z_L> and |A cap supp(X_L)| odd
    n = code.n
    V = set(i for i, c in enumerate(code.xl_str) if c in "XY")
    css = _css_leak(code, V, max_weight, max_group)
    assert css is not _TOO_LARGE, (
        f"the Z-side span of this code holds more than {max_group} elements; raise max_group")
    if css is not None:
        return css
    assert 2 ** len(code.stab_strings) <= max_group, (
        f"{len(code.stab_strings)} generators is {2 ** len(code.stab_strings)} group elements; "
        f"raise max_group or use a CSS code")
    G_S = _group(code.stab_strings, n)
    G_0 = _group(code.stab_strings + [code.zl_str], n)
    cap = max_weight if max_weight is not None else n
    min_order = None
    for g in G_S:
        zsupp = [i for i, c in enumerate(g) if c == "Z"]
        for r in range(1, min(len(zsupp), cap) + 1):
            if min_order is not None and r >= min_order:
                break
            for A in itertools.combinations(zsupp, r):
                if len(set(A) & V) % 2 == 1:
                    zstr = "".join("Z" if i in A else "I" for i in range(n))
                    if zstr in G_0:
                        min_order = r if min_order is None else min(min_order, r)
                        break
    return (min_order is not None), min_order


def population_leak(code, kraus, readout=None):
    d0 = code.syndrome_dist(code.apply(code.logical_state(0.0, 0.0), kraus))
    d1 = code.syndrome_dist(code.apply(code.logical_state(math.pi, 0.0), kraus))
    if readout is not None:
        from .expectations import apply_readout
        d0, d1 = apply_readout(d0, readout), apply_readout(d1, readout)
    return tvd(d0, d1), d0, d1


def l2_leak(code, kraus, readout=None):
    """The L2 distance between the two syndrome distributions, and the bounds it puts on the leak."""
    from .expectations import population_leak as expectation_leak
    _leak, d0, d1 = expectation_leak(code, kraus, readout)
    delta = np.asarray(d0, float) - np.asarray(d1, float)
    l1 = float(np.abs(delta).sum())
    l2 = float(np.sqrt((delta ** 2).sum()))
    if l2 == 0.0:
        return 0.0, 0.0, 0.0, float("nan")
    shape = l1 / l2
    return l2, l2 / 2, math.sqrt(math.ceil(shape ** 2)) * l2 / 2, shape


def phase_leak(code, kraus):
    base = code.syndrome_dist(code.apply(code.logical_state(math.pi / 2, 0.0), kraus))
    return max(float(np.abs(code.syndrome_dist(code.apply(code.logical_state(math.pi / 2, ph), kraus)) - base).max())
               for ph in (math.pi / 3, math.pi / 2, math.pi, 4 * math.pi / 3, 3 * math.pi / 2))


def pauli_null(code, p=0.1):
    d0 = code.syndrome_dist(code.apply(code.logical_state(0.0, 0.0), _ch.depolarizing(p)))
    d1 = code.syndrome_dist(code.apply(code.logical_state(math.pi, 0.0), _ch.depolarizing(p)))
    return tvd(d0, d1)


def eavesdropper_error(d0, d1, rounds):
    # Bhattacharyya bound on the error distinguishing |0_L> from |1_L> over `rounds` i.i.d. shots
    return 0.5 * float(np.sum(np.sqrt(np.asarray(d0) * np.asarray(d1)))) ** rounds


def measured_order(code, channel_fn, gammas=(0.01, 0.02, 0.04, 0.08)):
    pts = []
    for g in gammas:
        leak = population_leak(code, channel_fn(g))[0]
        pts.append((g, leak))
    slopes = []
    for i in range(len(pts) - 1):
        (g1, d1), (g2, d2) = pts[i], pts[i + 1]
        if d1 > 1e-13 and d2 > 1e-13:
            slopes.append(math.log(d2 / d1) / math.log(g2 / g1))
    return (float(np.mean(slopes)) if slopes else None), pts


class Report:
    def __init__(self, **kw):
        self.__dict__.update(kw)

    def __str__(self):
        L = []
        L.append(f"syndrome-leakage analysis: {self.code} under {self.channel}")
        L.append(f"  phase leak      : {self.phase:.2e}   ({'protected' if self.phase < 1e-9 else 'LEAKS'})")
        verdict = "LEAKS" if self.leaks else "protected"
        order = f", first at order gamma^{self.analytic_order}" if self.leaks else ""
        L.append(f"  population      : {verdict}{order}")
        L.append(f"    exact leak (gamma={self.gamma})     : {self.population:.3e}")
        if self.measured_slope is not None:
            L.append(f"    measured gamma-slope           : {self.measured_slope:.2f}   (analytic order {self.analytic_order})")
        L.append(f"    validates analytic order         : {self.validated}")
        L.append(f"  Pauli null (must be ~0)            : {self.pauli:.2e}")
        if self.leaks:
            L.append(f"  eavesdropper error, 10 / 100 rounds: {self.eve10:.2e} / {self.eve100:.2e}")
        L.append(f"  mode: exact (n={self.n} qubits) + analytic")
        return "\n".join(L)


def analyze(code, channel="amplitude_damping", gamma=0.2):
    # entry point; channel is a name in channels.LIBRARY or a Kraus list
    is_ad = isinstance(channel, str) and channel in ("amplitude_damping", "generalized_amplitude_damping")
    if isinstance(channel, str):
        fn = _ch.LIBRARY[channel]
        kraus = fn(gamma)
        cname = f"{channel}(gamma={gamma})"
    else:
        fn, kraus, cname = None, channel, "custom channel"
    code.verify_projectors()
    pop, d0, d1 = population_leak(code, kraus)
    ph = phase_leak(code, kraus)
    pnull = pauli_null(code)
    leaks, order = analytic_leak(code) if is_ad else (pop > 1e-9, None)
    slope, _ = measured_order(code, fn) if is_ad else (None, None)
    validated = (not leaks) if slope is None else (order is not None and abs(slope - order) < 0.3)
    return Report(code=code.name, channel=cname, n=code.n, gamma=gamma,
                  phase=ph, population=pop, pauli=pnull, leaks=leaks, analytic_order=order,
                  measured_slope=slope, validated=validated,
                  eve10=eavesdropper_error(d0, d1, 10), eve100=eavesdropper_error(d0, d1, 100))


def selftest():
    from .codes import STANDARD
    print("# self-test: exact simulation must validate the analytic leak order on every standard code")
    ok = True
    for name, maker in STANDARD.items():
        r = analyze(maker(), "amplitude_damping", gamma=0.2)
        good = r.validated and r.phase < 1e-9
        ok = ok and good
        ord_s = r.analytic_order if r.leaks else "-"
        slp = f"{r.measured_slope:.2f}" if r.measured_slope is not None else "-"
        print(f"#   {name:>14}: population {'LEAKS' if r.leaks else 'protected':>9}  analytic-order {str(ord_s):>2}  "
              f"measured-slope {slp:>5}  phase {r.phase:.0e}  -> {'ok' if good else 'FAIL'}")
    print(f"# all standard codes pass: {ok}")
    return ok


def _cli():
    import argparse
    from .codes import STANDARD
    ap = argparse.ArgumentParser(description="QEC syndrome information-leakage analyzer")
    ap.add_argument("code", nargs="?", choices=list(STANDARD), help="a standard code")
    ap.add_argument("channel", nargs="?", default="amplitude_damping", choices=list(_ch.LIBRARY))
    ap.add_argument("--gamma", type=float, default=0.2)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest or not a.code:
        selftest(); return
    print(analyze(STANDARD[a.code](), a.channel, a.gamma))


if __name__ == "__main__":
    _cli()
