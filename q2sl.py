# q2sl attest command line
import math
import runpy
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
def _codes():
    # read from the registries so the help can never drift from what assess accepts
    from syndrome_leakage.codes import STANDARD, STRINGS
    return tuple(STANDARD) + tuple(STRINGS)

COMMANDS = {
    "assess":      "analyses for a code name, file:<path>, qecdb:<id> or backend:<name>",
    "qecdb":       "search qecdb.org for CSS codes, e.g. q2sl qecdb n=9 k=1 d=3 family=surface",
    "leak":        "syndrome leak self-test (syndrome_leakage)",
    "fusion":      "residual min-entropy on synthetic payloads (entropy_fusion)",
    "eavesdrop":   "likelihood-ratio test on syndrome records (syndrome_leakage.eavesdrop)",
    "reconstruct": "the full reconstruction (needs torch, ldpc)",
    "example":     "the end to end run in scenario/one_posterior.py",
    "suite":       "run every suite module and print a SARIF summary",
    "all":         "run the fast numpy demos and the suite",
}


def _assess(args):
    pos = [a for a in (args or []) if not a.startswith("--")]
    target = pos[0] if pos else None
    if target is None or target in ("-h", "--help"):
        print("q2sl assess <target> [code]\n")
        print("  a code name:      " + ", ".join(_codes()))
        print("  a code file:      file:<path> or file:<Hx path>,<Hz path>  (0/1 text, .npz with hx and")
        print("                    hz, .alist, .mtx, or Pauli strings one per line)")
        print("  a qecdb.org code: qecdb:<id>  (find ids with q2sl qecdb n=9 k=1 d=3)")
        print("  a qiskit backend: backend:<name> [code]  (e.g. backend:manila steane; the leak at the")
        print("                    backend's T1 over one syndrome-extraction cycle, default code hamming_7)")
        return
    from syndrome_leakage.codes import STANDARD, STRINGS
    if target in STANDARD:
        return _assess_code(target)
    if target in STRINGS:
        return _report_code(STRINGS[target]())
    if target.startswith("file:"):
        return _assess_file(target[len("file:"):])
    if target.startswith("qecdb:"):
        return _assess_qecdb(target[len("qecdb:"):])
    if target.startswith("backend:"):
        return _assess_backend(target[len("backend:"):], pos[1] if len(pos) > 1 else "hamming_7")
    print(f"unknown target {target!r}; run 'q2sl assess' for the options")


def _assess_code(name, gamma=0.2, shots=40):
    import syndrome_leakage as sl
    from syndrome_leakage.analyze import population_leak, eavesdropper_error
    from syndrome_leakage import channels
    c = sl.codes.STANDARD[name]()
    r = sl.analyze(c, "amplitude_damping", gamma)
    w1 = sl.wasserstein.w1_leakage(c, "amplitude_damping", gamma)
    _, d0, d1 = population_leak(c, channels.LIBRARY["amplitude_damping"](gamma))
    eve = eavesdropper_error(d0, d1, shots)
    order = f"leaks at order gamma^{r.analytic_order}" if r.leaks else "protected"
    w1s = f"{w1['w1']:.3e}" + ("" if w1["exact"] else " (lower bound)")
    print(f"code: {c.name} under amplitude damping gamma={gamma}")
    print(f"  population {order}; phase {'protected' if r.phase < 1e-9 else 'leaks'}")
    print(f"  TVD {r.population:.3e} | W1 {w1s} | eavesdropper error over {shots} shots {eve:.3f}")
    part = _region_line(c, gamma) if r.leaks else None
    if part:
        print(f"  {part}")
    print(f"  read reliability from distinguishability: {0.5 + (0.5 - eve):.3f}")


def _assess_file(spec, gamma=0.2):
    from syndrome_leakage import load, protection
    from syndrome_leakage.analyze import analytic_leak
    paths = [p.strip() for p in spec.split(",")]
    try:
        code = load.load_code(*paths[:2]) if len(paths) > 1 else load.load_code(paths[0])
    except (OSError, AssertionError, ValueError) as e:
        print(f"could not read a code from {spec!r}: {e}")
        return 1
    return _report_code(code, gamma)


def _assess_qecdb(code_id, gamma=0.2):
    from syndrome_leakage import load
    try:
        code = load.from_qecdb(code_id.strip())
    except (LookupError, RuntimeError, ValueError, AssertionError) as e:
        print(f"could not load qecdb code {code_id!r}: {e}")
        return 1
    return _report_code(code, gamma)


def _qecdb(args=None):
    from syndrome_leakage import load
    query = {}
    for a in args or []:
        key, _, value = a.partition("=")
        if key in ("n", "k", "d") and value:
            lo, _, hi = value.partition("-")
            query[key] = (int(lo), int(hi or lo))
        elif key in ("family", "limit") and value:
            query[key] = int(value) if key == "limit" else value
        else:
            print("q2sl qecdb [n=9] [k=1] [d=3-5] [family=surface] [limit=25]")
            return 1
    try:
        rows = load.qecdb_search(**query)
    except (LookupError, RuntimeError) as e:
        print(e)
        return 1
    if not rows:
        print("no CSS code on qecdb.org matches")
        return 0
    for r in rows:
        print(f"  {r['_id']}  {r['name']:<14} {r.get('desc') or ''}")
    print(f"assess one with: q2sl assess qecdb:{rows[0]['_id']}")
    return 0


def _report_code(code, gamma=0.2):
    from syndrome_leakage import protection
    from syndrome_leakage.analyze import analytic_leak
    from syndrome_leakage.channels import amplitude_damping
    m = len(code.stab_strings)
    d = protection.distance(code)
    print(f"code: {code.name}, {code.n} qubits, {m} generators, k={code.k}, "
          f"distance {d if d is not None else 'beyond the search'}")
    for i, (z, x) in enumerate(code.logical_pairs):
        one = code if code.k == 1 else code.with_logical(i)
        label = "  " if code.k == 1 else f"  logical {i}: "
        print(f"{label}Z_L {z}")
        print(f"{' ' * len(label)}X_L {x}")
        try:
            leaks, order = analytic_leak(one)
            order_text = f"leaks at order gamma^{order}" if leaks else "protected"
        except AssertionError as e:
            order_text = _order_by_program(one, e)
        line = f"{' ' * len(label)}population {order_text}"
        if code.n <= 14:
            from syndrome_leakage import expectations as ex
            leak, _d0, _d1 = ex.population_leak(one, amplitude_damping(gamma))
            line += f", TVD {leak:.3e} at gamma={gamma}"
        print(line)
        if not order_text.startswith("protected"):
            part = _region_line(one, gamma)
            if part:
                print(f"{' ' * len(label)}{part}")
        if code.n > 14 and order_text != "protected":
            for text in _contracted(one, gamma, fit=not order_text.startswith("leaks")):
                print(f"{' ' * len(label)}{text}")
        if code.n <= 11:
            row = protection.protection_row(one, amplitude_damping(gamma))
            print(f"{' ' * len(label)}logical error after one round of recovery: "
                  + ", ".join(f"{k} {v:.3e}" for k, v in row.items()))
    if code.n > 14:
        print(f"  the exact leak needs {code.n} qubits simulated and the engine reaches about 14; past that,"
              f" the record of a few Z generators and half the L2 distance from the tensor network bound the leak from below")
    if code.k > 1:
        print("  the per-logical numbers belong to the basis printed above; another symplectic basis for "
              "the same code gives different ones")
    return 0


def _region_line(code, gamma, time_limit=120):
    # the fewest Z generators whose record leaks, and the distance of that record, which bounds the whole record's
    try:
        from syndrome_leakage import regions
        r = regions.smallest_leaking_set(code, time_limit=time_limit)
    except (ImportError, ValueError):
        return None
    if r["size"] is None:
        return None
    few = str(r["size"]) if r["optimal"] else f"at most {r['size']}"
    line = f"fewest Z generators that leak: {few}, generators {r['checks']}"
    try:
        leak = regions.z_region_leak(code, gamma, r["checks"])[0]
    except ValueError:
        return line
    return line + f"; TVD of their record {leak:.3e} at gamma={gamma}, at most the TVD of the whole record"


def _order_by_program(code, refusal, time_limit=120):
    # past the enumeration, the same order by integer programming, where scipy is installed
    try:
        from syndrome_leakage import zchecks
        r = zchecks.leak_order(code, time_limit=time_limit)
    except ImportError:
        return f"order out of reach: {refusal}; the integer program needs scipy"
    except ValueError as e:
        return f"order out of reach: {refusal}; {e}"
    if r["order"] is None:
        return "protected (integer program)"
    if r["optimal"]:
        return f"leaks at order gamma^{r['order']} (integer program)"
    return (f"leaks at order at most gamma^{r['order']}, at least gamma^{math.ceil(r['dual_bound'] - 1e-9)} "
            f"(integer program stopped at {time_limit} s)")


def _contracted(code, gamma, max_width=26, fit=False):
    # the L2 distance under amplitude damping, and with fit its exponent in gamma, which is the leak order
    import numpy as np
    from syndrome_leakage.channels import amplitude_damping
    try:
        from syndrome_leakage import tensor
        r = tensor.l2_leak(code, amplitude_damping(gamma), max_width=max_width)
    except ImportError:
        return ['L2 distance needs cotengra; pip install -e ".[tensor]"']
    except MemoryError as e:
        return [f"L2 distance out of reach: {e}"]
    except ValueError as e:
        return [f"L2 distance not computed: {e}"]
    lines = [f"L2 distance {r['l2']:.3e} at gamma={gamma}, so the leak is at least {r['lower']:.3e} "
             f"(tensor network, width {r['width']:.0f})"]
    if fit:
        gs = np.array([1e-4, 1e-5, 1e-6])
        l2 = [tensor.l2_leak(code, amplitude_damping(g), max_width=max_width)["l2"] for g in gs]
        if min(l2) > 0:
            slope = np.polyfit(np.log(gs), np.log(l2), 1)[0]
            text = f"L2 distance grows as gamma^{slope:.3f} over gamma 1e-4 to 1e-6"
            if slope >= 1 and abs(slope - round(slope)) < 0.05:
                text += f", so the leak order is {round(slope)}"
            lines.append(text)
        else:
            lines.append("L2 distance is zero to floating point at gamma 1e-4 to 1e-6")
    return lines


def _assess_backend(name, code="hamming_7"):
    from syndrome_leakage import hardware
    if not hardware.available():
        print('qiskit is not installed; pip install -e ".[hardware]"')
        return
    try:
        import qiskit_ibm_runtime  # noqa: F401
    except ImportError:
        print('backend lookup needs qiskit-ibm-runtime, which is not installed; '
              'pip install -e ".[hardware]"')
        return
    from syndrome_leakage.codes import STANDARD
    if code not in STANDARD:
        print(f"backend analysis takes one of " + ", ".join(STANDARD) + f"; not {code!r}")
        return
    backend = _load_backend(name)
    if backend is None:
        print(f"backend {name!r} not found in qiskit_ibm_runtime.fake_provider")
        return
    r = hardware.syndrome_leak_from_backend(backend, code=code)
    if r is None:
        print(f"backend {name!r} reports no T1")
        return
    print(f"backend: {r['backend']}  code: {r['code']}  qubits {r['qubits']}")
    print(f"  T1 {r['T1'] * 1e6:.1f} us, syndrome time {r['syndrome_time'] * 1e6:.3f} us, "
          f"gamma {r['gamma']:.3e}, syndrome TVD per shot {r['tvd']:.3e}")
    print(f"  Chernoff exponent {r['chernoff']:.3e} nats/shot, rounds to 1% error {r['rounds_for_1pct']}")
    for row in r["rows"]:
        print(f"  rounds {row['rounds']:>5}  attack error {row['attack_error']:.4f}  "
              f"Chernoff {row['chernoff_error']:.4f}  Bhattacharyya bound {row['bhattacharyya_bound']:.4f}")


def _load_backend(name):
    try:
        import qiskit_ibm_runtime.fake_provider as m
    except Exception:
        return None
    for cand in (name, "Fake" + name.capitalize() + "V2", "Fake" + name + "V2"):
        cls = getattr(m, cand, None)
        if cls is not None:
            try:
                return cls()
            except Exception:
                pass
    return None


def _leak(args=None):
    import syndrome_leakage as sl
    sl.selftest()


def _fusion(args=None):
    runpy.run_module("entropy_fusion", run_name="__main__")


def _eavesdrop(args=None):
    from syndrome_leakage import eavesdrop
    eavesdrop.main()


def _reconstruct(args=None):
    from reconstruction import scenario
    scenario.main()


def _showcase(args=None):
    path = os.path.join(HERE, "scenario", "one_posterior.py")
    if not os.path.exists(path):
        print("the end to end run lives in scenario/one_posterior.py, which ships with the repository "
              "and not with the installed package.")
        print("clone https://github.com/NikoKleion/q2sl-Attest and run it from there, or run "
              "`python scenario/one_posterior.py`.")
        return 1
    runpy.run_path(path, run_name="__main__")


def _suite(args=None):
    import suite
    from suite import report as R
    findings = suite.find()
    for f in findings:
        s, level, _ = R.severity(f)
        print(f"  {level:>7}  {f.target:>28}  {s:>4}  {f.headline[:52]}")
    sarif = R.to_sarif(findings)
    print(f"\n  SARIF 2.1.0: {len(sarif['runs'][0]['results'])} results")


def _all(args=None):
    for name, fn in (("leak", _leak), ("fusion", _fusion), ("suite", _suite)):
        print(f"\n===== {name} =====")
        fn()


DISPATCH = {"assess": _assess, "qecdb": _qecdb, "leak": _leak, "fusion": _fusion, "eavesdrop": _eavesdrop,
            "reconstruct": _reconstruct, "example": _showcase, "suite": _suite, "all": _all}


def main(argv=None):
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help", "help"):
        print("q2sl <command> [args]\n")
        for c, d in COMMANDS.items():
            print(f"  {c:>12}  {d}")
        return
    fn = DISPATCH.get(argv[0])
    if fn is None:
        print(f"unknown command {argv[0]!r}; run 'q2sl help'")
        return
    fn(argv[1:])


if __name__ == "__main__":
    main()
