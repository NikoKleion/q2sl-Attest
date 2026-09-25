# the documented runs, one command each; --save writes to results/
import math
import sys

import numpy as np

from . import codes as _codes
from .analyze import analytic_leak, population_leak, selftest as _selftest
from .channels import (amplitude_damping, channel_from_t1t2, coherent_diagonal, dephasing,
                       depolarizing)
from .core import tvd
from .css import code_distance, css_strings
from .eavesdrop import (bhattacharyya_coefficient, chernoff_exponent, ml_attack, repeated_extraction,
                        rounds_for_error, state_grid, worst_case_pair)
from .kl_check import ad_population_distance

STD = _codes.STANDARD
MARRAKESH = (163.3e-6, 77.0e-6, 1e-6)  # median T1, median T2, syndrome cycle
FEZ = (133.19e-6, 4000)  # T1 and shots of results/hardware_ibm_fez.json


def _order(code):
    leaks, o = analytic_leak(code)
    return o if leaks else None


def selftest():
    # analytic leak order against the exact simulation, every standard code
    _selftest()


def attack(trials=200000):
    # achieved eavesdropper error against the two closed-form bounds
    c = STD["hamming_7"]()
    leak, d0, d1 = population_leak(c, amplitude_damping(0.2))
    C, s = chernoff_exponent(d0, d1)
    print(f"{c.name}, amplitude_damping(0.2), {trials} trials")
    print(f"syndrome TVD {leak:.4e}, Chernoff {C:.6f} nats at s*={s:.3f}, "
          f"Bhattacharyya coefficient {bhattacharyya_coefficient(d0, d1):.6f}")
    print()
    print("  rounds   ML attack error      Chernoff rate   Bhattacharyya bound")
    for r in (1, 10, 50, 200):
        e = ml_attack(d0, d1, r, trials=trials, seed=r)
        se = math.sqrt(max(e * (1 - e), 1e-12) / trials)
        print(f"  {r:6d}   {e:.4f} +/- {se:.4f}     {0.5 * math.exp(-C * r):13.4f}   "
              f"{0.5 * bhattacharyya_coefficient(d0, d1) ** r:19.4f}")
    print()
    print(f"one-shot closed form (1 - TVD)/2 = {0.5 * (1 - leak):.4f}")
    print("rounds to 1 percent error, by code and gamma:")
    for g in (0.1, 0.2, 0.3):
        row = []
        for name in ("repetition", "code_4_1_2", "hamming_7"):
            _l, a, b = population_leak(STD[name](), amplitude_damping(g))
            row.append(f"{name} {rounds_for_error(a, b, 0.01)}")
        print(f"  gamma={g}: " + ", ".join(row))


def leak_order(gamma=1e-4):
    # analytic leak order against the codeword amplitude-damping population distance
    from . import expectations as ex
    from .hardware import shor_code
    print("code                       analytic order   AD population distance   match")
    for name in STD:
        c = STD[name]()
        a, k = _order(c), ad_population_distance(c)
        print(f"  {c.name:24} {str(a):>14}   {str(k):>22}   {'yes' if a == k else 'NO'}")
    print()
    print(f"the coefficient of the leading term, leak / gamma^order at gamma={gamma:g}, and the number of")
    print("minimal sets the order comes from")
    print(f"  {'code':26} {'order':>5} {'coefficient':>12} {'minimal sets':>13}")
    codes = [STD[n]() for n in ("repetition", "code_4_1_2", "hamming_7")]
    codes += [shor_code(), ex.surface_code_3()]
    for c in codes:
        w = _order(c)
        leak = ex.population_leak(c, amplitude_damping(gamma))[0]
        print(f"  {c.name:26} {w:5d} {leak / gamma ** w:12.4f} {_minimal_sets(c, w):13d}")
    print()
    print("  the coefficient is not the count: the surface code has eight minimal sets and a coefficient")
    print("  of seven, and on other codes the two differ by a factor of two in both directions")


def _minimal_sets(code, w):
    # sets of weight w carrying the order: Z_A in <S, Z_L>, inside a Z-support, odd overlap with X_L
    from .css import css_matrices
    _Hx, Hz = css_matrices(code.stab_strings)
    rows = [int("".join("1" if b else "0" for b in r[::-1]), 2) for r in Hz if r.any()]
    zl = int("".join("1" if ch == "Z" else "0" for ch in code.zl_str[::-1]), 2)

    def span(rs):
        out = [0]
        for r in rs:
            out += [x ^ r for x in out]
        return out

    supports = set(span(rows))
    vmask = sum(1 << i for i, ch in enumerate(code.xl_str) if ch in "XY")
    return sum(1 for a in span(rows + [zl]) if bin(a).count("1") == w
               and bin(a & vmask).count("1") % 2 == 1 and any(a & ~s == 0 for s in supports))


def coherent(theta=0.3):
    # coherent diagonal noise: which codes leak, how it scales, Pauli controls
    print(f"coherent_diagonal(theta={theta}), population TVD by code")
    for name in STD:
        leak, _, _ = population_leak(STD[name](), coherent_diagonal(theta))
        print(f"  {name:12} {leak:.4e}")
    c = STD["code_4_1_2"]()
    print("\n[[4,1,2]] scaling in theta")
    for t in (0.05, 0.1, 0.2, 0.3):
        leak, _, _ = population_leak(c, coherent_diagonal(t))
        print(f"  theta={t:.2f}  TVD {leak:.4e}  TVD/theta^2 {leak / t ** 2:.3f}")
    print("\nPauli controls on [[4,1,2]] (must be 0)")
    for nm, K in (("dephasing(0.1)", dephasing(0.1)), ("depolarizing(0.1)", depolarizing(0.1))):
        leak, _, _ = population_leak(c, K)
        print(f"  {nm}: {leak:.2e}")
    leak, d0, d1 = population_leak(c, coherent_diagonal(0.2))
    C, _s = chernoff_exponent(d0, d1)
    print(f"\nattack at theta=0.2: TVD {leak:.4e}, Chernoff {C:.4f} nats, "
          f"rounds to 1 percent {rounds_for_error(d0, d1, 0.01)}")


def structure():
    # leak order against code distance, CSS codes
    from . import expectations as ex
    from .hardware import shor_code
    from .css import css_matrices
    r = 3
    H = np.array([[(col >> (r - 1 - b)) & 1 for b in range(r)] for col in range(1, 2 ** r)], np.uint8).T
    print("projector path")
    for c in (STD["code_4_1_2"](), STD["hamming_7"](), STD["steane"](), css_strings(H, H, "Hamming-CSS self-dual")):
        print(f"  {c.name:26} d={code_distance(*css_matrices(c.stab_strings))} order {_order(c)}")
    print("\nstructure alone, no state simulated: the rotated surface code family")
    from .analyze import analytic_leak
    from . import protection as pr
    print(f"  {'code':26} {'n':>4} {'generators':>11} {'distance':>9} {'order':>6}")
    for dd in (3, 5):
        sc = css_strings(*ex.rotated_surface_code(dd), name=f"rotated surface d={dd}")
        print(f"  {sc.name:26} {sc.n:4d} {len(sc.stab_strings):11d} {pr.distance(sc):9d} "
              f"{analytic_leak(sc)[1]:6d}")
    print("\nnine qubits, strings-only codes with the expectation engine for the measured slope")
    for c in (shor_code(), ex.surface_code_3()):
        pts = [(g, ex.population_leak(c, amplitude_damping(g))[0]) for g in (0.01, 0.02, 0.04, 0.08)]
        slope = np.mean([math.log(d2 / d1) / math.log(g2 / g1) for (g1, d1), (g2, d2) in zip(pts, pts[1:])])
        print(f"  {c.name:26} d={code_distance(*css_matrices(c.stab_strings))} order {_order(c)} "
              f"measured slope {slope:.2f}")


def device(T1=MARRAKESH[0], T2=MARRAKESH[1], gate=MARRAKESH[2]):
    # leak and rounds-to-error at measured device parameters
    K = channel_from_t1t2(T1, T2, gate)
    print(f"T1={T1 * 1e6:.1f}us T2={T2 * 1e6:.1f}us, syndrome cycle {gate * 1e6:.1f}us")
    print(f"{'code':12} {'pop TVD/round':>14} {'rounds->1%':>12}")
    for name in STD:
        leak, d0, d1 = population_leak(STD[name](), K)
        print(f"{name:12} {leak:14.3e} {str(rounds_for_error(d0, d1, 0.01)):>12}")
    lc, _, _ = population_leak(STD["repetition"](), channel_from_t1t2(1e9, T2, gate))
    print(f"control, T1 to infinity: repetition TVD {lc:.2e}")


def held_memory(rounds=300):
    # held, actively corrected memory against the re-preparation regime
    K = channel_from_t1t2(*MARRAKESH)
    print(f"held state with recovery, {rounds} rounds, device parameters")
    for name, o in (("repetition", 1), ("code_4_1_2", 2), ("hamming_7", 3)):
        res = repeated_extraction(STD[name](), kraus_override=K, rounds=rounds, correct=True)
        tot = sum(r["chernoff"] for r in res["rounds"])
        print(f"  {name:12} leak order {o}: partial sum {tot:.4e} nats at {rounds} rounds, "
              f"attack bound {0.5 * math.exp(-tot):.4f}")
    _l, d0, d1 = population_leak(STD["repetition"](), K)
    print(f"  re-preparation, repetition: rounds to 1 percent {rounds_for_error(d0, d1, 0.01)}")
    rp = repeated_extraction(STD["repetition"](), kraus_override=dephasing(0.05), rounds=20, correct=True)
    print(f"  Pauli control, max per-round TVD {max(r['tvd'] for r in rp['rounds']):.2e}")
    print("  sums are partial sums at the stated round count")


def worst_pair():
    # leak between the population pair against other logical pairs
    grid = state_grid(9, 8)
    print(f"Bloch grid of {len(grid)} logical states, Chernoff maximised over all pairs")
    for ch, param in (("amplitude_damping", 0.2), ("coherent_diagonal", 0.3)):
        print(f"{ch} ({param})")
        for name in ("repetition", "code_4_1_2", "hamming_7"):
            w = worst_case_pair(STD[name](), gamma=param, channel=ch, grid=grid)
            print(f"  {name:12} worst {w['chernoff']:.4f}  axis {w['axis_chernoff']:.4f}  "
                  f"axis is worst: {w['axis_is_worst']}")


def pauli_boundary():
    # syndrome state-dependence for Pauli against non-Pauli channels
    grid = state_grid(7, 6)
    print(f"max pairwise syndrome TVD over {len(grid)} logical states")
    print(f"{'code':12} {'depolarizing':>14} {'dephasing':>12} {'amp damping':>13} {'coherent diag':>15}")
    for name in ("repetition", "code_4_1_2", "hamming_7"):
        c = STD[name]()
        row = []
        for K in (depolarizing(0.1), dephasing(0.1), amplitude_damping(0.2), coherent_diagonal(0.3)):
            ds = [c.syndrome_dist(c.apply(c.logical_state(t, p), K)) for t, p in grid]
            row.append(max(tvd(ds[i], ds[j]) for i in range(len(ds)) for j in range(i + 1, len(ds))))
        print(f"{name:12} {row[0]:14.2e} {row[1]:12.2e} {row[2]:13.3e} {row[3]:15.3e}")


def audit():
    # the attack against its closed form, and the Chernoff grid
    print("one-shot attack against the closed form (1 - TVD)/2, 200000 trials")
    for name, K in (("repetition", amplitude_damping(0.2)), ("code_4_1_2", amplitude_damping(0.2)),
                    ("code_4_1_2 coherent", coherent_diagonal(0.3)), ("hamming_7", amplitude_damping(0.2))):
        c = STD[name.split()[0]]()
        leak, d0, d1 = population_leak(c, K)
        sim = ml_attack(d0, d1, 1, trials=200000, seed=1)
        se = math.sqrt(max(sim * (1 - sim), 1e-12) / 200000)
        print(f"  {name:20} theory {0.5 * (1 - leak):.5f}  simulated {sim:.5f} +/- {se:.5f}")
    c = STD["hamming_7"]()
    _l, d0, d1 = population_leak(c, amplitude_damping(0.2))
    print("Chernoff grid resolution")
    for g in (51, 201, 2001):
        C, s = chernoff_exponent(d0, d1, grid=g)
        print(f"  grid {g:5d}: C={C:.8f} s*={s:.4f}")


def kl_blocks():
    # Knill-Laflamme branch matrices (Leung et al. 1997) beside the syndrome-resolved blocks, amplitude damping
    import numpy as np
    from .analyze import analytic_leak
    from .approx_qec import branch_matrices, branch_spread, fitted_order, population_split, syndrome_blocks
    from .channels import amplitude_damping
    from .codes import STANDARD
    g = 0.1
    M = branch_matrices(STANDARD["code_4_1_2"](), amplitude_damping(g))[(0, 0, 0, 0)]
    ev = sorted(np.linalg.eigvalsh(M))
    ref = sorted([(1 - g) ** 2, 0.5 * (1 + (1 - g) ** 4)])
    print(f"[[4,1,2]] no-jump branch at gamma {g}, eigenvalues of P_C A^dagger A P_C")
    print(f"  this package          {ev[0]:.12f}  {ev[1]:.12f}")
    print(f"  Leung et al. eq. 39   {ref[0]:.12f}  {ref[1]:.12f}")
    print()
    print(f"{'code':12} {'analytic order':>15} {'syndrome slope':>15} {'branch spread slope':>20}")
    fmt = lambda v: "none" if v is None else f"{v:.2f}"
    for name in STANDARD:
        c = STANDARD[name]()
        _leaks, order = analytic_leak(c)
        syn = fitted_order(lambda x: population_split(syndrome_blocks(c, amplitude_damping(x))))
        br = fitted_order(lambda x: branch_spread(c, amplitude_damping(x))[0])
        print(f"{name:12} {str(order if order is not None else 'none'):>15} {fmt(syn):>15} {fmt(br):>20}")


def surface():
    # the rotated [[9,1,3]] surface code on the stabilizer-expectation engine, past the projector limit
    import time
    import numpy as np
    from . import expectations as ex
    from .analyze import analytic_leak
    from .channels import amplitude_damping, coherent_diagonal, depolarizing
    from .codes import STANDARD
    from .core import op
    worst = 0.0
    for name in STANDARD:
        c = STANDARD[name]()
        for K in (amplitude_damping(0.2), coherent_diagonal(0.3), depolarizing(0.1)):
            for t, phi in ((0.0, 0.0), (math.pi, 0.0), (math.pi / 2, math.pi / 3)):
                a = c.syndrome_dist(c.apply(c.logical_state(t, phi), K))
                worst = max(worst, float(np.abs(a - ex.syndrome_dist(c, K, t, phi)).max()))
    print(f"standard codes, max |projector engine - expectation engine| = {worst:.1e}")
    sc = ex.surface_code_3()
    _leaks, order = analytic_leak(sc)
    print(f"{sc.name}: n {sc.n}, {len(sc.stab_strings)} stabilizers, Z_L {sc.zl_str}, X_L {sc.xl_str}, "
          f"analytic leak order {order}")
    print()
    print(f"{'channel':22} {'population TVD':>15}")
    rows = []
    for label, K in (("amplitude damping 0.1", amplitude_damping(0.1)), ("amplitude damping 0.2", amplitude_damping(0.2)),
                     ("coherent diagonal 0.3", coherent_diagonal(0.3)), ("depolarizing 0.1", depolarizing(0.1))):
        leak, d0, d1 = ex.population_leak(sc, K)
        rows.append((label, K, d0, d1))
        print(f"{label:22} {leak:15.6e}")
    # control: evolve the density matrix and apply the projector factors one at a time, nothing stored
    label, K, d0, d1 = rows[1]
    n = sc.n
    G = [op(s) for s in sc.stab_strings]
    I = np.eye(2 ** n)
    diff = 0.0
    for t, d in ((0.0, d0), (math.pi, d1)):
        psi = ex.logical_state_vector(sc, t, 0.0)
        T = np.outer(psi, psi.conj()).reshape([2] * n + [2] * n)
        for q in range(n):
            out = np.zeros_like(T)
            for Km in K:
                A = np.moveaxis(np.tensordot(Km, T, axes=([1], [q])), 0, q)
                out += np.moveaxis(np.tensordot(Km.conj(), A, axes=([1], [n + q])), 0, n + q)
            T = out
        rho = T.reshape(2 ** n, 2 ** n)
        ctrl = []
        for s in range(2 ** len(G)):
            M = rho
            for j, g in enumerate(G):
                M = ((I + (-1) ** ((s >> j) & 1) * g) @ M) / 2
            ctrl.append(float(np.real(np.trace(M))))
        ctrl = np.clip(np.array(ctrl), 0, None)
        diff = max(diff, float(np.abs(d - ctrl / ctrl.sum()).max()))
    print(f"control at {label}: max difference {diff:.1e}")
    ref = ex.stabilizer_expectations(sc, depolarizing(0.1), 0.0, 0.0)
    spread = max(float(np.abs(ex.stabilizer_expectations(sc, depolarizing(0.1), t, phi) - ref).max())
                 for t, phi in ((math.pi, 0.0), (math.pi / 2, 0.0), (math.pi / 2, math.pi / 3)))
    print(f"depolarizing, spread of stabilizer expectations over logical states {spread:.1e}")


def _decoder_section(codes, trials=200, seed=0):
    # the single-qubit table against BP+OSD, on the states and then on syndromes alone
    from . import expectations as ex
    from . import protection as pr
    from .css import css_strings
    try:
        import ldpc  # noqa: F401
    except ImportError:
        print("the decoder comparison needs ldpc")
        return
    print("recovery from the built-in table against BP+OSD of ldpc, depolarizing 0.01")
    print(f"  {'code':26} {'d':>2} {'table':>11} {'BP+OSD':>11} {'ratio':>7}")
    for c in codes:
        if pr.distance(c) < 3:
            continue
        a = max(pr.protection_row(c, depolarizing(0.01)).values())
        try:
            dec = pr.css_decoder(c)
        except ValueError:
            print(f"  {c.name:26} {pr.distance(c):2d} {a:11.3e} {'-':>11} {'not CSS':>7}")
            continue
        b = max(pr.protection_row(c, depolarizing(0.01), decoder=dec).values())
        print(f"  {c.name:26} {pr.distance(c):2d} {a:11.3e} {b:11.3e} {b / a:7.3f}")
    print()
    print(f"random errors on the d=5 rotated surface code, {trials} of each weight, syndromes only")
    big = css_strings(*ex.rotated_surface_code(5), name="rotated surface d=5")
    tab, osd = pr.table_decoder(big), pr.css_decoder(big)
    rng = np.random.default_rng(seed)
    print(f"  {'weight':>6} {'decoder':>8} {'corrected':>10} {'logical':>8} {'unexplained':>12}")
    for w in (1, 2):
        errors = []
        for _ in range(trials):
            e = ["I"] * big.n
            for q, p in zip(rng.choice(big.n, w, replace=False), rng.choice(list("XYZ"), w)):
                e[q] = p
            errors.append("".join(e))
        for name, dec in (("table", tab), ("BP+OSD", osd)):
            tally = {"corrected": 0, "logical": 0, "detected": 0}
            for e in errors:
                tally[pr.outcome(big, e, dec(pr.syndrome_of(big, e)))] += 1
            print(f"  {w:6d} {name:>8} {tally['corrected']:10d} {tally['logical']:8d} "
                  f"{tally['detected']:12d}")
    print()
    print("  unexplained: the correction leaves a syndrome behind, which the table does when it holds no")
    print("  entry for what was measured")
    print()


def _slope(f, a, b):
    va, vb = f(a), f(b)
    return (va, math.log(vb / va) / math.log(b / a)) if va > 1e-14 else (va, None)


def protection():
    # logical error after one round of syndrome measurement and single-qubit recovery, beside the leak
    from . import expectations as ex
    from . import protection as pr
    from .core import op
    from .hardware import shor_code
    worst = 0.0
    for name in STD:
        c = STD[name]()
        table = pr.recovery_strings(c)
        for K in (amplitude_damping(0.2), depolarizing(0.1), coherent_diagonal(0.3)):
            for t, phi in ((0.0, 0.0), (math.pi, 0.0), (math.pi / 2, math.pi / 3)):
                psi = ex.logical_state_vector(c, t, phi)
                rho = c.apply(np.outer(psi, psi.conj()), K)
                out = np.zeros_like(rho)
                for b, P in c.PROJ.items():
                    R = op(table.get(b, "I" * c.n))
                    out += R @ P @ rho @ P @ R
                ref = 1.0 - float(np.real(np.vdot(psi, out @ psi)))
                worst = max(worst, abs(ref - pr.logical_error(c, K, t, phi)))
    print(f"standard codes, max |projector form - recovery-fidelity form| = {worst:.1e}")
    codes = [STD[n]() for n in ("repetition", "code_4_1_2", "hamming_7", "five_qubit", "steane")]
    codes += [shor_code(), ex.surface_code_3()]
    worst = 0.0
    for c in codes:
        if pr.distance(c) < 3:
            continue
        for t, phi in ((0.0, 0.0), (math.pi, 0.0), (math.pi / 2, math.pi / 3)):
            psi = ex.logical_state_vector(c, t, phi)
            for q in range(c.n):
                for p in "XYZ":
                    v = ex.apply_pauli_string(psi, "I" * q + p + "I" * (c.n - q - 1))
                    worst = max(worst, abs(1.0 - pr.recovered_fidelity(c, np.outer(v, v.conj()), psi)))
    print(f"distance-3 codes, every single-qubit Pauli on 0_L, 1_L and a general state: max infidelity {worst:.1e}")
    print()
    print("amplitude damping and depolarizing at strength 0.01; slope of log value against log strength, 0.01 to 0.02")
    print("logical error is the largest of 0_L, 1_L and +_L")
    print()
    print(f"  {'code':26} {'d':>2} {'leak AD':>10} {'slope':>6} {'logical AD':>11} {'slope':>6} "
          f"{'logical dep':>12} {'slope':>6}")
    for c in codes:
        leak, ls = _slope(lambda g: ex.population_leak(c, amplitude_damping(g))[0], 0.01, 0.02)
        lad, las = _slope(lambda g: max(pr.protection_row(c, amplitude_damping(g)).values()), 0.01, 0.02)
        ldp, lds = _slope(lambda p: max(pr.protection_row(c, depolarizing(p)).values()), 0.01, 0.02)
        f = lambda s: f"{s:6.2f}" if s is not None else "     -"
        print(f"  {c.name:26} {pr.distance(c):2d} {leak:10.3e} {f(ls)} {lad:11.3e} {f(las)} {ldp:12.3e} {f(lds)}")
    print()
    _decoder_section(codes)


def shor_circuits(T1=FEZ[0], delays=(20e-6, 50e-6, 100e-6), shot_counts=(FEZ[1], 4 * FEZ[1])):
    # Shor-code hardware circuits on the local Aer simulator, and the leak against the shot-noise floor
    from . import expectations as ex
    from . import hardware as hw
    from .protection import syndrome_of
    try:
        from qiskit import QuantumCircuit, QuantumRegister, transpile
        from qiskit.quantum_info import Statevector
        from qiskit_aer import AerSimulator
    except ImportError:
        print("requires qiskit and qiskit-aer")
        return
    code = hw.shor_code()
    worst = 0.0
    for prep, t in ((0, 0.0), (1, math.pi)):
        data = QuantumRegister(9, "d")
        qc = QuantumCircuit(data)
        qc.h(0)
        if prep:
            qc.z(0)
        hw.shor_encoder(qc, data)
        v = Statevector(qc).reverse_qargs().data
        worst = max(worst, abs(1.0 - abs(np.vdot(v, ex.logical_state_vector(code, t, 0.0)))))
    print(f"encoder output against the 0_L and 1_L state vectors: max 1 - |overlap| {worst:.1e}")
    cases = [None] + [(q, p) for q in range(9) for p in "XYZ"]
    circs, want = [], []
    for e in cases:
        cs, _ = hw.build_shor_circuits([0.0], prep_errors=None if e is None else {(0.0, 0): e, (0.0, 1): e})
        s = 0 if e is None else sum(b << j for j, b in enumerate(syndrome_of(code, "I" * e[0] + e[1] + "I" * (8 - e[0]))))
        circs += cs
        want += [s, s]
    sim = AerSimulator(method="statevector")
    res = sim.run(transpile(circs, sim), shots=200, seed_simulator=1).result()
    ok = sum(hw.syndrome_dist_from_counts(res.get_counts(i), 200, n_bits=8)[s] == 1.0 for i, s in enumerate(want))
    print(f"Aer, 200 shots each: {ok} of {len(want)} circuits (no error and 27 single-qubit Paulis, both preps) "
          f"give the predicted syndrome on every shot")
    cs, _ = hw.build_shor_circuits([0.0], checks="z")
    print(f"the six Z generators alone: {cs[0].num_qubits} qubits, "
          f"{sum(1 for i in cs[0].data if i.operation.num_qubits == 2)} two-qubit gates, against "
          f"{hw.build_shor_circuits([0.0])[0][0].num_qubits} and "
          f"{sum(1 for i in hw.build_shor_circuits([0.0])[0][0].data if i.operation.num_qubits == 2)}")
    print()
    print(f"population leak at T1 {T1 * 1e6:.2f} us against the null floor, the mean TVD of two samples of 0_L")
    print("statistics of one record: all eight generators, the six Z generators, the number that fired")
    print()
    head = f"  {'statistic':22} {'delay us':>8} {'gamma':>6} {'leak':>7} "
    print(head + " ".join(f"{'floor ' + str(n):>12}" for n in shot_counts))
    rng = np.random.default_rng(0)
    stats = (("full syndrome (256)", lambda d: d),
             ("Z generators (64)", lambda d: hw.marginal_dist(d, hw.SHOR_Z_BITS)),
             ("number that fired (7)", lambda d: hw.syndrome_weight_dist(hw.marginal_dist(d, hw.SHOR_Z_BITS))))
    for name, f in stats:
        for t in delays:
            g = 1.0 - math.exp(-t / T1)
            _leak, d0, d1 = ex.population_leak(code, amplitude_damping(g))
            a, b = f(d0), f(d1)
            floors = [np.mean([tvd(rng.multinomial(n, a) / n, rng.multinomial(n, a) / n) for _ in range(400)])
                      for n in shot_counts]
            print(f"  {name:22} {t * 1e6:8.1f} {g:6.3f} {tvd(a, b):7.4f} "
                  + " ".join(f"{fl:12.4f}" for fl in floors))
        print()
    print("  the X generators carry no population leak, so dropping them leaves the leak and lowers the floor")


def rounds(n_rounds=3000, checkpoints=(10, 100, 300, 1000, 3000)):
    # a held state through many rounds of extraction and recovery, as one map on the logical qubit
    from . import expectations as ex
    from . import protection as pr
    from . import rounds as rd
    from .eavesdrop import repeated_extraction
    from .hardware import shor_code
    K = channel_from_t1t2(*MARRAKESH)
    print("one round is a fixed channel on the logical qubit when the decoder explains every syndrome,")
    print("so the round count costs nothing once the map is built")
    print()
    print("control against the density-matrix path of eavesdrop.repeated_extraction, same recovery:")
    for name in ("repetition", "five_qubit"):
        c = STD[name]()
        dec = pr.table_decoder(c)
        old = repeated_extraction(c, kraus_override=K, rounds=6, correct=True)["rounds"]
        new = rd.iterate(c, K, rounds=6, decoder=dec)
        worst = max(abs(a["tvd"] - b["tvd"]) for a, b in zip(old, new))
        print(f"  {c.name:26} max difference over six rounds {worst:.1e}")
    c = STD["repetition"]()
    flat = max(r["tvd"] for r in rd.iterate(c, dephasing(0.05), rounds=20))
    print(f"  Pauli control, dephasing on the repetition code: max per-round distance {flat:.1e}")
    print()
    print(f"held state at the {MARRAKESH[0] * 1e6:.1f} us and {MARRAKESH[1] * 1e6:.1f} us device channel, "
          f"{n_rounds} rounds")
    print("cumulative is the sum of the per-round Chernoff exponents, in nats")
    print()
    head = f"  {'code':26} {'order':>5} {'round 1':>10} {'round 300':>10} " + " ".join(
        f"{'sum ' + str(c):>10}" for c in checkpoints)
    print("recovery is BP+OSD, which explains every syndrome on all five codes")
    print(head)
    codes = [STD["repetition"](), STD["code_4_1_2"](), STD["hamming_7"](), shor_code(), ex.surface_code_3()]
    for c in codes:
        rows = rd.iterate(c, K, rounds=n_rounds, decoder=pr.css_decoder(c))
        sums = [rows[min(cp, n_rounds) - 1]["cumulative"] for cp in checkpoints]
        print(f"  {c.name:26} {str(_order(c)):>5} {rows[0]['tvd']:10.3e} {rows[299]['tvd']:10.3e} "
              + " ".join(f"{s:10.3e}" for s in sums))
    print()
    print("the sum depends on the decoder, since a correction that damages the state ends the leak with it")
    print(f"  {'decoder':26} {'round 1':>10} {'round 300':>10} {'sum 300':>10} {'sum 3000':>10}")
    c = STD["repetition"]()
    for label, dec in (("single-qubit table", pr.table_decoder(c)), ("BP+OSD", pr.css_decoder(c)),
                       ("any coset representative", rd.coset_decoder(c))):
        rows = rd.iterate(c, K, rounds=n_rounds, decoder=dec)
        print(f"  {label:26} {rows[0]['tvd']:10.3e} {rows[299]['tvd']:10.3e} "
              f"{rows[299]['cumulative']:10.3e} {rows[-1]['cumulative']:10.3e}")
    print()
    print("  the order sets how fast the sum stops growing, at a fixed decoder")


def shor_device_model(shots=None):
    # the ibm_fez circuits under a model from its calibration, against the measurement
    import copy
    import json
    import os
    import pickle
    from . import estimate as es
    from . import gate_level as gl
    try:
        import qiskit_aer  # noqa: F401
    except ImportError:
        print("requires qiskit and qiskit-aer")
        return
    root = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
    R = json.load(open(os.path.join(root, "hardware_shor_ibm_fez.json")))
    meta, counts = R["meta"], R["counts"]
    job = json.load(open(os.path.join(root, "hardware_shor_ibm_fez_job.json")))
    circs = gl.load_circuits(os.path.join(root, "hardware_shor_ibm_fez_circuits.qpy"))
    twins = gl.load_circuits(os.path.join(root, "hardware_shor_ibm_fez_twins.qpy"),
                             os.path.join(root, "hardware_shor_ibm_fez_twins_keys.json"))
    calfile = json.load(open(os.path.join(root, "ibm_fez_calibration.json")))
    cal, cmeta = calfile["calibration"], calfile["meta"]
    shots = shots or meta["shots"]
    labels = [tuple(x) for x in meta["labels"]]
    delays = sorted({t for t, _p in labels})
    print(f"the {len(circs)} circuits of job {meta['job_id']} as executed on {meta['backend']}, "
          f"{circs[0].count_ops().get('cz', 0)} cz each")
    print(f"calibration read {cmeta['fetched_utc'][:10]}, last updated {cmeta['calibration_last_update'][:19]}; "
          f"the job ran {meta['submitted_utc'][:10]}")
    print("every gate, every scheduled idle and every readout carries its own qubit's numbers")
    print(f"{shots} shots per circuit, seeds spaced by more than the shot count")
    print()
    lays = job["layouts"]
    distinct = {}
    for lab, lay in zip(labels, lays):
        distinct.setdefault(tuple(sorted(lay["data_qubits"] + lay["ancilla_qubits"])), []).append(lab)
    print(f"the transpiler placed the circuits on {len(distinct)} layouts:")
    for q, labs in distinct.items():
        print("  " + ", ".join(f"{t * 1e6:.0f} us |{p}_L>" for t, p in labs))
    print("so at every nonzero delay the two logical states ran on different qubits, and the measured")
    print("distance holds a layout difference as well as a state difference")
    print()
    job_day = copy.deepcopy(cal)
    for i, q in enumerate(job["data_qubits"]):
        job_day["qubits"][str(q)]["t1"], job_day["qubits"][str(q)]["t2"] = job["t1"][i], job["t2"][i]
    for i, q in enumerate(job["ancilla_qubits"]):
        job_day["qubits"][str(q)]["readout"] = job["ancilla_readout_error"][i]
    dev = {lab: es.dist_from_counts(counts[f"{lab[0]}_{lab[1]}"], 64)[0] for lab in labels}
    allc = list(circs) + list(twins.values())
    keys = list(labels) + list(twins.keys())
    models = {}
    for name, c in (("today", cal), ("job day", job_day)):
        sim_counts = gl.run_round(allc, c, shots, gl.spaced_seeds(len(allc), shots))
        models[name] = {k: es.dist_from_counts(x, 64)[0] for k, x in zip(keys, sim_counts)}
        json.dump({"meta": {"shots": shots, "calibration": name},
                   "counts": {str(k): x for k, x in zip(keys, sim_counts)}},
                  open(os.path.join(root, f"hardware_shor_model_{name.replace(' ', '_')}.json"), "w"))
    print("the measured distance against the model of the same two circuits, and the model's split of it:")
    print("  state, the two logical states on one layout (the twin of |0_L> on the |1_L> qubits), and")
    print("  layout, one logical state on the two layouts")
    print()
    for name in ("today", "job day"):
        M = models[name]
        print(f"  calibration: {name}")
        print(f"  {'delay us':>8} {'measured':>9} {'model':>8} {'ratio':>7} {'state':>8} {'layout':>8} "
              f"{'trivial meas':>13} {'trivial model':>14}")
        for t in delays:
            m = tvd(dev[(t, 0)], dev[(t, 1)])
            g = tvd(M[(t, 0)], M[(t, 1)])
            if t == 0:
                print(f"  {t * 1e6:8.1f} {m:9.4f} {g:8.4f} {'-':>7} {'-':>8} {'-':>8} "
                      f"{dev[(t, 0)][0]:13.4f} {M[(t, 0)][0]:14.4f}")
                continue
            tw = M[f"{t}_0_on_1"]
            state = tvd(tw, M[(t, 1)])
            layout = tvd(M[(t, 0)], tw)
            print(f"  {t * 1e6:8.1f} {m:9.4f} {g:8.4f} {m / g:7.3f} {state:8.4f} {layout:8.4f} "
                  f"{dev[(t, 0)][0]:13.4f} {M[(t, 0)][0]:14.4f}")
        print()
    print("  per-syndrome agreement, max |P_model - P_measured| over the 64 outcomes and both states")
    print(f"  {'delay us':>8} {'today':>8} {'job day':>8}")
    for t in delays:
        w = [max(float(np.abs(models[name][(t, p)] - dev[(t, p)]).max()) for p in (0, 1))
             for name in ("today", "job day")]
        print(f"  {t * 1e6:8.1f} {w[0]:8.4f} {w[1]:8.4f}")


def readout(delays=(20e-6, 50e-6, 100e-6), rates=(0.005, 0.01, 0.02, 0.05)):
    # error on the syndrome bits themselves: each bit flips with probability q
    from . import expectations as ex
    from . import hardware as hw
    code = hw.shor_code()
    rep = STD["repetition"]()

    def brute(d, q):
        m = int(round(math.log2(len(d))))
        out = np.zeros_like(d)
        for s in range(len(d)):
            for s2 in range(len(d)):
                w = 1.0
                for j in range(m):
                    w *= q if ((s >> j) & 1) != ((s2 >> j) & 1) else 1 - q
                out[s2] += d[s] * w
        return out

    d = rep.syndrome_dist(rep.apply(rep.logical_state(math.pi, 0.0), amplitude_damping(0.3)))
    worst = max(float(np.abs(ex.apply_readout(d, q) - brute(d, q)).max()) for q in (0.0, 0.01, 0.05, 0.2))
    print("a flip on syndrome bit j multiplies the expectation of every group element that contains "
          "generator j by 1 - 2q")
    print(f"control against bit flips applied one outcome at a time: max difference {worst:.1e}")
    print()

    def zleak(g, q, extra=0.0):
        gt = 1.0 - (1.0 - g) * (1.0 - extra)
        qs = [0.0, 0.0] + [q] * 6
        a = ex.syndrome_dist(code, amplitude_damping(gt), 0.0, 0.0, readout=qs)
        b = ex.syndrome_dist(code, amplitude_damping(gt), math.pi, 0.0, readout=qs)
        return tvd(hw.marginal_dist(a, hw.SHOR_Z_BITS), hw.marginal_dist(b, hw.SHOR_Z_BITS))

    print("leak after readout error, over the leak without it")
    print(f"  {'code':26} {'bits':>5} {'delay us':>8} " + " ".join(f"{'q=' + str(q):>9}" for q in rates))
    for t in delays:
        g = 1.0 - math.exp(-t / FEZ[0])
        base = population_leak(rep, amplitude_damping(g))[0]
        row = [population_leak(rep, amplitude_damping(g), readout=q)[0] / base for q in rates]
        print(f"  {rep.name:26} {2:5d} {t * 1e6:8.1f} " + " ".join(f"{r:9.4f}" for r in row))
    for t in delays:
        g = 1.0 - math.exp(-t / FEZ[0])
        base = zleak(g, 0.0)
        print(f"  {'Shor, six Z generators':26} {6:5d} {t * 1e6:8.1f} "
              + " ".join(f"{zleak(g, q) / base:9.4f}" for q in rates))
    print()
    print("  the factor does not move with the delay, which is the constant attenuation the hardware fit "
          "assumed")

    def implied(code_obj, target, bits=None):
        lo, hi = 0.0, 0.3
        g = 1.0 - math.exp(-50e-6 / FEZ[0])
        base = zleak(g, 0.0) if bits else population_leak(code_obj, amplitude_damping(g))[0]
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            r = (zleak(g, mid) if bits else population_leak(code_obj, amplitude_damping(g), readout=mid)[0]) / base
            lo, hi = (mid, hi) if r > target else (lo, mid)
        return 0.5 * (lo + hi)

    q_hw = implied(rep, 0.945)
    print()
    print(f"the hardware fit needed a constant attenuation of 0.945, which on two bits is q = {q_hw:.4f}")
    try:
        from qiskit_ibm_runtime.fake_provider import FakeFez
        target = FakeFez().target
        rd = [target["measure"][(q,)].error for q in range(156)]
        print(f"readout error of fake_fez: median {np.median(rd):.4f}, quartiles {np.percentile(rd, 25):.4f} "
              f"and {np.percentile(rd, 75):.4f}")
        print(f"so readout alone covers {np.median(rd) / q_hw * 100:.0f} percent of the flip rate the fit "
              f"implies")
    except ImportError:
        print("readout error of the backend needs qiskit-ibm-runtime")
    print()
    print("the `shor_noise` run measured 1.09, 0.89 and 0.77 of the ideal leak at the three delays, which")
    print("moves with the delay, so bit flips alone do not account for it. Extraction also damps the data")
    print("qubits, which raises the leak where it is small and compresses it where it is large:")
    print()
    print(f"  {'model':34} {'20 us':>8} {'50 us':>8} {'100 us':>8}")
    for label, extra, q in (("bit flips only, q=0.02", 0.0, 0.02),
                            ("extra damping 0.01, q=0.01", 0.01, 0.01),
                            ("extra damping 0.02, q=0.02", 0.02, 0.02),
                            ("extra damping 0.03, q=0.02", 0.03, 0.02)):
        cells = []
        for t in delays:
            g = 1.0 - math.exp(-t / FEZ[0])
            cells.append(zleak(g, q, extra) / zleak(g, 0.0))
        print(f"  {label:34} " + " ".join(f"{c:8.3f}" for c in cells))
    print(f"  {'shor_noise, full noise model':34} {1.09:8.2f} {0.89:8.2f} {0.77:8.2f}")


def disorder(delays=(20e-6, 50e-6, 100e-6), perms=200, seed=0):
    import itertools
    # one relaxation channel per qubit, from a device's own T1 and T2, against the uniform channel
    from . import expectations as ex
    from . import protection as pr
    from .channels import per_qubit_damping
    from .hardware import shor_code
    try:
        from qiskit_ibm_runtime.fake_provider import FakeFez
    except ImportError:
        print("requires qiskit-ibm-runtime")
        return
    target = FakeFez().target
    t1s = [target.qubit_properties[q].t1 for q in range(9)]
    all_t1 = [target.qubit_properties[q].t1 for q in range(156)]
    med = float(np.median(all_t1))
    print(f"T1 of the first nine qubits of fake_fez, us: {', '.join(f'{t * 1e6:.0f}' for t in t1s)}")
    print(f"median over all 156 qubits {med * 1e6:.1f} us, spread of the nine "
          f"{min(t1s) * 1e6:.0f} to {max(t1s) * 1e6:.0f} us")
    print()
    print("the uniform rows hold one channel on every qubit; the disordered rows hold one per qubit")
    print("the matched row is uniform at the mean of the nine gammas, so it differs only by the spread")
    print()
    print(f"  {'code':26} {'delay us':>8} {'uniform at median':>18} {'matched mean':>13} "
          f"{'disordered':>11} {'over matched':>13}")
    rng = np.random.default_rng(seed)
    for c in (shor_code(), ex.surface_code_3()):
        for t in delays:
            gam = [1.0 - math.exp(-t / x) for x in t1s]
            g_med = 1.0 - math.exp(-t / med)
            uni = ex.population_leak(c, amplitude_damping(g_med))[0]
            matched = ex.population_leak(c, amplitude_damping(float(np.mean(gam))))[0]
            dis = ex.population_leak(c, per_qubit_damping(gam))[0]
            print(f"  {c.name:26} {t * 1e6:8.1f} {uni:18.5f} {matched:13.5f} {dis:11.5f} "
                  f"{dis / matched:13.3f}")
        print()
    print(f"the same nine values dealt to the qubits in {perms} random orders, Shor at 50 us")
    c = shor_code()
    gam = [1.0 - math.exp(-50e-6 / x) for x in t1s]
    vals = []
    for _ in range(perms):
        vals.append(ex.population_leak(c, per_qubit_damping(list(rng.permutation(gam))))[0])
    vals = np.array(vals)
    matched = ex.population_leak(c, amplitude_damping(float(np.mean(gam))))[0]
    print(f"  min {vals.min():.5f}  median {np.median(vals):.5f}  max {vals.max():.5f}  "
          f"spread {vals.max() - vals.min():.5f}, {(vals.max() - vals.min()) / np.median(vals) * 100:.1f}"
          f" percent of the median")
    print(f"  uniform at the mean gamma {matched:.5f}, above {(vals < matched).mean() * 100:.0f} percent "
          f"of the orders")
    print()
    print("the measured code: three data qubits drawn from the device, against the median-T1 model")
    print(f"  {'delay us':>8} {'median model':>13} {'min':>8} {'median':>8} {'max':>8} {'spread/median':>14}")
    rep = STD["repetition"]()
    for t in delays:
        uni = population_leak(rep, amplitude_damping(1.0 - math.exp(-t / med)))[0]
        vals = []
        for _ in range(300):
            pick = [all_t1[q] for q in rng.choice(len(all_t1), 3, replace=False)]
            vals.append(population_leak(rep, per_qubit_damping([1.0 - math.exp(-t / x) for x in pick]))[0])
        vals = np.array(vals)
        print(f"  {t * 1e6:8.1f} {uni:13.4f} {vals.min():8.4f} {np.median(vals):8.4f} {vals.max():8.4f} "
              f"{(vals.max() - vals.min()) / np.median(vals):14.2f}")
    three = [1.0 - math.exp(-50e-6 / x) for x in t1s[:3]]
    orders = {round(population_leak(rep, per_qubit_damping(list(p)))[0], 10)
              for p in itertools.permutations(three)}
    print(f"  the same three values in all six orders give {len(orders)} value(s): "
          f"{', '.join(f'{v:.4f}' for v in sorted(orders))}")
    print()
    print("logical error after recovery, same nine qubits at 50 us")
    print(f"  {'code':26} {'matched mean':>13} {'disordered':>11} {'over matched':>13}")
    for c in (shor_code(), ex.surface_code_3()):
        gam = [1.0 - math.exp(-50e-6 / x) for x in t1s]
        a = max(pr.protection_row(c, amplitude_damping(float(np.mean(gam)))).values())
        b = max(pr.protection_row(c, per_qubit_damping(gam)).values())
        print(f"  {c.name:26} {a:13.5f} {b:11.5f} {b / a:13.3f}")


def estimator(reps=60, delays=(20e-6, 50e-6), shot_counts=(FEZ[1], 4 * FEZ[1])):
    # the finite-sample tools of estimate.py against known leaks: calibration, power, achieved distance
    from . import estimate as es
    from . import expectations as ex
    from . import hardware as hw
    code = hw.shor_code()
    z = lambda d: hw.marginal_dist(d, hw.SHOR_Z_BITS)
    print("statistics of one record against a known leak, Shor code at the ibm_fez T1")
    print("power is the fraction of draws with a permutation p below 0.05, at the stated leak")
    print("false positive is the same test on two draws of one state, which must sit near 0.05")
    print()
    print(f"  {'delay us':>8} {'shots':>7} {'statistic':22} {'leak':>7} {'power':>6} "
          f"{'false positive':>15} {'corrected':>10}")
    for t in delays:
        g = 1.0 - math.exp(-t / FEZ[0])
        _l, d0, d1 = ex.population_leak(code, amplitude_damping(g))
        stats = (("full syndrome (256)", d0, d1, None),
                 ("Z generators (64)", z(d0), z(d1), None),
                 ("number that fired (7)", z(d0), z(d1), hw.syndrome_weight_dist))
        for shots in shot_counts:
            for name, a, b, tr in stats:
                rng = np.random.default_rng(11)
                hits = fps = 0
                ach = []
                for i in range(reps):
                    x, y = rng.multinomial(shots, a) / shots, rng.multinomial(shots, b) / shots
                    r = es.leak_from_dists(x, shots, y, shots, transform=tr, boots=200, splits=4, seed=i,
                                           null_reps=8)
                    hits += r["p_value"] < 0.05
                    ach.append(r["distance_corrected"])
                    x2, y2 = rng.multinomial(shots, a) / shots, rng.multinomial(shots, a) / shots
                    r2 = es.leak_from_dists(x2, shots, y2, shots, transform=tr, boots=200, splits=4,
                                            seed=100 + i, null_reps=0)
                    fps += r2["p_value"] < 0.05
                leak = tvd(tr(a) if tr else a, tr(b) if tr else b)
                print(f"  {t * 1e6:8.1f} {shots:7d} {name:22} {leak:7.4f} {hits / reps:6.2f} "
                      f"{fps / reps:15.2f} {np.mean(ach):9.4f}")
        print()
    print("the achieved distance sits below the leak when the shots are few, since the rule is fitted on"
          " half of them")


def shor_noise(delays=(0.0, 20e-6, 50e-6, 100e-6), shot_counts=(FEZ[1], 4 * FEZ[1]), backend_name="fez"):
    # the Shor Z-check circuits under a fake backend's calibration, the rehearsal for a device run
    from . import expectations as ex
    from . import hardware as hw
    try:
        from qiskit import QuantumCircuit, transpile
        from qiskit_aer import AerSimulator
        from qiskit_aer.noise import NoiseModel, ReadoutError, depolarizing_error, thermal_relaxation_error
        from qiskit_ibm_runtime.fake_provider import FakeFez
    except ImportError:
        print("requires qiskit, qiskit-aer and qiskit-ibm-runtime")
        return
    backend = FakeFez()
    target = backend.target

    def compress(tqc):
        used = sorted({tqc.find_bit(q).index for inst in tqc.data for q in inst.qubits})
        pos = {q: i for i, q in enumerate(used)}
        out = QuantumCircuit(len(used), tqc.num_clbits)
        for inst in tqc.data:
            out.append(inst.operation, [pos[tqc.find_bit(q).index] for q in inst.qubits],
                       [tqc.find_bit(c).index for c in inst.clbits])
        return out, used

    def noise_for(used, delay_s):
        nm = NoiseModel(basis_gates=["rz", "sx", "x", "cz", "delay", "id"])
        pos = {q: i for i, q in enumerate(used)}
        for q in used:
            p = target.qubit_properties[q]
            t2 = min(p.t2, 2 * p.t1)
            for g in ("sx", "x"):
                props = target[g].get((q,))
                if props is None:
                    continue
                err = depolarizing_error(max(props.error or 0.0, 0.0), 1)
                nm.add_quantum_error(err.compose(thermal_relaxation_error(p.t1, t2, props.duration or 32e-9)),
                                     g, [pos[q]])
            m = target["measure"].get((q,))
            if m is not None and m.error:
                nm.add_readout_error(ReadoutError([[1 - m.error, m.error], [m.error, 1 - m.error]]), [pos[q]])
            if delay_s > 0:
                nm.add_quantum_error(thermal_relaxation_error(p.t1, t2, delay_s), "delay", [pos[q]])
        for (a, b), props in target["cz"].items():
            if a in pos and b in pos:
                dur = props.duration or 1e-7
                pa, pb = target.qubit_properties[a], target.qubit_properties[b]
                err = depolarizing_error(max(props.error or 0.0, 0.0), 2)
                relax = thermal_relaxation_error(pa.t1, min(pa.t2, 2 * pa.t1), dur).expand(
                    thermal_relaxation_error(pb.t1, min(pb.t2, 2 * pb.t1), dur))
                nm.add_quantum_error(err.compose(relax), "cz", [pos[a], pos[b]])
        return nm

    circs, labels = hw.build_shor_circuits(delays, checks="z")
    tqcs = transpile(circs, backend=backend, optimization_level=2, seed_transpiler=7)
    print(f"Shor Z-check circuits on the calibration of {backend.name}, one round, "
          f"{tqcs[0].num_qubits}-qubit device")
    print(f"transpiled depth {tqcs[0].depth()}, {sum(1 for i in tqcs[0].data if i.operation.name == 'cz')} cz, "
          f"on {len(compress(tqcs[0])[1])} qubits")
    print("|1_L> differs from |0_L> by a Z on one data qubit, a frame change, so the two preparations "
          "carry the same gate count")
    print("seeds are spaced by more than the shot count; close seeds share Aer's random numbers")
    print()
    dists = {}
    for shots in shot_counts:
        for k, ((t, prep), tqc) in enumerate(zip(labels, tqcs)):
            small, used = compress(tqc)
            sim = AerSimulator(method="statevector", noise_model=noise_for(used, t))
            for r in ((0, 1) if prep == 0 else (0,)):
                seed = 1_000_003 * (3 * k + r + 1) + (0 if shots == shot_counts[0] else 680_002_040)
                res = sim.run(small, shots=shots, seed_simulator=seed).result()
                dists[(shots, t, prep, r)] = hw.syndrome_dist_from_counts(res.get_counts(), shots, n_bits=6)
    code = hw.shor_code()
    for shots in shot_counts:
        print(f"{shots} shots per circuit")
        print(f"  {'delay us':>8} {'leak':>7} {'null':>7} {'ideal':>7} {'measured/ideal':>15} "
              f"{'fired leak':>11} {'fired null':>11} {'trivial 0_L':>12}")
        for t in delays:
            d0, d1, d0b = dists[(shots, t, 0, 0)], dists[(shots, t, 1, 0)], dists[(shots, t, 0, 1)]
            g = 1.0 - math.exp(-t / FEZ[0])
            _l, e0, e1 = ex.population_leak(code, amplitude_damping(g))
            ideal = tvd(hw.marginal_dist(e0, hw.SHOR_Z_BITS), hw.marginal_dist(e1, hw.SHOR_Z_BITS))
            leak = tvd(d0, d1)
            ratio = f"{leak / ideal:15.3f}" if ideal > 1e-9 else f"{'-':>15}"
            w0, w1, w0b = (hw.syndrome_weight_dist(x) for x in (d0, d1, d0b))
            print(f"  {t * 1e6:8.1f} {leak:7.4f} {tvd(d0, d0b):7.4f} {ideal:7.4f} {ratio} "
                  f"{tvd(w0, w1):11.4f} {tvd(w0, w0b):11.4f} {d0[0]:12.4f}")
        print()
    print("  null: two runs of |0_L> at that delay, the floor a device run measures at zero leak")
    print("  fired: the number of generators that fired, a statistic of the same record")


def _tn_codes():
    # constructions used by the tensor run: toric, hypergraph product and bivariate bicycle
    def rep(n):
        H = np.zeros((n - 1, n), np.uint8)
        for i in range(n - 1):
            H[i, i] = H[i, i + 1] = 1
        return H

    def hgp(H1, H2):
        m1, n1 = H1.shape
        m2, n2 = H2.shape
        Hx = np.hstack([np.kron(H1, np.eye(n2, dtype=np.uint8)), np.kron(np.eye(m1, dtype=np.uint8), H2.T)])
        Hz = np.hstack([np.kron(np.eye(n1, dtype=np.uint8), H2), np.kron(H1.T, np.eye(m2, dtype=np.uint8))])
        return Hx % 2, Hz % 2

    def toric(L):
        n = 2 * L * L
        h = lambda i, j: (i % L) * L + (j % L)
        v = lambda i, j: L * L + (i % L) * L + (j % L)
        Hx = np.zeros((L * L, n), np.uint8)
        Hz = np.zeros((L * L, n), np.uint8)
        for i in range(L):
            for j in range(L):
                for q in (h(i, j), h(i, j - 1), v(i, j), v(i - 1, j)):
                    Hx[i * L + j, q] ^= 1
                for q in (h(i, j), h(i + 1, j), v(i, j), v(i, j + 1)):
                    Hz[i * L + j, q] ^= 1
        return Hx, Hz

    def bb(l, m, a, b):
        # Bravyi et al. 2024: x = S_l (x) I_m, y = I_l (x) S_m, Hx = [A | B], Hz = [B^T | A^T]
        S = lambda k: np.roll(np.eye(k, dtype=np.uint8), 1, axis=1)
        x, y = np.kron(S(l), np.eye(m, dtype=np.uint8)), np.kron(np.eye(l, dtype=np.uint8), S(m))
        mono = lambda v, p: np.eye(l * m, dtype=np.uint8) if v == "1" else \
            np.linalg.matrix_power(x if v == "x" else y, p) % 2
        A, B = (sum(mono(v, p) for v, p in t) % 2 for t in (a, b))
        return np.hstack([A, B]) % 2, np.hstack([B.T, A.T]) % 2

    return rep, hgp, toric, bb


def tensor_network():
    # the L2 syndrome distance by tensor network contraction
    import importlib.util
    if importlib.util.find_spec("cotengra") is None:
        print("requires cotengra")
        return
    from . import expectations as ex
    from . import tensor as tn
    from .analyze import l2_leak as exact_l2
    from .css import css_strings
    rep, hgp, toric, bb = _tn_codes()

    t2 = css_strings(*toric(2), "toric L=2")
    small = [STD["repetition"](), STD["code_4_1_2"](), STD["hamming_7"](), _codes.shor_code(),
             ex.surface_code_3(), t2.with_logical(0), t2.with_logical(1)]
    theta = 0.3
    xrot = [np.cos(theta) * np.eye(2) - 1j * np.sin(theta) * np.array([[0, 1], [1, 0]])]
    worst = 0.0
    print("against the exact engine, L2 distance at amplitude damping 0.05")
    print(f"  {'code':30} {'exact':>12} {'contracted':>12}")
    for code in small:
        for K in (amplitude_damping(0.2), amplitude_damping(0.05), xrot):
            e, t = exact_l2(code, K)[0], tn.l2_leak(code, K)["l2"]
            worst = max(worst, abs(t - e) / e if e > 1e-12 else abs(t - e))
        e, t = exact_l2(code, amplitude_damping(0.05))[0], tn.l2_leak(code, amplitude_damping(0.05))["l2"]
        print(f"  {code.name:30} {e:12.6e} {t:12.6e}")
    print(f"  every code at damping 0.2 and 0.05 and under a coherent X rotation of {theta} agrees to 1e-12: "
          f"{worst < 1e-12}")
    print(f"  a Pauli channel gives {tn.l2_leak(_codes.shor_code(), depolarizing(0.1))['l2']:.1f} on Shor")
    print()

    code = STD["hamming_7"]()
    print(f"precision at small damping: L2 / gamma^3 on {code.name}, which must settle to a constant")
    print(f"  {'gamma':>7} {'exact engine':>14} {'contraction':>14}")
    for g in (1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10):
        K = amplitude_damping(g)
        e = f"{exact_l2(code, K)[0] / g ** 3:14.6f}" if g >= 1e-6 else f"{'':>14}"
        print(f"  {g:7.0e} {e} {tn.l2_leak(code, K)['l2'] / g ** 3:14.6f}")
    print()

    reach = [css_strings(*ex.rotated_surface_code(d), f"rotated surface d={d}") for d in (5, 7, 9, 11)]
    reach += [css_strings(*hgp(rep(r), rep(r)), f"hypergraph product rep{r}") for r in (4, 5)]
    reach += [css_strings(*toric(3), "toric L=3").with_logical(0)]
    print("past the exact engine, which stops at about 14 qubits: L2 at damping 0.05, and the exponent of")
    print("L2 against damping fitted over 1e-4, 1e-5 and 1e-6")
    print(f"  {'code':28} {'n':>4} {'k':>3} {'L2 at 0.05':>12} {'exponent':>9}")
    gs = np.array([1e-4, 1e-5, 1e-6])
    for code in reach:
        l2 = tn.l2_leak(code, amplitude_damping(0.05))["l2"]
        fit = [tn.l2_leak(code, amplitude_damping(g))["l2"] for g in gs]
        slope = np.polyfit(np.log(gs), np.log(fit), 1)[0]
        print(f"  {code.name:28} {code.n:4d} {code.k:3d} {l2:12.4e} {slope:9.3f}")
    out = [css_strings(*ex.rotated_surface_code(13), "rotated surface d=13"),
           css_strings(*toric(4), "toric L=4").with_logical(0),
           css_strings(*bb(3, 3, [("1", 0), ("x", 1), ("y", 2)], [("1", 0), ("y", 1), ("x", 2)]),
                       "BB [[18,4,4]]").with_logical(0)]
    for code in out:
        wide = tn.contraction_width(code, amplitude_damping(0.05)) > 26
        print(f"  {code.name:28} {code.n:4d} {code.k:3d}   out of reach: width past 26: {wide}")


def sampled_relaxation():
    # syndrome records under T1 and T2 relaxation from a stabilizer simulator
    import importlib.util
    if importlib.util.find_spec("stim") is None:
        print("requires stim")
        return
    from . import expectations as ex
    from . import sampled as sp
    rep, hgp, toric, bb = _tn_codes()

    def choi(kraus):
        units = [np.array(u, complex).reshape(2, 2) for u in np.eye(4)]
        return np.array([sum(k @ u @ k.conj().T for k in kraus) for u in units])

    worst = 0.0
    settings = ((100e-6, 60e-6, 20e-6), (163.3e-6, 77.0e-6, 1e-6), (163.3e-6, 77.0e-6, 120e-6),
                (50e-6, 50e-6, 30e-6))
    for T1, T2, t in settings:
        K = sp.mixture_channel(*sp.relaxation_mixture(t, T1, T2))
        worst = max(worst, float(np.abs(choi(K) - choi(channel_from_t1t2(T1, T2, t))).max()))
    try:
        sp.relaxation_mixture(1e-6, 50e-6, 80e-6)
        refused = False
    except ValueError:
        refused = True
    print(f"reset and Z mixture against channel_from_t1t2 at {len(settings)} settings, largest Choi difference "
          f"below 1e-14: {worst < 1e-14}")
    print(f"T2 > T1 refused: {refused}")
    print()

    pr, pz = sp.relaxation_mixture(30e-6, 100e-6, 60e-6)
    K = sp.mixture_channel(pr, pz)
    shots = 200000
    print(f"records against the exact engine: T1 100 us, T2 60 us, idle 30 us, so reset {pr:.4f} and Z {pz:.4f};")
    print(f"{shots} shots per state; distances between distributions of the full syndrome")
    print(f"  {'code':28} {'exact leak':>10} {'sampled':>9} {'0_L off':>9} {'1_L off':>9} {'floor':>8} "
          f"{'p 0_L':>7} {'p 1_L':>7}")
    t2 = css_strings(*toric(2), "toric L=2")
    small = [STD["repetition"](), STD["code_4_1_2"](), STD["five_qubit"](), _codes.shor_code(),
             ex.surface_code_3(), t2.with_logical(1)]
    ps = []
    for code in small:
        r0, r1 = sp.sample_records(code, pr, pz, shots, seed=12)
        d0, d1 = ex.syndrome_dist(code, K, 0.0), ex.syndrome_dist(code, K, math.pi)
        c0 = np.bincount(sp.statistic(r0, code, "syndrome"), minlength=len(d0))
        c1 = np.bincount(sp.statistic(r1, code, "syndrome"), minlength=len(d1))
        e0, e1 = c0 / shots, c1 / shots
        f = max(sp.sample_floor(d0, shots), sp.sample_floor(d1, shots))
        p0, p1 = sp.fit_p_value(c0, d0, seed=1), sp.fit_p_value(c1, d1, seed=2)
        ps += [p0, p1]
        print(f"  {code.name:28} {tvd(d0, d1):10.4f} {tvd(e0, e1):9.4f} {tvd(e0, d0):9.4f} {tvd(e1, d1):9.4f} "
              f"{f:8.4f} {p0:7.3f} {p1:7.3f}")
    print("  p is the chance that a sample drawn from the exact distribution lies as far from it; the smallest")
    print(f"  of the {len(ps)} is {min(ps):.3f}, and {sum(p < 0.05 for p in ps)} fall below 0.05")
    print()

    T1, T2 = 163.3e-6, 77.0e-6
    shots = 20000
    print(f"past the exact engines: T1 163.3 us and T2 77.0 us, {shots} shots per state, read through the number")
    print("of Z generators that read -1; the exact leak where the state vector reaches, the contracted L2 / 2")
    print("where the tensor network reaches")
    print(f"  {'code':28} {'n':>4} {'k':>3} {'idle us':>7} {'reset':>7} {'tvd':>7} {'null':>7} {'p':>6} "
          f"{'achieved':>9} {'exact':>9} {'L2 / 2':>9}")
    big = [css_strings(*ex.rotated_surface_code(d), f"rotated surface d={d}") for d in (3, 5, 7)]
    big += [css_strings(*bb(6, 6, [("x", 3), ("y", 1), ("y", 2)], [("y", 3), ("x", 1), ("x", 2)]),
                        "BB [[72,12,6]]").with_logical(0),
            css_strings(*bb(12, 6, [("x", 3), ("y", 1), ("y", 2)], [("y", 3), ("x", 1), ("x", 2)]),
                        "BB [[144,12,12]]").with_logical(0)]
    contract = importlib.util.find_spec("cotengra") is not None
    if contract:
        from . import tensor as tn
    for code in big:
        for t in (20e-6, 60e-6, 120e-6):
            pr, pz = sp.relaxation_mixture(t, T1, T2)
            K = sp.mixture_channel(pr, pz)
            r0, r1 = sp.sample_records(code, pr, pz, shots, seed=21)
            out = sp.leak_from_records(r0, r1, code, "z_weight", boots=1000, splits=10, seed=0, null_reps=10)
            exact = f"{ex.population_leak(code, K)[0]:9.4f}" if code.n <= 14 else f"{'':>9}"
            lower = f"{'':>9}"
            if contract and code.k == 1:
                lower = f"{tn.l2_leak(code, K)['lower']:9.2e}"
            print(f"  {code.name:28} {code.n:4d} {code.k:3d} {t * 1e6:7.0f} {pr:7.4f} {out['tvd']:7.4f} "
                  f"{out['null']:7.4f} {out['p_value']:6.3f} {out['distance_corrected']:9.4f} {exact} {lower}")


def z_checks():
    # the Z checks under relaxation as a classical process, and the leak order by integer programming
    import importlib.util
    from . import expectations as ex
    from . import sampled as sp
    from . import zchecks as zc
    from .estimate import achieved_distance, permutation_test
    rep, hgp, toric, bb = _tn_codes()
    t2 = css_strings(*toric(2), "toric L=2")
    small = [STD["repetition"](), STD["code_4_1_2"](), STD["hamming_7"](), STD["steane"](), _codes.shor_code(),
             ex.surface_code_3(), t2.with_logical(0), t2.with_logical(1)]

    def z_marginal(code, kraus, bit):
        zpos = zc.css_parts(code)["zpos"]
        out = {}
        for s, w in enumerate(ex.syndrome_dist(code, kraus, math.pi * bit)):
            key = tuple((s >> j) & 1 for j in zpos)
            out[key] = out.get(key, 0.0) + w
        return out

    def gap(a, b):
        return max(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in set(a) | set(b))

    T1, t = 100e-6, 30e-6
    reset = 1 - math.exp(-t / T1)
    worst = max(gap(z_marginal(code, channel_from_t1t2(T1, T2, t), bit),
                    z_marginal(code, sp.mixture_channel(reset, 0.0), bit))
                for code in small for T2 in (60e-6, 100e-6, 150e-6, 200e-6) for bit in (0, 1))
    print("Z-check syndromes under relaxation at T1 100 us, idle 30 us, and T2 of 60, 100, 150 and 200 us,")
    print(f"against reset alone, on {len(small)} codes: largest difference below 1e-12: {worst < 1e-12}")
    print()

    print("the classical process s = Hz (x AND R) against the exact engine at reset 0.3, and 200000 sampled")
    print("shots against the enumeration of the process; p is the chance a sample of the process lies as far")
    print(f"  {'code':28} {'exact to 1e-12':>14} {'p 0_L':>7} {'p 1_L':>7}")
    for code in small:
        ok = all(gap(z_marginal(code, sp.mixture_channel(0.3, 0.0), b), zc.exact_distribution(code, 0.3, b))
                 < 1e-12 for b in (0, 1))
        s0, s1 = zc.sample_syndromes(code, 0.3, 200000, seed=4)
        ps = []
        for bit, smp in ((0, s0), (1, s1)):
            exact = zc.exact_distribution(code, 0.3, bit)
            keys = sorted(exact)
            index = {k: i for i, k in enumerate(keys)}
            counts = np.zeros(len(keys))
            for row, c in zip(*np.unique(smp, axis=0, return_counts=True)):
                counts[index[tuple(int(v) for v in row)]] += c
            ps.append(sp.fit_p_value(counts, np.array([exact[k] for k in keys]), seed=bit))
        print(f"  {code.name:28} {str(ok):>14} {ps[0]:7.3f} {ps[1]:7.3f}")
    bb72 = css_strings(*bb(6, 6, [("x", 3), ("y", 1), ("y", 2)], [("y", 3), ("x", 1), ("x", 2)]), "BB [[72,12,6]]")
    bb144 = css_strings(*bb(12, 6, [("x", 3), ("y", 1), ("y", 2)], [("y", 3), ("x", 1), ("x", 2)]),
                        "BB [[144,12,12]]")
    if importlib.util.find_spec("stim") is not None:
        code = bb72.with_logical(0)
        zmask = np.array(["Z" in g for g in code.stab_strings])
        f0, f1 = zc.sample_syndromes(code, 0.3, 20000, seed=5)
        r0, r1 = sp.sample_records(code, 0.3, 0.0, 20000, seed=6)
        pv = []
        for a, b in ((f0, r0[:, zmask]), (f1, r1[:, zmask])):
            m = a.shape[1] + 1
            pv.append(permutation_test(np.bincount(a.sum(1), minlength=m), np.bincount(b.sum(1), minlength=m),
                                       boots=1000, seed=len(pv))["p_value"])
        print(f"  against the Stim sampler on {code.name}, reset 0.3, 20000 shots each, number of Z checks that")
        print(f"  fire: permutation p {pv[0]:.3f} for 0_L and {pv[1]:.3f} for 1_L")
    print()

    print("leak order by integer programming against the enumeration of the Z-side span")
    print(f"  {'code':28} {'integer program':>15} {'enumeration':>11}")
    checks = small + [css_strings(*ex.rotated_surface_code(5), "rotated surface d=5"),
                      css_strings(*toric(3), "toric L=3").with_logical(0),
                      css_strings(*toric(4), "toric L=4").with_logical(0),
                      css_strings(*hgp(rep(4), rep(4)), "hypergraph product rep4"),
                      css_strings(*bb(3, 3, [("1", 0), ("x", 1), ("y", 2)], [("1", 0), ("y", 1), ("x", 2)]),
                                  "BB [[18,4,4]]").with_logical(0)]
    same = True
    for code in checks:
        got = zc.leak_order(code)["order"]
        leaks, order = analytic_leak(code)
        want = order if leaks else None
        same = same and got == want
        print(f"  {code.name:28} {str(got):>15} {str(want):>11}")
    print(f"  equal on every code: {same}")
    print()

    print("Z distance by integer programming, against distances known in closed form and published")
    for code, known in ((css_strings(*ex.rotated_surface_code(7), "rotated surface d=7"), 7),
                        (css_strings(*toric(4), "toric L=4"), 4),
                        (css_strings(*hgp(rep(5), rep(5)), "hypergraph product rep5"), 5), (bb72, 6)):
        r = zc.z_distance(code)
        print(f"  {code.name:28} {r['distance']:3d}  proven {r['optimal']}  known {known}")
    for code in (bb72, bb144):
        r = zc.leak_order(code.with_logical(0))
        print(f"  leak order of {code.name}, logical 0: {r['order']}, proven {r['optimal']}")
    print()

    if importlib.util.find_spec("ldpc") is not None:
        print("which statistic of the Z syndrome carries the leak: exact distances by enumeration, the decoded")
        print("parity from BP+OSD")
        print(f"  {'code':28} {'reset':>5} {'full':>8} {'z_weight':>9} {'parity':>8} {'both':>8}")
        for code in (_codes.shor_code(), ex.surface_code_3(), t2.with_logical(0)):
            for pr in (0.1, 0.3, 0.5):
                stats = [None, lambda s: s.sum(1),
                         lambda s, pr=pr, code=code: zc.decoded_parity(code, s, pr),
                         lambda s, pr=pr, code=code: 2 * s.sum(1) + zc.decoded_parity(code, s, pr)]
                vals = []
                for f in stats:
                    d0, d1 = (zc.exact_distribution(code, pr, b, statistic=f) for b in (0, 1))
                    vals.append(0.5 * sum(abs(d0.get(k, 0.0) - d1.get(k, 0.0)) for k in set(d0) | set(d1)))
                print(f"  {code.name:28} {pr:5.1f} {vals[0]:8.5f} {vals[1]:9.5f} {vals[2]:8.5f} {vals[3]:8.5f}")
        print()

    print("the leak read from records through the number of Z checks that fire, T1 163.3 us, T2 not entering")
    print(f"  {'code':28} {'order':>5} {'idle us':>7} {'reset':>7} {'shots':>9} {'tvd':>8} {'null':>8} "
          f"{'p':>6} {'achieved':>9}")
    d5 = css_strings(*ex.rotated_surface_code(5), "rotated surface d=5")
    d7 = css_strings(*ex.rotated_surface_code(7), "rotated surface d=7")
    runs = [(ex.surface_code_3(), 3, t_us, 4_000_000) for t_us in (20, 60, 120)]
    runs += [(d5, 5, t_us, 4_000_000) for t_us in (20, 60, 120)]
    runs += [(d7, 7, 120, 10_000_000), (bb72.with_logical(0), 6, 120, 40_000_000),
             (bb144.with_logical(0), 12, 120, 20_000_000)]
    for i, (code, order, t_us, shots) in enumerate(runs):
        pr = 1 - math.exp(-t_us / 163.3)
        c0, c1 = zc.sample_counts(code, pr, shots, "z_weight", seed=100 + i)
        pt = permutation_test(c0, c1, boots=400, seed=1)
        ad = achieved_distance(c0, c1, splits=8, seed=2, null_reps=8)
        print(f"  {code.name:28} {order:5d} {t_us:7d} {pr:7.4f} {shots:9d} {pt['tvd']:8.5f} {pt['null']:8.5f} "
              f"{pt['p_value']:6.3f} {ad['distance_corrected']:9.5f}")


RUNS = {"selftest": selftest, "attack": attack, "leak_order": leak_order, "coherent": coherent,
        "structure": structure, "device": device, "held_memory": held_memory,
        "worst_pair": worst_pair, "pauli_boundary": pauli_boundary, "audit": audit, "kl_blocks": kl_blocks,
        "surface": surface, "protection": protection, "shor_circuits": shor_circuits,
        "shor_noise": shor_noise, "estimator": estimator, "disorder": disorder, "readout": readout, "rounds": rounds,
        "shor_device_model": shor_device_model, "tensor": tensor_network,
        "sampled": sampled_relaxation, "zchecks": z_checks}


def _incomplete(text):
    # a run whose whole output is a missing-dependency notice produced nothing worth saving
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    return not lines or all(ln.startswith("requires ") for ln in lines)


def main(argv):
    names = [a for a in argv if not a.startswith("-")] or list(RUNS)
    save = "--save" in argv
    import contextlib
    import io
    import os
    out_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
    for n in names:
        if n not in RUNS:
            print(f"unknown run {n}; choose from {', '.join(RUNS)}")
            return 1
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            RUNS[n]()
        text = buf.getvalue()
        print(f"----- {n} -----")
        print(text)
        if save:
            if _incomplete(text):
                print(f"not saving {n}.txt: the run could not do its work, and writing it would "
                      f"overwrite the saved output with a stub")
                continue
            os.makedirs(out_dir, exist_ok=True)
            with open(os.path.join(out_dir, f"{n}.txt"), "w", encoding="utf-8") as f:
                f.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
