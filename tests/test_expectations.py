# the stabilizer-expectation engine against the projector engine, and the rotated surface code
import math

import numpy as np

import syndrome_leakage as sl
from syndrome_leakage import expectations as ex
from syndrome_leakage.analyze import analytic_leak
from syndrome_leakage.channels import (amplitude_damping, coherent_diagonal, depolarizing,
                                       generalized_amplitude_damping)
from syndrome_leakage.css import code_distance

STD = sl.codes.STANDARD
STATES = [(0.0, 0.0), (math.pi, 0.0), (math.pi / 2, 0.0), (math.pi / 2, math.pi / 3), (1.1, 0.4)]


def test_pauli_product_phases():
    assert ex.pauli_product("X", "Y") == (1j, "Z")
    assert ex.pauli_product("Y", "X") == (-1j, "Z")
    assert ex.pauli_product("Z", "X") == (1j, "Y")
    assert ex.pauli_product("XZ", "XI") == (1.0, "IZ")
    assert ex.pauli_product("XX", "ZZ") == (-1.0, "YY")


def test_group_elements_of_a_stabilizer_group_are_hermitian():
    for name in ("steane", "five_qubit", "code_4_1_2"):
        code = STD[name]()
        for phase, _ in ex.group_elements(code.stab_strings, code.n):
            assert phase in (1.0, -1.0), (name, phase)


def test_logical_states_match_the_projector_engine():
    for name, maker in STD.items():
        code = maker()
        for t, phi in STATES:
            psi = ex.logical_state_vector(code, t, phi)
            rho = code.logical_state(t, phi)
            assert np.allclose(np.outer(psi, psi.conj()), rho, atol=1e-12), (name, t, phi)


def test_syndrome_dist_matches_the_projector_engine():
    channels = [amplitude_damping(0.2), coherent_diagonal(0.3),
                generalized_amplitude_damping(0.3, 0.8), depolarizing(0.1)]
    for name, maker in STD.items():
        code = maker()
        for K in channels:
            for t, phi in STATES:
                a = code.syndrome_dist(code.apply(code.logical_state(t, phi), K))
                b = ex.syndrome_dist(code, K, t, phi)
                assert np.abs(a - b).max() < 1e-12, (name, t, phi)


def test_walsh_hadamard_inverts():
    v = np.arange(8, dtype=float)
    assert np.allclose(ex.walsh_hadamard(ex.walsh_hadamard(v)) / 8, v)


def test_rotated_surface_code_checks():
    Hx, Hz = ex.rotated_surface_code(3)
    assert not np.any((Hx @ Hz.T) % 2)
    assert code_distance(Hx, Hz) == 3
    sc = ex.surface_code_3()
    assert sc.n == 9 and len(sc.stab_strings) == 8


def test_surface_code_leaks():
    sc = ex.surface_code_3()
    assert analytic_leak(sc) == (True, 3)
    leak_01, _, _ = ex.population_leak(sc, amplitude_damping(0.1))
    leak_02, _, _ = ex.population_leak(sc, amplitude_damping(0.2))
    assert math.isclose(leak_01, 5.317650e-03, rel_tol=1e-5)
    assert math.isclose(leak_02, 3.128320e-02, rel_tol=1e-5)
    # every Z check has even weight and Z_L odd weight: no transversal Z rotation leaks
    # (Hu, Liang and Calderbank 2022, example 1)
    coherent, _, _ = ex.population_leak(sc, coherent_diagonal(0.3))
    assert coherent < 1e-12


def test_surface_code_pauli_channel_is_state_independent():
    # Wagner, Kampermann, Bruss and Kliesch 2022, eq. 40, on a code past the projector limit
    sc = ex.surface_code_3()
    K = depolarizing(0.1)
    ref = ex.stabilizer_expectations(sc, K, 0.0, 0.0)
    for t, phi in STATES[1:]:
        assert np.abs(ex.stabilizer_expectations(sc, K, t, phi) - ref).max() < 1e-12


def test_conjugate_pair_carries_no_population_leak():
    # amplitude damping separates logical states along the Z_L axis only
    from syndrome_leakage.hardware import shor_code
    from syndrome_leakage.core import tvd
    for c in (STD["repetition"](), STD["steane"](), STD["hamming_7"](), shor_code(), ex.surface_code_3()):
        for g in (0.05, 0.2):
            K = amplitude_damping(g)
            plus = ex.syndrome_dist(c, K, math.pi / 2, 0.0)
            minus = ex.syndrome_dist(c, K, math.pi / 2, math.pi)
            assert tvd(plus, minus) < 1e-12, (c.name, g)


def test_one_channel_per_qubit_reproduces_the_uniform_path():
    from syndrome_leakage import protection as pr
    from syndrome_leakage.channels import per_qubit_damping
    from syndrome_leakage.hardware import shor_code
    for c in (STD["steane"](), shor_code(), ex.surface_code_3()):
        uniform = ex.syndrome_dist(c, amplitude_damping(0.2), 0.0, 0.0)
        listed = ex.syndrome_dist(c, per_qubit_damping([0.2] * c.n), 0.0, 0.0)
        assert np.abs(uniform - listed).max() < 1e-15, c.name
        a = pr.logical_error(c, amplitude_damping(0.05))
        b = pr.logical_error(c, per_qubit_damping([0.05] * c.n))
        assert abs(a - b) < 1e-15, c.name
    c = STD["steane"]()
    rho = c.logical_state(0.0, 0.0)
    assert np.abs(c.apply(rho, amplitude_damping(0.2))
                  - c.apply(rho, per_qubit_damping([0.2] * 7))).max() < 1e-15


def test_disorder_moves_the_leak_and_depends_on_placement():
    from syndrome_leakage.channels import per_qubit_damping
    from syndrome_leakage.hardware import shor_code
    c = shor_code()
    gam = [0.05, 0.35, 0.2, 0.2, 0.2, 0.2, 0.05, 0.35, 0.2]
    matched = ex.population_leak(c, amplitude_damping(float(np.mean(gam))))[0]
    assert abs(ex.population_leak(c, per_qubit_damping(gam))[0] - matched) > 1e-3
    rng = np.random.default_rng(0)
    vals = [ex.population_leak(c, per_qubit_damping(list(rng.permutation(gam))))[0] for _ in range(20)]
    assert max(vals) - min(vals) > 1e-3, (min(vals), max(vals))


def test_readout_transform_matches_bit_flips():
    from syndrome_leakage.core import tvd

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

    c = STD["repetition"]()
    d = c.syndrome_dist(c.apply(c.logical_state(math.pi, 0.0), amplitude_damping(0.3)))
    for q in (0.0, 0.01, 0.05, 0.2):
        assert np.abs(ex.apply_readout(d, q) - brute(d, q)).max() < 1e-14, q
    assert np.abs(ex.apply_readout(d, 0.0) - d).max() < 1e-14
    assert abs(ex.apply_readout(d, 0.5).max() - 1 / len(d)) < 1e-12


def test_readout_attenuates_the_leak_by_a_factor_flat_in_gamma():
    from syndrome_leakage.hardware import shor_code
    for c in (STD["repetition"](), shor_code()):
        ratios = []
        for g in (0.1, 0.3, 0.5):
            base = ex.population_leak(c, amplitude_damping(g))[0]
            ratios.append(ex.population_leak(c, amplitude_damping(g), readout=0.02)[0] / base)
        assert max(ratios) - min(ratios) < 1e-9, (c.name, ratios)
        assert 0.5 < ratios[0] < 1.0, ratios


def test_rotated_surface_code_at_any_odd_distance():
    from syndrome_leakage.analyze import analytic_leak
    from syndrome_leakage.css import _rank, code_distance, css_strings
    Hx3, Hz3 = ex.rotated_surface_code(3)
    old_z = sorted([(0, 1, 3, 4), (4, 5, 7, 8), (2, 5), (3, 6)])
    old_x = sorted([(1, 2, 4, 5), (3, 4, 6, 7), (0, 1), (7, 8)])
    sets = lambda H: sorted(tuple(np.flatnonzero(r).tolist()) for r in H)
    assert sets(Hx3) == old_x and sets(Hz3) == old_z
    for d in (3, 5):
        Hx, Hz = ex.rotated_surface_code(d)
        assert not np.any((Hx @ Hz.T) % 2)
        assert d * d - _rank(Hx) - _rank(Hz) == 1
        assert code_distance(Hx, Hz) == d
        assert analytic_leak(css_strings(Hx, Hz, f"d={d}"))[1] == d


def test_analytic_order_matches_the_published_values():
    from syndrome_leakage.analyze import analytic_leak
    from syndrome_leakage.hardware import shor_code
    want = {"repetition": 1, "code_4_1_2": 2, "hamming_7": 3, "five_qubit": None, "steane": None}
    for name, order in want.items():
        assert analytic_leak(STD[name]())[1] == order, name
    assert analytic_leak(shor_code())[1] == 3
    assert analytic_leak(ex.surface_code_3())[1] == 3
