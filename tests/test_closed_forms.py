# the exact engine against closed forms
import math

import numpy as np

import syndrome_leakage as sl
from syndrome_leakage.analyze import population_leak
from syndrome_leakage.channels import amplitude_damping, coherent_diagonal, generalized_amplitude_damping
from syndrome_leakage.eavesdrop import chernoff_exponent
from syndrome_leakage.hlc import syndrome_probabilities

STD = sl.codes.STANDARD


def test_code_4_1_2_coherent_rotation():
    # one logical state never flips XXXX: TVD = sin^2(2 theta), Chernoff = -2 ln cos(2 theta)
    for theta in (0.1, 0.2, 0.3):
        leak, d0, d1 = population_leak(STD["code_4_1_2"](), coherent_diagonal(theta))
        assert abs(leak - math.sin(2 * theta) ** 2) < 1e-12, theta
        C, _ = chernoff_exponent(d0, d1)
        assert math.isclose(C, -2 * math.log(math.cos(2 * theta)), rel_tol=1e-9), theta


def test_code_4_1_2_amplitude_damping():
    # TVD = 4 g^2 (1-g)^2; two syndromes carry 2 g^2 (1-g)^2 for one state only
    for gamma in (0.05, 0.1, 0.2, 0.4):
        leak, d0, d1 = population_leak(STD["code_4_1_2"](), amplitude_damping(gamma))
        assert abs(leak - 4 * gamma ** 2 * (1 - gamma) ** 2) < 1e-12, gamma
        C, _ = chernoff_exponent(d0, d1)
        assert math.isclose(C, -math.log(1 - 2 * gamma ** 2 * (1 - gamma) ** 2), rel_tol=1e-9), gamma


def test_chernoff_ignores_float_residue():
    C, s = chernoff_exponent(np.array([0.848, 0.152, 0.0]), np.array([1.0, 1e-16, 0.0]))
    assert math.isclose(C, -math.log(0.848), rel_tol=1e-12)
    assert s == 1.0


def test_hlc_eq_93_steane():
    # Hu, Liang and Calderbank eq. 93: trivial syndrome (7 cos 4t + 25)/32, every other (1 - cos 4t)/32
    for theta in (0.1, 0.3, math.pi / 4, 1.0):
        p = syndrome_probabilities(STD["steane"](), theta, 1.0)
        trivial = tuple([0] * len(next(iter(p))))
        assert abs(p[trivial] - (7 * math.cos(4 * theta) + 25) / 32) < 1e-12, theta
        for key, value in p.items():
            if key != trivial:
                assert abs(value - (1 - math.cos(4 * theta)) / 32) < 1e-12, (theta, key)


def test_generalized_amplitude_damping_at_one_half_is_pauli():
    for name in ("repetition", "code_4_1_2", "hamming_7", "steane", "five_qubit"):
        leak, _, _ = population_leak(STD[name](), generalized_amplitude_damping(0.3, 0.5))
        assert leak < 1e-12, name


def test_generalized_amplitude_damping_p_swap():
    # X on every qubit swaps p and 1 - p when it is the logical X and every Z check has even weight
    for name in ("repetition", "steane"):
        c = STD[name]()
        _, a0, a1 = population_leak(c, generalized_amplitude_damping(0.3, 0.8))
        _, b0, b1 = population_leak(c, generalized_amplitude_damping(0.3, 0.2))
        assert np.allclose(a0, b1, atol=1e-12) and np.allclose(a1, b0, atol=1e-12), name


def test_values_published_by_shen_and_zhong():
    # arXiv:2609.09334: the distance 3 rotated surface code at damping 1e-4, one round, 6.9981e-12; and eq. 2
    # at one round, 1 - (1 - gamma)^n - gamma^n for the repetition code
    from syndrome_leakage import expectations as ex
    assert abs(ex.population_leak(ex.surface_code_3(), amplitude_damping(1e-4))[0] - 6.9981e-12) < 0.5e-16
    for g in (0.2, 0.05, 1e-4):
        leak = population_leak(STD["repetition"](), amplitude_damping(g))[0]
        assert abs(leak - (1 - (1 - g) ** 3 - g ** 3)) < 1e-14, g
