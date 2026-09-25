# structure run: leak order across codes built through the strings-only CSS path, past the projector ceiling
import numpy as np

import syndrome_leakage as sl
from syndrome_leakage.analyze import analytic_leak
from syndrome_leakage.css import css_strings, code_distance


def _order(code):
    leaks, o = analytic_leak(code)
    return o if leaks else None


def _hamming_H(r):
    ncol = 2 ** r - 1
    return np.array([[(col >> (r - 1 - bit)) & 1 for bit in range(r)] for col in range(1, ncol + 1)], np.uint8).T


def test_hamming_7_distance_and_order():
    from syndrome_leakage.css import css_matrices
    c = sl.codes.STANDARD["hamming_7"]()
    assert code_distance(*css_matrices(c.stab_strings)) == 2
    assert _order(c) == 3
    s = sl.codes.STANDARD["steane"]()
    assert code_distance(*css_matrices(s.stab_strings)) == 3
    assert _order(s) is None


def test_self_dual_hamming_css_is_protected():
    H = _hamming_H(3)
    sc = css_strings(H, H, "Hamming-CSS[[7,1,3]]")
    assert code_distance(H, H) == 3
    assert _order(sc) is None


def test_string_path_reaches_shor_and_it_leaks():
    Hx = np.array([[1, 1, 1, 1, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 1, 1, 1, 1]], np.uint8)
    Hz = np.array([[1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 1, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 1]], np.uint8)
    shor = css_strings(Hx, Hz, "Shor[[9,1,3]]")
    assert code_distance(Hx, Hz) == 3
    assert _order(shor) == 3


def test_the_leading_coefficient_is_not_the_count_of_minimal_sets():
    # pinned so a future derivation has targets, and so the coincidence is not mistaken for a rule
    from syndrome_leakage import expectations as ex
    from syndrome_leakage.channels import amplitude_damping
    from syndrome_leakage.experiments import _minimal_sets
    from syndrome_leakage.hardware import shor_code
    g = 1e-4
    want = {"3-qubit repetition": (1, 3.0, 3), "[[4,1,2]]": (2, 4.0, 4), "Hamming [[7,1,2]]": (3, 4.0, 4),
            "Shor [[9,1,3]]": (3, 27.0, 27), "rotated surface [[9,1,3]]": (3, 7.0, 8)}
    codes = [sl.codes.STANDARD[n]() for n in ("repetition", "code_4_1_2", "hamming_7")]
    codes += [shor_code(), ex.surface_code_3()]
    for c in codes:
        w, coeff, sets = want[c.name]
        assert _order(c) == w, c.name
        fitted = ex.population_leak(c, amplitude_damping(g))[0] / g ** w
        assert abs(fitted - coeff) < 0.01 * coeff, (c.name, fitted, coeff)
        assert _minimal_sets(c, w) == sets, c.name
    assert want["rotated surface [[9,1,3]]"][1] != want["rotated surface [[9,1,3]]"][2]


def test_l2_bounds_bracket_the_leak():
    # the L2 distance bounds the leak from both sides, on every code that leaks
    from syndrome_leakage.analyze import l2_leak
    from syndrome_leakage.channels import amplitude_damping
    from syndrome_leakage.codes import STANDARD, shor_code
    from syndrome_leakage.expectations import population_leak
    for code in [STANDARD[n]() for n in ("repetition", "code_4_1_2", "hamming_7")] + [shor_code()]:
        for gamma in (0.05, 1e-3):
            kraus = amplitude_damping(gamma)
            leak = population_leak(code, kraus)[0]
            l2, lo, hi, shape = l2_leak(code, kraus)
            assert lo <= leak <= hi * (1 + 1e-9), (code.name, gamma, leak, lo, hi)
            assert 1.0 <= shape <= 2 ** (code.n / 2), (code.name, gamma, shape)


def test_the_shape_factor_is_flat_in_the_damping_strength():
    # it is what lets a measurement at a convenient strength be reused at a smaller one, and it is
    # also the canary for the precision floor: it wanders once the leak is a difference of equals
    from syndrome_leakage.analyze import l2_leak
    from syndrome_leakage.channels import amplitude_damping
    from syndrome_leakage.codes import STANDARD
    code = STANDARD["hamming_7"]()
    shapes = [l2_leak(code, amplitude_damping(g))[3] for g in (0.05, 1e-2, 1e-3, 1e-4)]
    assert max(shapes) - min(shapes) < 0.05, shapes
