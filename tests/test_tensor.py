# the tensor network L2 distance: equal to the exact engine wherever that engine reaches, exactly zero
# under a Pauli channel, and free of the exact engine's precision floor at small damping
import math

import numpy as np

from _harness import needs
from syndrome_leakage import expectations as ex
from syndrome_leakage.analyze import analytic_leak
from syndrome_leakage.analyze import l2_leak as exact_l2
from syndrome_leakage.channels import amplitude_damping, depolarizing
from syndrome_leakage.codes import STANDARD, shor_code
from syndrome_leakage.css import css_strings


def _tensor():
    needs("cotengra")
    from syndrome_leakage import tensor
    return tensor


def _toric(L):
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


def test_equals_the_exact_engine_on_every_reachable_code():
    tn = _tensor()
    toric = css_strings(*_toric(2), "toric L=2")
    codes = [STANDARD[n]() for n in ("repetition", "code_4_1_2", "hamming_7")]
    codes += [shor_code(), ex.surface_code_3(), toric.with_logical(0), toric.with_logical(1)]
    for code in codes:
        for gamma in (0.2, 0.05):
            K = amplitude_damping(gamma)
            exact = exact_l2(code, K)[0]
            got = tn.l2_leak(code, K)["l2"]
            assert abs(got - exact) <= 1e-12 * exact, (code.name, gamma, got, exact)


def test_a_protected_code_and_a_pauli_channel_give_zero():
    tn = _tensor()
    assert tn.l2_leak(STANDARD["steane"](), amplitude_damping(0.2))["l2"] < 1e-15
    for code in (shor_code(), ex.surface_code_3()):
        assert tn.l2_leak(code, depolarizing(0.1))["l2"] == 0.0


def test_the_lower_bound_is_below_the_leak():
    tn = _tensor()
    for code in (shor_code(), ex.surface_code_3()):
        K = amplitude_damping(0.1)
        assert tn.l2_leak(code, K)["lower"] <= ex.population_leak(code, K)[0]


def test_no_precision_floor_at_small_damping():
    # L2 / gamma^order must settle to a constant; the difference of two distributions loses it below
    # damping 1e-4, the contraction computes the difference as one object and does not
    tn = _tensor()
    code = STANDARD["hamming_7"]()
    w = analytic_leak(code)[1]
    ratios = [tn.l2_leak(code, amplitude_damping(g))["l2"] / g ** w for g in (1e-6, 1e-8, 1e-10)]
    assert max(ratios) - min(ratios) < 1e-5 * max(ratios), ratios
    assert math.isfinite(ratios[0]) and ratios[0] > 0


def test_a_code_that_is_not_css_is_refused():
    tn = _tensor()
    try:
        tn.l2_leak(STANDARD["five_qubit"](), amplitude_damping(0.1))
    except ValueError as e:
        assert "not CSS" in str(e)
    else:
        raise AssertionError("a non-CSS code was accepted")


def test_the_width_limit_refuses_rather_than_exhausts_memory():
    tn = _tensor()
    code = ex.surface_code_3()
    width = tn.contraction_width(code, amplitude_damping(0.1))
    try:
        tn.l2_leak(code, amplitude_damping(0.1), max_width=width - 1)
    except MemoryError as e:
        assert "past the limit" in str(e)
    else:
        raise AssertionError("the width limit was not enforced")


def test_a_channel_that_moves_x_parts_uses_the_general_network():
    # a coherent X rotation turns Z into Y, so the three copies keep separate X variables
    tn = _tensor()
    from syndrome_leakage.tensor import _locks_x, local_table
    theta = 0.3
    K = [np.cos(theta) * np.eye(2) - 1j * np.sin(theta) * np.array([[0, 1], [1, 0]])]
    assert not _locks_x(local_table(K))
    assert _locks_x(local_table(amplitude_damping(0.2)))
    for code in (STANDARD["code_4_1_2"](), shor_code(), ex.surface_code_3()):
        exact = exact_l2(code, K)[0]
        assert abs(tn.l2_leak(code, K)["l2"] - exact) <= 1e-12 * max(exact, 1e-300), code.name


def test_the_plan_is_reproducible_and_leaves_the_global_generator_alone():
    # the search is seeded, so two searches on one network agree, and the caller's random stream is intact
    import random
    tn = _tensor()
    code = ex.surface_code_3()
    random.seed(7)
    expected = random.random()
    random.seed(7)
    trees = []
    for _ in range(2):
        tn._PLANS.clear()
        inputs, arrays, _scale = tn._network(code, amplitude_damping(0.1))
        trees.append(tn._plan(inputs, arrays, None))
    assert random.random() == expected
    assert trees[0].contraction_width() == trees[1].contraction_width()
    assert trees[0].contraction_cost() == trees[1].contraction_cost()


def test_the_command_line_reads_the_order_from_the_contraction():
    # past the order search, the report fits the exponent of the L2 distance; a protected code gives zero
    _tensor()
    import q2sl
    lines = q2sl._contracted(css_strings(*ex.rotated_surface_code(5), "surface d=5"), 0.2, fit=True)
    assert lines[1].endswith("so the leak order is 5"), lines
    assert q2sl._contracted(STANDARD["steane"](), 0.2, fit=True)[1].startswith("L2 distance is zero")
    assert q2sl._contracted(STANDARD["five_qubit"](), 0.2)[0].startswith("L2 distance not computed")
