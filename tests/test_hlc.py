# syndrome probabilities from the generator coefficients of Hu, Liang and Calderbank, against the exact engine
import math

import numpy as np

import syndrome_leakage as sl
from syndrome_leakage.channels import coherent_diagonal
from syndrome_leakage.hlc import generator_coefficients, syndrome_probabilities

STD = sl.codes.STANDARD
CSS = ("repetition", "code_4_1_2", "steane", "hamming_7")


def test_generator_coefficients_are_normalised():
    # Theorem 6 of Hu, Liang and Calderbank: the squared coefficients sum to 1
    for name in CSS:
        code = STD[name]()
        for theta in (0.1, 0.7):
            A = generator_coefficients(code, theta)
            assert abs(sum(abs(a) ** 2 for a in A.values()) - 1.0) < 1e-12, name


def test_eq_91_matches_exact_simulation():
    for name in CSS:
        code = STD[name]()
        keys = list(code.PROJ.keys())
        for theta in (0.1, 0.3, 0.7):
            for t in (0.0, math.pi, math.pi / 2, 1.1):
                exact = code.syndrome_dist(code.apply(code.logical_state(t, 0.4), coherent_diagonal(theta)))
                p = syndrome_probabilities(code, theta, math.cos(t))
                closed = np.array([p.get(k, 0.0) for k in keys])
                assert np.allclose(closed, exact, atol=1e-12), (name, theta, t)


def test_steane_cross_terms_vanish():
    code = STD["steane"]()
    A = generator_coefficients(code, 0.3)
    keys = {k for k, _g in A}
    assert all(abs((A.get((k, 0), 0) * np.conj(A.get((k, 1), 0))).real) < 1e-12 for k in keys)
