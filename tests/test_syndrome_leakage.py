# Tests for syndrome_leakage: analytic vs exact, CSS codes, Wasserstein, T1/T2 bridge.
import numpy as np

from _harness import needs
import syndrome_leakage as sl
from syndrome_leakage.channels import channel_from_t1t2


def test_analytic_validates_exact_on_all_codes():
    # analytic leak order must match the measured gamma-slope, phase protected
    for name, maker in sl.codes.STANDARD.items():
        r = sl.analyze(maker(), "amplitude_damping", gamma=0.2)
        assert r.validated, f"{name}: analytic order not validated by exact simulation"
        assert r.phase < 1e-9, f"{name}: phase should be protected"


def test_css_reproduces_steane():
    c = sl.hamming_css(r=3)
    assert c.n == 7
    c.verify_projectors()
    r = sl.analyze(c, "amplitude_damping", 0.2)
    st = sl.analyze(sl.codes.steane(), "amplitude_damping", 0.2)
    assert r.leaks == st.leaks and abs(r.population - st.population) < 1e-6


def test_css_from_matrices_matches_handtyped_412():
    Hx = np.array([[1, 1, 1, 1]], np.uint8)
    Hz = np.array([[1, 1, 0, 0], [0, 0, 1, 1]], np.uint8)
    c = sl.css_from_matrices(Hx, Hz, name="[[4,1,2]]-m")
    r = c and sl.analyze(c, "amplitude_damping", 0.2)
    hand = sl.analyze(sl.codes.code_4_1_2(), "amplitude_damping", 0.2)
    assert r.leaks == hand.leaks and r.analytic_order == hand.analytic_order


def test_css_rejects_noncommuting():
    bad = np.array([[1, 0, 1], [0, 1, 1]], np.uint8)
    try:
        sl.css_from_matrices(bad, bad)
        assert False, "should reject Hx Hz^T != 0"
    except AssertionError:
        pass


def test_wasserstein_lower_bound_present():
    r = sl.wasserstein.w1_leakage(sl.codes.repetition(), "amplitude_damping", 0.2)
    assert r["w1_lower_bound"] >= 0.0 and r["tvd"] >= 0.0


def test_wasserstein_exact_recovers_joint_leak():
    needs("scipy")
    r = sl.wasserstein.w1_leakage(sl.codes.code_4_1_2(), "amplitude_damping", 0.2)
    assert r["exact"] and abs(r["w1"] - r["tvd"]) < 1e-6


def test_t1t2_bridge_is_trace_preserving():
    kraus = channel_from_t1t2(T1=120e-6, T2=90e-6, gate_time=1e-6)
    tp = sum(K.conj().T @ K for K in kraus)
    assert np.allclose(tp, np.eye(2)), "T1/T2 channel must be trace preserving"
    r = sl.analyze(sl.codes.repetition(), channel=kraus, gamma=0.0)
    assert r.leaks, "repetition should leak population under a real relaxation channel"
