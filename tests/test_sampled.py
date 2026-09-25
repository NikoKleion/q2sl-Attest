# syndrome records sampled under T1 and T2 relaxation: the mixture is the relaxation channel, the records
# follow the exact distribution, a Pauli channel leaves the two logical states identical, and a seed fixes
# the records
import math

import numpy as np

from _harness import needs
from syndrome_leakage import expectations as ex
from syndrome_leakage import sampled
from syndrome_leakage.channels import channel_from_t1t2
from syndrome_leakage.codes import STANDARD
from syndrome_leakage.core import tvd


def _choi(kraus):
    units = [np.array(u, complex).reshape(2, 2) for u in np.eye(4)]
    return np.array([sum(k @ u @ k.conj().T for k in kraus) for u in units])


def test_the_mixture_is_the_relaxation_channel():
    for T1, T2, t in ((100e-6, 60e-6, 20e-6), (163.3e-6, 77.0e-6, 1e-6), (50e-6, 50e-6, 30e-6)):
        K = sampled.mixture_channel(*sampled.relaxation_mixture(t, T1, T2))
        assert np.abs(_choi(K) - _choi(channel_from_t1t2(T1, T2, t))).max() < 1e-14
    try:
        sampled.relaxation_mixture(1e-6, 50e-6, 80e-6)
    except ValueError as e:
        assert "T2 > T1" in str(e)
    else:
        raise AssertionError("T2 > T1 was accepted")
    pr, pz = sampled.relaxation_mixture(1e-6, [50e-6, 80e-6], [40e-6, 80e-6])
    assert pr.shape == (2,) and pz[1] == 0.0


def test_the_records_follow_the_exact_distribution():
    needs("stim")
    pr, pz = sampled.relaxation_mixture(30e-6, 100e-6, 60e-6)
    K = sampled.mixture_channel(pr, pz)
    shots = 4000
    for code in (STANDARD["code_4_1_2"](), STANDARD["five_qubit"](), ex.surface_code_3()):
        r0, r1 = sampled.sample_records(code, pr, pz, shots, seed=3)
        for bit, r in ((0, r0), (1, r1)):
            exact = ex.syndrome_dist(code, K, math.pi * bit)
            counts = np.bincount(sampled.statistic(r, code, "syndrome"), minlength=len(exact))
            assert tvd(counts / shots, exact) < 2 * sampled.sample_floor(exact, shots), (code.name, bit)
            assert sampled.fit_p_value(counts, exact) > 0.001, (code.name, bit)


def test_a_pauli_channel_leaves_the_two_states_identical():
    # Z errors keep both logical states in eigenstates of every generator, with the same signs
    needs("stim")
    code = ex.surface_code_3()
    rng = np.random.default_rng(0)
    shots, n, m = 300, code.n, len(code.stab_strings)
    reset = np.zeros((shots, n), bool)
    z = rng.random((shots, n)) < 0.2
    coins = rng.random((shots, n + m)) < 0.5
    r0 = sampled._run(code, 0, reset, z, coins)
    assert r0.any()
    assert np.array_equal(r0, sampled._run(code, 1, reset, z, coins))


def test_a_seed_fixes_the_records():
    needs("stim")
    code = STANDARD["code_4_1_2"]()
    a = sampled.sample_records(code, 0.3, 0.1, 200, seed=5)
    b = sampled.sample_records(code, 0.3, 0.1, 200, seed=5)
    c = sampled.sample_records(code, 0.3, 0.1, 200, seed=6)
    assert all(np.array_equal(x, y) for x, y in zip(a, b))
    assert not np.array_equal(a[0], c[0])


def test_a_statistic_of_the_records_bounds_the_leak_from_below():
    needs("stim")
    code = ex.surface_code_3()
    pr, pz = sampled.relaxation_mixture(60e-6, 100e-6, 60e-6)
    leak = ex.population_leak(code, sampled.mixture_channel(pr, pz))[0]
    r0, r1 = sampled.sample_records(code, pr, pz, 6000, seed=1)
    out = sampled.leak_from_records(r0, r1, code, "z_weight", boots=300, splits=6, null_reps=6)
    assert out["p_value"] < 0.01
    assert out["distance_corrected"] <= leak + 0.02


def test_per_qubit_parameters_and_the_weight_statistics():
    needs("stim")
    code = STANDARD["code_4_1_2"]()
    r0, _r1 = sampled.sample_records(code, [0.0, 0.0, 0.5, 0.5], 0.0, 50, seed=0, readout=[0.0] * 3)
    assert r0.shape == (50, 3)
    try:
        sampled.sample_records(code, [0.1, 0.2], 0.0, 5)
    except ValueError as e:
        assert "p_reset" in str(e)
    else:
        raise AssertionError("a p_reset of the wrong length was accepted")
    rec = np.ones((2, 3), np.uint8)
    z = sum(1 for s in code.stab_strings if set(s) <= {"I", "Z"})
    assert list(sampled.statistic(rec, code, "z_weight")) == [z, z]
    assert list(sampled.statistic(rec, code, "weight")) == [3, 3]
