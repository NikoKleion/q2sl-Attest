# finite-sample tools: calibration under the null, recovery of a known leak, the plug-in bias
import math

import numpy as np

from syndrome_leakage import estimate as es
from syndrome_leakage import expectations as ex
from syndrome_leakage import hardware as hw
from syndrome_leakage.channels import amplitude_damping
from syndrome_leakage.core import tvd


def _pair(delay=50e-6):
    code = hw.shor_code()
    g = 1.0 - math.exp(-delay / 133.19e-6)
    _l, d0, d1 = ex.population_leak(code, amplitude_damping(g))
    return hw.marginal_dist(d0, hw.SHOR_Z_BITS), hw.marginal_dist(d1, hw.SHOR_Z_BITS)


def test_dist_from_counts():
    d, n = es.dist_from_counts({"000001": 3, "000010": 1}, 64)
    assert n == 4 and d[1] == 0.75 and d[2] == 0.25 and abs(d.sum() - 1) < 1e-12


def test_permutation_test_is_calibrated():
    p0, _p1 = _pair()
    rng = np.random.default_rng(0)
    ps = []
    for i in range(40):
        a = rng.multinomial(4000, p0) / 4000
        b = rng.multinomial(4000, p0) / 4000
        ps.append(es.leak_from_dists(a, 4000, b, 4000, boots=200, splits=2, seed=i)["p_value"])
    ps = np.array(ps)
    assert (ps < 0.05).mean() <= 0.15, (ps < 0.05).mean()
    assert 0.3 < ps.mean() < 0.7, ps.mean()


def test_zero_leak_reads_as_zero():
    p0, _p1 = _pair()
    rng = np.random.default_rng(1)
    raw, corrected = [], []
    for i in range(5):
        a = rng.multinomial(4000, p0) / 4000
        b = rng.multinomial(4000, p0) / 4000
        r = es.leak_from_dists(a, 4000, b, 4000, boots=200, splits=4, seed=2 + i)
        raw.append(r["tvd"])
        corrected.append(r["distance_corrected"])
    assert min(raw) > 0.02, raw                      # the plug-in distance never reads zero
    assert np.mean(corrected) < 0.02, corrected      # the corrected achieved distance does


def test_achieved_distance_recovers_a_known_leak():
    p0, p1 = _pair()
    true = tvd(p0, p1)
    rng = np.random.default_rng(3)
    a = rng.multinomial(64000, p0) / 64000
    b = rng.multinomial(64000, p1) / 64000
    r = es.leak_from_dists(a, 64000, b, 64000, boots=200, splits=8, seed=4)
    assert abs(r["distance"] - true) < 0.02, (r["distance"], true)
    assert abs(r["distance_corrected"] - true) < 0.02, (r["distance_corrected"], true)
    assert r["distance_ci"][0] <= true <= r["distance_ci"][1]
    assert r["p_value"] < 0.01


def test_null_floor_and_shots_for_leak():
    p0, p1 = _pair(20e-6)
    assert es.null_floor(p0, 16000, trials=60, seed=0) < es.null_floor(p0, 4000, trials=60, seed=0)
    n = es.shots_for_leak(p0, p1, factor=2.0, trials=60)
    assert n is not None and es.null_floor(p0, n, trials=60, seed=0) * 2 <= tvd(p0, p1)


def test_counts_in_the_shapes_results_arrive_in():
    a, n = es.dist_from_counts({"00 01": 3, "0x2": 1, 3: 4}, 4)
    assert n == 8 and list(a) == [0.0, 0.375, 0.125, 0.5]
    b, n = es.dist_from_counts(["01", "0x1", 1, "11"], 4)
    assert n == 4 and list(b) == [0.0, 0.75, 0.0, 0.25]
    c, n = es.dist_from_memory(np.array([0, 1, 1, 3]), 4)
    assert n == 4 and list(c) == [0.25, 0.5, 0.0, 0.25]
    assert es.outcome_index("0b101") == 5 and es.outcome_index("1_0 1") == 5
