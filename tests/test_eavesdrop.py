# Tests for the syndrome eavesdropper experiment: the achieved attack, its bounds, and the controls.
import numpy as np

from syndrome_leakage import eavesdrop as ev
from syndrome_leakage.codes import STANDARD


def test_identical_distributions_give_chance_error():
    d = np.array([0.25, 0.25, 0.25, 0.25])
    assert abs(ev.ml_attack(d, d, rounds=50, trials=2000) - 0.5) < 0.05
    assert ev.rounds_for_error(d, d) is None


def test_disjoint_distributions_are_perfectly_distinguished():
    d0 = np.array([1.0, 0.0]); d1 = np.array([0.0, 1.0])
    assert ev.ml_attack(d0, d1, rounds=1, trials=1000) == 0.0


def test_attack_error_decreases_with_rounds():
    res = ev.attack_point(STANDARD["hamming_7"](), gamma=0.2, rounds=(1, 10, 50, 200), trials=3000, seed=2)
    errs = [r["attack_error"] for r in res["rows"]]
    assert errs == sorted(errs, reverse=True)
    assert errs[0] > 0.45 and errs[-1] < 0.35


def test_attack_is_at_least_as_strong_as_the_bounds():
    res = ev.attack_point(STANDARD["hamming_7"](), gamma=0.2, rounds=(10, 50, 200), trials=3000, seed=3)
    for r in res["rows"]:
        assert r["attack_error"] <= r["bhattacharyya_bound"] + 0.02
        assert r["attack_error"] <= r["chernoff_error"] + 0.02


def test_bhattacharyya_bound_is_loose_at_practical_round_counts():
    res = ev.attack_point(STANDARD["hamming_7"](), gamma=0.2, rounds=(200,), trials=3000, seed=4)
    row = res["rows"][0]
    assert row["bhattacharyya_bound"] - row["attack_error"] > 0.05


def test_chernoff_exponent_matches_bhattacharyya_when_symmetric():
    d0 = np.array([0.7, 0.3]); d1 = np.array([0.3, 0.7])
    C, s = ev.chernoff_exponent(d0, d1)
    assert abs(s - 0.5) < 0.02
    assert abs(C - (-np.log(ev.bhattacharyya_coefficient(d0, d1)))) < 1e-6


def test_pauli_channel_leaks_nothing():
    err, t = ev.pauli_control(rounds=200, trials=3000)
    assert t < 1e-12
    assert abs(err - 0.5) < 0.05


def test_protected_codes_are_indistinguishable():
    for name in ("steane", "five_qubit"):
        res = ev.attack_point(STANDARD[name](), gamma=0.2, rounds=(200,), trials=2000, seed=5)
        assert res["tvd"] < 1e-12, f"{name} should not leak under amplitude damping"
        assert abs(res["rows"][0]["attack_error"] - 0.5) < 0.06
        assert res["rounds_for_1pct"] is None


def test_repetition_code_leaks_quickly():
    res = ev.attack_point(STANDARD["repetition"](), gamma=0.2, rounds=(10,), trials=2000, seed=6)
    assert res["tvd"] > 0.1
    assert res["rounds_for_1pct"] is not None and res["rounds_for_1pct"] < 50
    assert res["rows"][0]["attack_error"] < 0.05


def test_leakage_grows_with_noise():
    need = [ev.attack_point(STANDARD["hamming_7"](), gamma=g, rounds=(10,), trials=500, seed=7)["rounds_for_1pct"]
            for g in (0.1, 0.2, 0.3)]
    assert need[0] > need[1] > need[2]


def test_report_runs():
    text = ev.report(trials=500)
    assert "ML attack error" in text and "Chernoff" in text


def test_population_axis_is_the_worst_case_pair():
    for name in ("repetition", "code_4_1_2", "hamming_7"):
        w = ev.worst_case_pair(STANDARD[name](), gamma=0.2)
        assert w["axis_is_worst"], f"{name}: a non-axis pair leaks more than |0_L> vs |1_L>"
        assert abs(w["tvd"] - w["axis_tvd"]) < 1e-9


def test_worst_case_finds_nothing_on_protected_codes():
    for name in ("steane", "five_qubit"):
        w = ev.worst_case_pair(STANDARD[name](), gamma=0.2)
        assert w["tvd"] < 1e-12


def test_state_grid_covers_poles_once():
    g = ev.state_grid()
    poles = [s for s in g if s[0] == 0.0 or abs(s[0] - np.pi) < 1e-12]
    assert len(poles) == 2


def test_recovery_table_is_complete_and_trivial_on_zero_syndrome():
    c = STANDARD["hamming_7"]()
    tbl = ev.recovery_table(c)
    assert len(tbl) == len(c.PROJ)
    assert np.allclose(tbl[0], np.eye(c.dim))
    for R in tbl:
        assert np.allclose(R @ R.conj().T, np.eye(c.dim))


def test_repeated_extraction_pauli_noise_stays_flat_at_zero():
    c = STANDARD["hamming_7"]()
    for corrected in (False, True):
        res = ev.repeated_extraction(c, gamma=0.1, channel="depolarizing", rounds=5, correct=corrected)
        for r in res["rounds"]:
            assert r["tvd"] < 1e-12, f"Pauli noise must not make the syndrome state-dependent: {r}"


def test_repeated_extraction_reports_every_round():
    c = STANDARD["hamming_7"]()
    res = ev.repeated_extraction(c, gamma=0.2, rounds=3)
    assert [r["round"] for r in res["rounds"]] == [1, 2, 3]
    assert res["rounds"][0]["tvd"] > 0.0


def test_held_state_leak_decays_and_total_is_bounded():
    c = STANDARD["hamming_7"]()
    for corrected in (False, True):
        res = ev.repeated_extraction(c, gamma=0.2, rounds=25, correct=corrected)
        tv = [r["tvd"] for r in res["rounds"]]
        assert tv == sorted(tv, reverse=True), "per-round leak on a held state must not grow"
        assert tv[-1] < 0.01 * tv[0], "the leak must decay by orders of magnitude"
        total = sum(r["chernoff"] for r in res["rounds"])
        assert total < 0.01, f"cumulative distinguishing information must stay bounded, got {total}"


def test_repeated_extraction_is_deterministic():
    c = STANDARD["hamming_7"]()
    a = ev.repeated_extraction(c, gamma=0.2, rounds=4)
    b = ev.repeated_extraction(c, gamma=0.2, rounds=4)
    assert [r["tvd"] for r in a["rounds"]] == [r["tvd"] for r in b["rounds"]]
