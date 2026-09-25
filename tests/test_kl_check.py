# leak_order run: the analytic leak order against the codeword amplitude-damping distance
import syndrome_leakage as sl
from syndrome_leakage.analyze import analytic_leak
from syndrome_leakage.kl_check import ad_population_distance

# name: (analytic leak order, amplitude-damping population distance)
EXPECTED = {
    "repetition": (1, 1),
    "code_4_1_2": (2, 2),
    "hamming_7": (3, 3),
    "steane": (None, 3),
    "five_qubit": (None, 5),
}


def _analytic_order(code):
    leaks, order = analytic_leak(code)
    return order if leaks else None


def test_correspondence_table():
    for name, (exp_order, exp_kl) in EXPECTED.items():
        code = sl.codes.STANDARD[name]()
        assert _analytic_order(code) == exp_order, f"{name}: analytic order changed"
        assert ad_population_distance(code) == exp_kl, f"{name}: AD population distance changed"


def test_analytic_is_not_the_ad_distance():
    # the two measures diverge on the codes whose population difference is syndrome-invariant
    diverge = [name for name in sl.codes.STANDARD
               if _analytic_order(sl.codes.STANDARD[name]()) != ad_population_distance(sl.codes.STANDARD[name]())]
    assert "steane" in diverge and "five_qubit" in diverge
