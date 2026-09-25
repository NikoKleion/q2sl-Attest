# the Z checks under relaxation as a classical process: equal to the exact engine, blind to T2, sampled
# faithfully, and the integer-program order equal to the enumeration wherever the enumeration reaches
import math

import numpy as np

from _harness import needs
from syndrome_leakage import expectations as ex
from syndrome_leakage import sampled, zchecks
from syndrome_leakage.analyze import analytic_leak
from syndrome_leakage.channels import channel_from_t1t2
from syndrome_leakage.codes import STANDARD, shor_code
from syndrome_leakage.css import css_strings


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


def _z_marginal(code, kraus, bit):
    zpos = zchecks.css_parts(code)["zpos"]
    out = {}
    for s, w in enumerate(ex.syndrome_dist(code, kraus, math.pi * bit)):
        key = tuple((s >> j) & 1 for j in zpos)
        out[key] = out.get(key, 0.0) + w
    return out


def _gap(a, b):
    return max(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in set(a) | set(b))


def _small():
    toric = css_strings(*_toric(2), "toric L=2")
    return [STANDARD["hamming_7"](), STANDARD["steane"](), ex.surface_code_3(), toric.with_logical(1)]


def test_the_classical_process_is_the_exact_z_marginal():
    for code in _small():
        for bit in (0, 1):
            exact = _z_marginal(code, sampled.mixture_channel(0.3, 0.0), bit)
            assert _gap(exact, zchecks.exact_distribution(code, 0.3, bit)) < 1e-12, code.name


def test_the_z_checks_do_not_see_t2():
    code = ex.surface_code_3()
    T1, t = 100e-6, 30e-6
    reset = 1 - math.exp(-t / T1)
    for T2 in (60e-6, 150e-6, 200e-6):
        for bit in (0, 1):
            assert _gap(_z_marginal(code, channel_from_t1t2(T1, T2, t), bit),
                        _z_marginal(code, sampled.mixture_channel(reset, 0.0), bit)) < 1e-12


def test_the_sampler_follows_the_enumeration():
    for code in (ex.surface_code_3(), STANDARD["hamming_7"]()):
        s0, s1 = zchecks.sample_syndromes(code, 0.3, 40000, seed=2)
        for bit, s in ((0, s0), (1, s1)):
            exact = zchecks.exact_distribution(code, 0.3, bit)
            keys = sorted(exact)
            index = {k: i for i, k in enumerate(keys)}
            counts = np.zeros(len(keys))
            for row, c in zip(*np.unique(s, axis=0, return_counts=True)):
                counts[index[tuple(int(v) for v in row)]] += c
            assert sampled.fit_p_value(counts, np.array([exact[k] for k in keys])) > 0.001, (code.name, bit)


def test_counts_are_fixed_by_the_seed():
    code = ex.surface_code_3()
    a = zchecks.sample_counts(code, 0.3, 3000, seed=4, chunk=1000)
    b = zchecks.sample_counts(code, 0.3, 3000, seed=4, chunk=1000)
    assert all(np.array_equal(x, y) for x, y in zip(a, b)) and a[0].sum() == 3000


def test_a_code_that_is_not_css_is_refused():
    try:
        zchecks.css_parts(STANDARD["five_qubit"]())
    except ValueError as e:
        assert "not CSS" in str(e)
    else:
        raise AssertionError("a non-CSS code was accepted")


def test_the_integer_program_order_equals_the_enumeration():
    needs("scipy")
    codes = [STANDARD[n]() for n in ("repetition", "code_4_1_2", "hamming_7", "steane")]
    codes += [shor_code(), ex.surface_code_3(), css_strings(*ex.rotated_surface_code(5), "surface d=5"),
              css_strings(*_toric(3), "toric L=3").with_logical(0)]
    for code in codes:
        leaks, order = analytic_leak(code)
        got = zchecks.leak_order(code)
        assert got["order"] == (order if leaks else None), (code.name, got["order"], order)
        assert got["optimal"]


def test_the_integer_program_distance():
    needs("scipy")
    assert zchecks.z_distance(css_strings(*ex.rotated_surface_code(7), "surface d=7"))["distance"] == 7
    assert zchecks.z_distance(css_strings(*_toric(4), "toric L=4"))["distance"] == 4


def test_the_decoded_parity_against_its_exact_value():
    needs("ldpc")
    code = ex.surface_code_3()
    exact = [zchecks.exact_distribution(code, 0.5, b, statistic=lambda s: zchecks.decoded_parity(code, s, 0.5))
             for b in (0, 1)]
    s0, s1 = zchecks.sample_syndromes(code, 0.5, 100000, seed=9)
    for q, s in zip(exact, (s0, s1)):
        p1 = q.get(1, 0.0)
        got = zchecks.decoded_parity(code, s, 0.5).mean()
        assert abs(got - p1) < 4 * math.sqrt(p1 * (1 - p1) / 100000)


def test_matching_and_bposd_agree_on_the_surface_code():
    needs("ldpc")
    needs("pymatching")
    code = ex.surface_code_3()
    s0, _s1 = zchecks.sample_syndromes(code, 0.3, 2000, seed=1)
    assert np.array_equal(zchecks.decoded_parity(code, s0, 0.3, "bposd"),
                          zchecks.decoded_parity(code, s0, 0.3, "matching"))


def test_the_fast_sampler_agrees_with_the_stim_sampler():
    needs("stim")
    from syndrome_leakage.estimate import permutation_test
    code = ex.surface_code_3()
    zmask = np.array(["Z" in s for s in code.stab_strings])
    f0, f1 = zchecks.sample_syndromes(code, 0.3, 5000, seed=3)
    r0, r1 = sampled.sample_records(code, 0.3, 0.0, 5000, seed=4)
    for a, b in ((f0, r0[:, zmask]), (f1, r1[:, zmask])):
        m = a.shape[1] + 1
        pt = permutation_test(np.bincount(a.sum(1), minlength=m), np.bincount(b.sum(1), minlength=m),
                              boots=300, seed=0)
        assert pt["p_value"] > 0.001


def test_the_command_line_takes_the_order_from_the_integer_program():
    needs("scipy")
    import q2sl
    code = css_strings(*ex.rotated_surface_code(7), "surface d=7")
    assert q2sl._order_by_program(code, "refused") == "leaks at order gamma^7 (integer program)"
    assert q2sl._order_by_program(STANDARD["steane"](), "refused") == "protected (integer program)"
