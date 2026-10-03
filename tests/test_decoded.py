# the leak through a decoder's output
import itertools
import json
import math
import os

import numpy as np

import syndrome_leakage as sl
from syndrome_leakage import decoded as dc
from syndrome_leakage import estimate as es
from syndrome_leakage import expectations as ex
from syndrome_leakage import protection as pr
from syndrome_leakage import rounds as rd
from syndrome_leakage.channels import amplitude_damping, depolarizing
from syndrome_leakage.codes import shor_code
from syndrome_leakage.core import tvd

STD = sl.codes.STANDARD
RESULTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")


def _codes():
    return [STD[n]() for n in ("repetition", "code_4_1_2", "hamming_7", "five_qubit", "steane")] + \
        [shor_code(), ex.surface_code_3()]


def _bits(s, m):
    return tuple((s >> j) & 1 for j in range(m))


def test_least_weight_decoder_against_every_pauli():
    for name in ("repetition", "code_4_1_2", "five_qubit", "hamming_7"):
        c = STD[name]()
        best = {}
        for e in itertools.product("IXYZ", repeat=c.n):
            s, w = pr.syndrome_of(c, "".join(e)), sum(ch != "I" for ch in e)
            best[s] = min(best.get(s, c.n + 1), w)
        dec = dc.min_weight_decoder(c)
        assert len(best) == 2 ** len(c.stab_strings), name
        for s, w in best.items():
            r = dec(s)
            assert pr.syndrome_of(c, r) == s and sum(ch != "I" for ch in r) == w, (name, s, r)


def test_a_complete_decoder_keeps_the_whole_leak():
    # a correction that carries its syndrome differs for every syndrome, so its distribution is the syndrome's
    for c in _codes():
        m = len(c.stab_strings)
        leak, d0, d1 = ex.population_leak(c, amplitude_damping(0.2))
        for make in (rd.coset_decoder, dc.min_weight_decoder):
            dec = make(c)
            assert rd.is_complete(c, dec), c.name
            lab = dc.view_labels(c, dec)
            assert len(set(lab["correction"])) == 2 ** m, c.name
            assert abs(dc.view_leaks(d0, d1, lab)["correction"] - leak) < 1e-15, c.name


def test_no_view_exceeds_the_syndrome_leak():
    for c in _codes():
        for g in (0.2, 0.05):
            leak, d0, d1 = ex.population_leak(c, amplitude_damping(g))
            for make in (pr.table_decoder, rd.coset_decoder, dc.min_weight_decoder):
                v = dc.view_leaks(d0, d1, dc.view_labels(c, make(c)))
                assert v["syndrome"] == leak
                assert max(v.values()) <= leak + 1e-15, (c.name, g, v)


def test_a_pauli_channel_leaks_through_no_view():
    for c in _codes():
        _, d0, d1 = ex.population_leak(c, depolarizing(0.1))
        for make in (pr.table_decoder, dc.min_weight_decoder):
            assert max(dc.view_leaks(d0, d1, dc.view_labels(c, make(c))).values()) < 1e-14, c.name


def test_the_single_qubit_table_merges_syndromes():
    # 22 corrections for 256 syndromes on the Shor code, and three eighths of the leak
    c = shor_code()
    lab = dc.view_labels(c, pr.table_decoder(c))
    assert len(set(lab["correction"])) == 22
    for g in (0.2, 0.05):
        leak, d0, d1 = ex.population_leak(c, amplitude_damping(g))
        assert abs(dc.view_leaks(d0, d1, lab)["correction"] / leak - 0.375) < 1e-9


def test_frame_spectrum_against_direct_sums():
    c = ex.surface_code_3()
    m = len(c.stab_strings)
    _, d0, d1 = ex.population_leak(c, amplitude_damping(0.05))
    dec = dc.min_weight_decoder(c)
    spec = dc.frame_spectrum(d0, d1, dc.view_labels(c, dec)["frame_z"])
    for a in (0, 1, 37, 160, 255):
        rep = dc.representative(c, c.zl_str, a)
        bit = np.array([int(pr._anticommutes(dec(_bits(s, m)), rep)) for s in range(2 ** m)])
        assert abs(spec[a] - tvd(dc.push(d0, bit), dc.push(d1, bit))) < 1e-15, a
    assert len(dc.type_masks(c, "Z")) == 16 and len(dc.type_masks(c, "X")) == 16


def test_shor_leak_is_the_closed_form_and_sits_in_one_bit():
    # the parity of the number of X corrections is the frame bit of Z on every qubit
    c = shor_code()
    frame = dc.view_labels(c, dc.min_weight_decoder(c))["frame_z"]
    every = next(a for a in dc.type_masks(c, "Z") if dc.representative(c, c.zl_str, a) == "Z" * c.n)
    for per in ([0.2] * 9, [0.05] * 9, [0.01] * 9, [0.02 + 0.01 * q for q in range(9)]):
        leak, d0, d1 = ex.population_leak(c, [amplitude_damping(x) for x in per])
        assert abs(leak - dc.shor_closed_form(per)) < 1e-14
        spec = dc.frame_spectrum(d0, d1, frame)
        assert abs(spec[every] - leak) < 1e-14
        if len(set(per)) == 1:
            assert abs(spec[0] / leak - 1 / 27) < 1e-9


def test_every_view_keeps_the_leak_order():
    c = shor_code()
    lab = dc.view_labels(c, dc.min_weight_decoder(c))
    a = dc.view_leaks(*ex.population_leak(c, amplitude_damping(0.005))[1:], lab)
    b = dc.view_leaks(*ex.population_leak(c, amplitude_damping(0.01))[1:], lab)
    for view in ("syndrome",) + dc.VIEWS:
        assert abs(math.log2(b[view] / a[view]) - 3.0) < 0.05, view


def test_views_of_the_z_record():
    c = shor_code()
    zbits = [j for j, s in enumerate(c.stab_strings) if set(s) <= {"I", "Z"}]
    dec = dc.min_weight_decoder(c)
    lab = dc.view_labels(c, dec, checks=zbits)
    assert len(lab["correction"]) == 64 and len(set(lab["correction"])) == 64
    assert sorted(set(lab["weight"])) == [0, 1, 2, 3]
    for r in range(64):
        bits = [0] * len(c.stab_strings)
        for j, g in enumerate(zbits):
            bits[g] = (r >> j) & 1
        assert set(dec(tuple(bits))) <= {"I", "X"}
    d = np.random.default_rng(0).dirichlet(np.ones(64))
    assert np.allclose(dc.transform(lab["weight"])(d), dc.push(d, lab["weight"])) and abs(dc.push(d, lab["acted"]).sum() - 1) < 1e-12


def test_pinned_run_keeps_its_distance_in_the_parity_bit():
    # ibm_fez job daquif6ekp0c73arbd70: 64 outcomes against one bit, at the three delays above zero
    R = json.load(open(os.path.join(RESULTS, "hardware_shor_pinned_ibm_fez.json")))
    c = shor_code()
    zbits = [j for j, s in enumerate(c.stab_strings) if set(s) <= {"I", "Z"}]
    parity = dc.view_labels(c, dc.min_weight_decoder(c), checks=zbits)["weight"] % 2
    for t, record, bit in ((2e-05, 0.0438, 0.0367), (5e-05, 0.1673, 0.1671), (0.0001, 0.2337, 0.2337)):
        d0 = es.dist_from_counts(R["counts"][f"{t}_0"], 64)[0]
        d1 = es.dist_from_counts(R["counts"][f"{t}_1"], 64)[0]
        assert abs(tvd(d0, d1) - record) < 5e-5
        assert abs(tvd(dc.push(d0, parity), dc.push(d1, parity)) - bit) < 5e-5


def test_exact_z_distribution_against_the_enumeration():
    from syndrome_leakage import zchecks as zc
    for c in (shor_code(), ex.surface_code_3()):
        mz = len(zc.css_parts(c)["Hz"])
        for p in (0.3, 0.5, 0.52, [0.05 + 0.04 * q for q in range(c.n)]):
            for bit in (0, 1):
                ref = np.zeros(2 ** mz)
                for key, w in zc.exact_distribution(c, p, bit).items():
                    ref[sum(int(b) << j for j, b in enumerate(key))] += w
                assert np.abs(zc.exact_z_distribution(c, p, bit) - ref).max() < 1e-13, (c.name, p, bit)


def test_least_weight_x_against_every_x_error():
    from syndrome_leakage import zchecks as zc
    c = ex.surface_code_3()
    Hz = zc.css_parts(c)["Hz"]
    weight, error = zc.least_weight_x(c)
    best = {}
    for e in itertools.product((0, 1), repeat=c.n):
        s = sum(int(b) << j for j, b in enumerate(Hz @ np.array(e) % 2))
        best[s] = min(best.get(s, c.n), sum(e))
    assert all(weight[s] == w for s, w in best.items()) and len(best) == len(weight)
    for s, mask in enumerate(error):
        e = np.array([(int(mask) >> q) & 1 for q in range(c.n)])
        assert sum(int(b) << j for j, b in enumerate(Hz @ e % 2)) == s and e.sum() == weight[s]


def test_distance_five_z_record():
    # order 5, and the one bit that holds the whole leak at distance 3 holds less than the fired count here
    from syndrome_leakage import zchecks as zc
    from syndrome_leakage.css import css_strings
    c = css_strings(*ex.rotated_surface_code(5), name="rotated surface d=5")
    leak = lambda p: tvd(zc.exact_z_distribution(c, p, 0), zc.exact_z_distribution(c, p, 1))
    assert 4.7 < math.log2(leak(0.02) / leak(0.01)) < 5.0
    weight, _ = zc.least_weight_x(c)
    fired = np.array([bin(s).count("1") for s in range(len(weight))])
    d0, d1 = zc.exact_z_distribution(c, 0.3075, 0), zc.exact_z_distribution(c, 0.3075, 1)
    by_fired = tvd(dc.push(d0, fired), dc.push(d1, fired))
    by_parity = tvd(dc.push(d0, weight % 2), dc.push(d1, weight % 2))
    assert abs(by_fired - 0.00236) < 1e-5 and by_parity < by_fired < 0.2 * tvd(d0, d1)
    s3 = ex.surface_code_3()
    w3, _ = zc.least_weight_x(s3)
    a, b = zc.exact_z_distribution(s3, 0.3075, 0), zc.exact_z_distribution(s3, 0.3075, 1)
    assert abs(tvd(dc.push(a, w3 % 2), dc.push(b, w3 % 2)) - tvd(a, b)) < 1e-14
