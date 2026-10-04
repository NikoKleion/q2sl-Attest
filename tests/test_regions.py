# the leak of a part of the record: chosen generators, the smallest sets that leak, one qubit made ideal
import itertools

import numpy as np

import syndrome_leakage as sl
from _harness import needs
from syndrome_leakage import expectations as ex
from syndrome_leakage import regions as rg
from syndrome_leakage import zchecks
from syndrome_leakage.channels import amplitude_damping, depolarizing
from syndrome_leakage.codes import shor_code
from syndrome_leakage.core import tvd
from syndrome_leakage.css import css_strings
from syndrome_leakage.hardware import marginal_dist

STD = sl.codes.STANDARD


def _subsets(items):
    for size in range(1, len(items) + 1):
        yield from itertools.combinations(items, size)


def _z(code):
    return zchecks.css_parts(code)["zpos"]


def _surface(d):
    return css_strings(*ex.rotated_surface_code(d), name=f"rotated surface d={d}")


def test_region_equals_the_marginal_of_the_whole_record():
    # the products of the chosen generators alone against the whole record summed over the other bits
    cases = [(shor_code(), amplitude_damping(0.1)), (ex.surface_code_3(), amplitude_damping(0.2)),
             (STD["five_qubit"](), amplitude_damping(0.15)), (STD["code_4_1_2"](), depolarizing(0.1)),
             (STD["code_4_1_2"](), [amplitude_damping(g) for g in (0.05, 0.1, 0.2, 0.3)])]
    for code, kraus in cases:
        m = len(code.stab_strings)
        _, d0, d1 = ex.population_leak(code, kraus)
        for c in _subsets(range(m)):
            got = rg.region_leak(code, kraus, list(c))
            assert abs(got[0] - tvd(marginal_dist(d0, c), marginal_dist(d1, c))) < 1e-12, (code.name, c)
            assert np.allclose(got[1], marginal_dist(d0, c), atol=1e-12)
    whole = rg.region_leak(shor_code(), amplitude_damping(0.1), "all")[0]
    assert abs(whole - ex.population_leak(shor_code(), amplitude_damping(0.1))[0]) < 1e-12


def test_region_with_readout_error_equals_the_marginal():
    code, kraus, q = ex.surface_code_3(), amplitude_damping(0.2), np.linspace(0.01, 0.08, 8)
    _, d0, d1 = ex.population_leak(code, kraus, readout=q)
    for c in ((4, 6), (0, 4, 5, 7), (4, 5, 6, 7)):
        assert abs(rg.region_leak(code, kraus, list(c), readout=q)[0]
                   - tvd(marginal_dist(d0, c), marginal_dist(d1, c))) < 1e-12


def test_a_larger_region_never_leaks_less():
    for code in (shor_code(), ex.surface_code_3()):
        m = len(code.stab_strings)
        _, d0, d1 = ex.population_leak(code, amplitude_damping(0.1))
        leak = {c: tvd(marginal_dist(d0, c), marginal_dist(d1, c)) for c in _subsets(range(m))}
        for c, v in leak.items():
            for extra in set(range(m)) - set(c):
                assert leak[tuple(sorted(c + (extra,)))] >= v - 1e-13, (code.name, c, extra)


def test_x_generators_alone_carry_nothing_under_damping():
    for code in (shor_code(), ex.surface_code_3()):
        assert rg.region_leak(code, amplitude_damping(0.3), "x")[0] < 1e-13
        assert rg.region_leak(code, amplitude_damping(0.3), "z")[0] > 1e-3


def test_smallest_leaking_sets_of_shor_and_the_surface_code():
    for g in (0.1, 0.3):
        t = (g * (1 - g)) ** 3
        r = rg.leaking_sets(shor_code(), amplitude_damping(g))
        assert r["size"] == 3 and len(r["sets"]) == 8 and abs(r["whole"] - (1 - (1 - g) ** 3 - g ** 3) ** 3) < 1e-12
        blocks = [{2, 3}, {4, 5}, {6, 7}]
        for c, v in r["sets"]:
            assert all(len(set(c) & b) == 1 for b in blocks) and abs(v - 8 * t) < 1e-12, (c, v)
        r = rg.leaking_sets(ex.surface_code_3(), amplitude_damping(g))
        assert r["size"] == 2 and [c for c, _v in r["sets"]][2] == (4, 5)
        assert {c for c, _v in r["sets"][:2]} == {(4, 6), (5, 7)}
        assert np.allclose([v for _c, v in r["sets"]], [4 * t, 4 * t, 2 * t], atol=1e-12)
    assert rg.leaking_sets(STD["steane"](), amplitude_damping(0.2))["size"] is None
    assert rg.leaking_sets(shor_code(), amplitude_damping(0.1), max_size=2)["size"] is None


def test_z_region_equals_the_exact_engine_at_unequal_rates():
    for code in (shor_code(), ex.surface_code_3(), STD["code_4_1_2"]()):
        p = np.linspace(0.03, 0.4, code.n)
        chans = [amplitude_damping(x) for x in p]
        for c in _subsets(_z(code)):
            a, b = rg.z_region_leak(code, p, list(c)), rg.region_leak(code, chans, list(c))
            assert abs(a[0] - b[0]) < 1e-12 and np.allclose(a[1], b[1], atol=1e-12), (code.name, c)


def test_z_region_of_the_distance_5_code_equals_its_whole_z_record():
    code = _surface(5)
    zs = _z(code)
    full = [zchecks.exact_z_distribution(code, 0.1, b) for b in (0, 1)]
    assert abs(rg.z_region_leak(code, 0.1, "z")[0] - tvd(*full)) < 1e-13
    for c in ((18, 19, 23), (12, 13, 16, 17, 20, 21), (14, 15, 22), (12, 16)):
        where = [zs.index(g) for g in c]
        ref = tvd(marginal_dist(full[0], where), marginal_dist(full[1], where))
        assert abs(rg.z_region_leak(code, 0.1, list(c))[0] - ref) < 1e-13, c


def test_z_region_reaches_a_code_of_121_qubits():
    # a column of six Z checks of the distance 11 code, 22 qubits touched: a lower bound on the whole record
    code = _surface(11)
    leak, d0, d1 = rg.z_region_leak(code, 0.1, [60, 61, 62, 63, 64, 110])
    assert abs(leak - 2.008e-10) < 5e-13 and abs(d0.sum() - 1) < 1e-12 and len(d0) == 64
    assert rg.z_region_leak(code, 0.1, [60, 61, 62, 63, 64])[0] < 1e-15
    try:
        rg.z_region_dist(code, 0.1, 0, "z")
    except ValueError as e:
        assert "too many" in str(e)
    else:
        raise AssertionError("the whole Z record of the distance 11 code was accepted")


def test_z_region_rejects_an_x_generator():
    try:
        rg.z_region_leak(shor_code(), 0.1, [0, 2])
    except ValueError as e:
        assert "not of Z type" in str(e)
    else:
        raise AssertionError("an X generator was accepted")


def test_a_region_leaks_exactly_when_a_product_covers_a_logical():
    needs("scipy")
    for code in (shor_code(), ex.surface_code_3()):
        seen = set()
        for c in _subsets(_z(code)):
            leaks = rg.z_region_leak(code, 0.1, list(c))[0] > 1e-12
            order = rg.region_order(code, list(c))["order"]
            assert leaks == (order is not None), (code.name, c, order)
            seen.add((leaks, order))
        assert (True, 3) in seen and (False, None) in seen


def test_smallest_leaking_set_by_integer_program():
    needs("scipy")
    for code, size in ((shor_code(), 3), (ex.surface_code_3(), 2), (_surface(5), 3), (STD["steane"](), None),
                       (STD["repetition"](), 1)):
        r = rg.smallest_leaking_set(code)
        assert r["size"] == size and r["optimal"], (code.name, r)
        if size is not None:
            assert rg.z_region_leak(code, 0.1, r["checks"])[0] > 1e-12
            assert all(rg.z_region_leak(code, 0.1, list(c))[0] < 1e-12
                       for c in itertools.combinations(_z(code), size - 1)) or size == 1
    assert rg.region_order(_surface(5), [18, 19, 23])["order"] == 5


def test_qubit_drops():
    g = amplitude_damping(0.1)
    r = rg.qubit_drops(ex.surface_code_3(), g)
    assert abs(r["leak"] - ex.population_leak(ex.surface_code_3(), g)[0]) < 1e-12
    assert int(np.argmax(r["drop"])) == 4 and min(r["drop"]) > 0 and abs(sum(r["drop"]) - r["leak"]) > 1e-3
    r = rg.qubit_drops(shor_code(), g)
    assert np.allclose(r["drop"], r["drop"][0], atol=1e-12) and r["drop"][0] > 0
    # one check per block of Shor touches qubits 0, 1, 4, 5, 6, 7: the other three change nothing
    r = rg.qubit_drops(shor_code(), g, checks=[2, 5, 6])
    assert all(abs(r["drop"][q]) < 1e-13 for q in (2, 3, 8)) and all(r["drop"][q] > 1e-4 for q in (0, 1, 4, 5, 6, 7))
    z = rg.z_qubit_drops(shor_code(), 0.1, [2, 5, 6])
    assert np.allclose(z["drop"], r["drop"], atol=1e-12) and abs(z["leak"] - r["leak"]) < 1e-12


def test_a_drop_changes_sign_at_strong_damping():
    # Shor: the leak is F^3 and F^2 F' with one qubit noiseless, F = 1 - (1-g)^3 - g^3, F' = 1 - (1-g)^2; F = F' at 0.5
    for g, sign in ((0.3, 1), (0.5, 0), (0.6, -1)):
        F, Fq = 1 - (1 - g) ** 3 - g ** 3, 1 - (1 - g) ** 2
        r = rg.qubit_drops(shor_code(), amplitude_damping(g))
        assert np.allclose(r["drop"], F * F * (F - Fq), atol=1e-12) and np.sign(round(r["drop"][0], 12)) == sign
    r = rg.qubit_drops(ex.surface_code_3(), amplitude_damping(0.45))
    assert r["drop"][4] > 0.02 and sum(x < -1e-4 for x in r["drop"]) == 8


def test_css_region_equals_the_exact_engine_on_sets_of_x_and_z_generators():
    for code in (shor_code(), ex.surface_code_3(), STD["code_4_1_2"](), STD["steane"]()):
        p = np.linspace(0.03, 0.4, code.n)
        chans = [amplitude_damping(x) for x in p]
        for c in _subsets(range(len(code.stab_strings))):
            a, b = rg.css_region_leak(code, p, list(c)), rg.region_leak(code, chans, list(c))
            assert abs(a[0] - b[0]) < 1e-12 and np.allclose(a[1], b[1], atol=1e-12), (code.name, c)
    v = np.random.default_rng(0).normal(size=64)
    assert np.allclose(rg._walsh(v), np.real(ex.walsh_hadamard(v)), atol=1e-12)


def test_css_region_under_t1_and_t2():
    # X and Y shrink by exp(-t / T2), so a set that holds X generators depends on T2 and a set of Z generators does not
    from syndrome_leakage.channels import channel_from_t1t2
    code, t = ex.surface_code_3(), 40e-6
    T1, T2 = np.linspace(80e-6, 250e-6, 9), np.linspace(60e-6, 120e-6, 9)
    chans = [channel_from_t1t2(a, b, t) for a, b in zip(T1, T2)]
    p, coh = 1 - np.exp(-t / T1), np.exp(-t / T2)
    for c in ((0, 4, 6), "all", (2, 5, 7), "z"):
        want = rg.region_leak(code, chans, c if isinstance(c, str) else list(c))[0]
        assert abs(rg.css_region_leak(code, p, c if isinstance(c, str) else list(c), coherence=coh)[0] - want) < 1e-12
    assert abs(rg.css_region_leak(code, p, "z", coherence=coh)[0] - rg.css_region_leak(code, p, "z")[0]) < 1e-15
    assert rg.css_region_leak(code, p, "all", coherence=coh)[0] < rg.css_region_leak(code, p, "all")[0] - 1e-3
    assert rg.css_region_leak(code, p, "x", coherence=coh)[0] < 1e-15


def test_css_region_past_the_exact_engine():
    code = _surface(5)
    zs = _z(code)
    xs = [j for j in range(len(code.stab_strings)) if j not in zs]
    z_only = rg.css_region_leak(code, 0.05, zs)[0]
    assert abs(z_only - rg.z_region_leak(code, 0.05, "z")[0]) < 1e-15 and abs(z_only - 1.065691e-05) < 5e-12
    with_x = rg.css_region_leak(code, 0.05, zs + xs[:4], max_terms=2 ** 30)[0]
    assert abs(with_x - 1.147151e-05) < 5e-12
    big, col = _surface(11), [60, 61, 62, 63, 64, 110]
    p = np.random.default_rng(11).uniform(0.02, 0.4, big.n)
    assert abs(rg.css_region_leak(big, p, col)[0] - rg.z_region_leak(big, p, col)[0]) < 1e-15
    for bad in (lambda: rg.css_region_dist(code, 1.0, 0, zs), lambda: rg.css_region_dist(code, 0.05, 0, "all")):
        try:
            bad()
        except ValueError:
            continue
        raise AssertionError("an input outside the engine's reach was accepted")


def test_css_region_leak_keeps_its_digits_at_weak_damping():
    # the difference of the two characteristic functions term by term, against the two records subtracted
    code = ex.surface_code_3()
    p = np.linspace(0.03, 0.4, 9)
    for c in ([4, 6], [0, 4, 6], "all", "z"):
        leak, d0, d1 = rg.css_region_leak(code, p, c)
        assert abs(leak - tvd(d0, d1)) < 1e-13 and abs(leak - rg.css_region_tv(code, p, c)) < 1e-18
    # leak / gamma^d tends to 7 on the distance 3 code, 27 on Shor and 4 on [[4,1,2]], for the Z record and the whole
    for code, d, limit in ((code, 3, 7.0), (shor_code(), 3, 27.0), (STD["code_4_1_2"](), 2, 4.0)):
        for rec in ("z", "all"):
            assert abs(rg.css_region_tv(code, 1e-6, rec) / 1e-6 ** d - limit) < 1e-3, (code.name, rec)
    # and to 52 on the Z record of the distance 5 code, the number of its Z logicals of weight 5
    ratios = [rg.css_region_tv(_surface(5), g, "z") / g ** 5 for g in (1e-3, 1e-4, 1e-5)]
    assert abs(ratios[0] - 51.5662) < 1e-3 and abs(ratios[2] - 52.0) < 0.01 and ratios[0] < ratios[1] < ratios[2]


def test_a_qubit_off_every_least_weight_logical_has_a_negative_drop_at_weak_damping():
    for name, order in (("repetition", 1), ("code_4_1_2", 2), ("hamming_7", 3)):
        code = STD[name]()
        got, cover = rg.logical_cover(code)
        drops = rg.qubit_drops(code, amplitude_damping(1e-3))["drop"]
        assert got == order and all((d > 0) == bool(c) for d, c in zip(drops, cover)), (name, drops, cover)
    assert rg.logical_cover(STD["hamming_7"]())[1].tolist() == [True, True, True, False, True, True, True]
    assert rg.logical_cover(STD["steane"]()) [0] is None


def test_shor_drop_in_closed_form_at_unequal_rates():
    # drop of qubit q = gamma_q (1 - gamma_a - gamma_b) times the other two blocks' factors, a and b its block mates
    rng = np.random.default_rng(12)
    F = lambda g: 1 - np.prod(1 - g) - np.prod(g)
    for _ in range(6):
        g = rng.uniform(0.02, 0.95, 9)
        drops = rg.qubit_drops(shor_code(), [amplitude_damping(x) for x in g])["drop"]
        for q in range(9):
            blk = 3 * (q // 3)
            mates = [i for i in range(blk, blk + 3) if i != q]
            others = np.prod([F(g[b:b + 3]) for b in (0, 3, 6) if b != blk])
            assert abs(drops[q] - others * g[q] * (1 - g[mates].sum())) < 1e-12, (q, g)


def test_noise_on_fewer_qubits_than_a_logical_gives_no_leak():
    code = shor_code()
    for noisy in ((0,), (0, 3), (0, 1, 2), (0, 3, 6), (1, 4, 8)):
        chans = [amplitude_damping(0.2 if q in noisy else 0.0) for q in range(code.n)]
        leak = ex.population_leak(code, chans)[0]
        assert (leak > 1e-6) == (len({q // 3 for q in noisy}) == 3), (noisy, leak)
