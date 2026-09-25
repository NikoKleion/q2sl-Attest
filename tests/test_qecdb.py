# the qecdb.org loader: parsing is tested offline on two trimmed records saved from qecdb.org in
# tests/fixtures; the live fetch runs only when Q2SL_NETWORK=1
import json
import os

import numpy as np

from _harness import skip
from syndrome_leakage import expectations as ex
from syndrome_leakage import load
from syndrome_leakage.analyze import analytic_leak
from syndrome_leakage.channels import amplitude_damping

HERE = os.path.dirname(os.path.abspath(__file__))


def _fixture(name):
    return json.load(open(os.path.join(HERE, "fixtures", name), encoding="utf-8"))


def _commute(a, b):
    return sum(1 for p, q in zip(a, b) if p != "I" and q != "I" and p != q) % 2 == 0


def test_surface_record_reproduces_the_package_surface_code():
    code = load.code_from_qecdb_record(_fixture("qecdb_surface_9_1_3.json"))
    assert (code.n, code.k, len(code.stab_strings)) == (9, 1, 8)
    assert analytic_leak(code)[1] == 3
    leak = ex.population_leak(code, amplitude_damping(0.2))[0]
    ours = ex.population_leak(ex.surface_code_3(), amplitude_damping(0.2))[0]
    assert abs(leak - ours) < 1e-12, (leak, ours)


def test_the_record_logicals_are_kept():
    code = load.code_from_qecdb_record(_fixture("qecdb_surface_9_1_3.json"))
    assert code.logical_pairs == [("IIIIIIZZZ", "IIXIIXIIX")]


def test_logicals_are_sorted_by_type_not_position():
    # the record lists (XXII, ZIZI) then (ZZII, XIXI): the type order flips between pairs
    code = load.code_from_qecdb_record(_fixture("qecdb_4_2_2.json"))
    assert code.k == 2
    for i, (z, x) in enumerate(code.logical_pairs):
        assert set(z) <= {"I", "Z"} and set(x) <= {"I", "X"}
        assert not _commute(z, x)
        for j, (z2, x2) in enumerate(code.logical_pairs):
            if i != j:
                assert _commute(z, x2) and _commute(z2, x)


def test_a_non_css_record_is_refused():
    rec = {"_id": "test", "name": "[[5, 1, 3]]", "k": 1, "H": "XZZXI IXZZX XIXZZ ZXIXZ", "L": "XXXXX ZZZZZ"}
    try:
        load.code_from_qecdb_record(rec)
    except ValueError as e:
        assert "not CSS" in str(e)
    else:
        raise AssertionError("a non-CSS record was accepted")


def test_live_fetch_and_search():
    if os.environ.get("Q2SL_NETWORK") != "1":
        skip("set Q2SL_NETWORK=1 to reach qecdb.org")
    rows = load.qecdb_search(n=9, k=1, d=3)
    assert rows and all(r["n"] == 9 and r["k"] == 1 and r["d"] == 3 for r in rows)
    code = load.from_qecdb(rows[0]["_id"])
    assert code.n == 9 and np.isfinite(analytic_leak(code)[1])
