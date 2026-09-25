# repeated extraction as one channel on the logical qubit
import numpy as np

import syndrome_leakage as sl
from syndrome_leakage import expectations as ex
from syndrome_leakage import protection as pr
from syndrome_leakage import rounds as rd
from syndrome_leakage.channels import amplitude_damping, channel_from_t1t2, dephasing, per_qubit_damping
from syndrome_leakage.core import op
from syndrome_leakage.eavesdrop import repeated_extraction

STD = sl.codes.STANDARD


def test_pauli_action_matches_the_matrix():
    rng = np.random.default_rng(0)
    for s in ("XZI", "YYZ", "ZZII", "XYZI"):
        d = 2 ** len(s)
        rho = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        act, G = rd.pauli_action(s), op(s)
        assert np.abs(rd.left(rho, act) - G @ rho).max() < 1e-12
        assert np.abs(rd.right(rho, act) - rho @ G).max() < 1e-12


def test_coset_decoder_explains_every_syndrome():
    for name in ("repetition", "code_4_1_2", "five_qubit", "steane"):
        c = STD[name]()
        assert rd.is_complete(c, rd.coset_decoder(c)), name
    assert rd.is_complete(ex.surface_code_3(), rd.coset_decoder(ex.surface_code_3()))
    assert not rd.is_complete(STD["steane"](), pr.table_decoder(STD["steane"]()))


def test_against_the_density_matrix_path():
    K = amplitude_damping(0.2)
    for name in ("repetition", "five_qubit"):
        c = STD[name]()
        dec = pr.table_decoder(c)
        assert rd.is_complete(c, dec), name
        old = repeated_extraction(c, kraus_override=K, rounds=5, correct=True)["rounds"]
        new = rd.iterate(c, K, rounds=5, decoder=dec)
        for a, b in zip(old, new):
            assert abs(a["tvd"] - b["tvd"]) < 1e-12, (name, a["round"])
            assert abs(a["chernoff"] - b["chernoff"]) < 1e-9, (name, a["round"])


def test_pauli_channel_never_leaks_over_rounds():
    for name in ("repetition", "code_4_1_2"):
        c = STD[name]()
        rows = rd.iterate(c, dephasing(0.05), rounds=10)
        assert max(r["tvd"] for r in rows) < 1e-12, name


def test_higher_order_sums_settle_sooner():
    K = channel_from_t1t2(163.3e-6, 77.0e-6, 1e-6)
    tail = {}
    for name in ("repetition", "code_4_1_2", "hamming_7"):
        c = STD[name]()
        rows = rd.iterate(c, K, rounds=400, decoder=rd.coset_decoder(c))
        early, late = rows[99]["cumulative"], rows[399]["cumulative"]
        assert late >= early
        tail[name] = (late - early) / late
    assert tail["repetition"] > tail["code_4_1_2"] > tail["hamming_7"]
    assert tail["repetition"] > 0.1 and tail["hamming_7"] < 1e-3


def test_one_channel_per_qubit_runs():
    c = STD["repetition"]()
    a = rd.iterate(c, amplitude_damping(0.1), rounds=3)
    b = rd.iterate(c, per_qubit_damping([0.1] * 3), rounds=3)
    for x, y in zip(a, b):
        assert abs(x["tvd"] - y["tvd"]) < 1e-12
    mixed = rd.iterate(c, per_qubit_damping([0.05, 0.1, 0.2]), rounds=3)
    assert abs(mixed[0]["tvd"] - a[0]["tvd"]) > 1e-6


def test_rounds_on_a_code_with_two_logical_qubits():
    from syndrome_leakage.css import css_from_matrices
    c = css_from_matrices(np.array([[1, 1, 1, 1]], np.uint8), np.array([[1, 1, 1, 1]], np.uint8), "[[4,2,2]]")
    assert c.k == 2
    for i in range(c.k):
        ci = c.with_logical(i)
        rows = rd.iterate(ci, amplitude_damping(0.1), rounds=4, decoder=rd.coset_decoder(ci))
        assert all(r["tvd"] > 0 for r in rows)
        assert rows[0]["tvd"] > rows[-1]["tvd"]
    mixed = rd.iterate(c.with_logical(0), per_qubit_damping([0.05, 0.1, 0.15, 0.2]), rounds=3,
                       decoder=rd.coset_decoder(c))
    flat = rd.iterate(c.with_logical(0), amplitude_damping(0.1), rounds=3, decoder=rd.coset_decoder(c))
    assert abs(mixed[0]["tvd"] - flat[0]["tvd"]) > 1e-6
