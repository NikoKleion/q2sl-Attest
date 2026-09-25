# logical error after recovery: projector form, single-qubit correction, distance, scaling
import math

import numpy as np

import syndrome_leakage as sl
from syndrome_leakage import expectations as ex
from syndrome_leakage import protection as pr
from syndrome_leakage.channels import amplitude_damping, coherent_diagonal, depolarizing
from syndrome_leakage.core import op
from syndrome_leakage.css import code_distance, css_matrices
from syndrome_leakage.hardware import shor_code

STD = sl.codes.STANDARD
STATES = [(0.0, 0.0), (math.pi, 0.0), (math.pi / 2, math.pi / 3)]


def _distance_3():
    return [STD["five_qubit"](), STD["steane"](), shor_code(), ex.surface_code_3()]


def test_recovery_fidelity_form_matches_projectors():
    for name in STD:
        c = STD[name]()
        table = pr.recovery_strings(c)
        for K in (amplitude_damping(0.2), depolarizing(0.1), coherent_diagonal(0.3)):
            for t, phi in STATES:
                psi = ex.logical_state_vector(c, t, phi)
                rho = c.apply(np.outer(psi, psi.conj()), K)
                out = np.zeros_like(rho)
                for bits, P in c.PROJ.items():
                    R = op(table.get(bits, "I" * c.n))
                    out += R @ P @ rho @ P @ R
                ref = 1.0 - float(np.real(np.vdot(psi, out @ psi)))
                assert abs(ref - pr.logical_error(c, K, t, phi)) < 1e-12, (name, t, phi)


def test_every_single_qubit_pauli_corrected_at_distance_3():
    for c in _distance_3():
        for t, phi in STATES:
            psi = ex.logical_state_vector(c, t, phi)
            for q in range(c.n):
                for p in "XYZ":
                    v = ex.apply_pauli_string(psi, "I" * q + p + "I" * (c.n - q - 1))
                    assert abs(1.0 - pr.recovered_fidelity(c, np.outer(v, v.conj()), psi)) < 1e-12, (c.name, q, p)


def test_distance():
    want = {"repetition": 1, "code_4_1_2": 2, "hamming_7": 2, "five_qubit": 3, "steane": 3}
    for name, d in want.items():
        assert pr.distance(STD[name]()) == d, name
    for c in (STD["code_4_1_2"](), STD["hamming_7"](), STD["steane"](), shor_code(), ex.surface_code_3()):
        assert pr.distance(c) == code_distance(*css_matrices(c.stab_strings)), c.name


def test_hamming_7_weight_two_logical():
    c = STD["hamming_7"]()
    assert not any(pr.syndrome_of(c, "XIIIXII"))
    assert abs(abs(np.vdot(c.V1, op("XIIIXII") @ c.V0)) - 1.0) < 1e-12


def test_logical_error_slope():
    for c in (STD["steane"](), shor_code()):
        for mk in (amplitude_damping, depolarizing):
            a = max(pr.protection_row(c, mk(1e-3)).values())
            b = max(pr.protection_row(c, mk(2e-3)).values())
            assert abs(math.log(b / a) / math.log(2) - 2.0) < 0.05, (c.name, mk.__name__)
    c = STD["code_4_1_2"]()
    a = max(pr.protection_row(c, depolarizing(1e-3)).values())
    b = max(pr.protection_row(c, depolarizing(2e-3)).values())
    assert abs(math.log(b / a) / math.log(2) - 1.0) < 0.05


def test_no_noise_no_error():
    for c in _distance_3():
        assert max(pr.protection_row(c, [np.eye(2)]).values()) < 1e-12


def test_distance_past_weight_four():
    from syndrome_leakage.css import css_strings
    c5 = css_strings(*ex.rotated_surface_code(5), name="surface d=5")
    assert pr.distance(c5) == 5
    c7 = css_strings(*ex.rotated_surface_code(7), name="surface d=7")
    assert pr.distance(c7, max_coset=2 ** 10) is None       # both searches refused, no hang


def test_default_recovery_is_the_table_decoder():
    from syndrome_leakage.channels import depolarizing as dep
    for c in (STD["steane"](), shor_code()):
        a = pr.protection_row(c, dep(0.01))
        b = pr.protection_row(c, dep(0.01), decoder=pr.table_decoder(c))
        for k in a:
            assert abs(a[k] - b[k]) < 1e-15, (c.name, k)


def test_outcome_scores_a_correction():
    c = STD["steane"]()
    e = "XIIIIII"
    assert pr.outcome(c, e, e) == "corrected"
    assert pr.outcome(c, e, "I" * 7) == "detected"
    assert pr.outcome(c, e, "IXXXXXX") in ("logical", "corrected")
    assert pr.outcome(c, c.xl_str, "I" * 7) == "logical"


def test_bposd_decoder_beats_the_table_and_reaches_weight_two():
    from _harness import needs
    needs("ldpc")
    from syndrome_leakage.channels import depolarizing as dep
    from syndrome_leakage.css import css_strings
    for c in (STD["steane"](), shor_code()):
        dec = pr.css_decoder(c)
        for q in range(c.n):
            for p in "XYZ":
                e = "I" * q + p + "I" * (c.n - q - 1)
                assert pr.outcome(c, e, dec(pr.syndrome_of(c, e))) == "corrected", (c.name, q, p)
        table = max(pr.protection_row(c, dep(0.01)).values())
        osd = max(pr.protection_row(c, dep(0.01), decoder=dec).values())
        assert osd < table, (c.name, osd, table)
    big = css_strings(*ex.rotated_surface_code(5), name="d=5")
    dec, tab = pr.css_decoder(big), pr.table_decoder(big)
    rng = np.random.default_rng(0)
    wins = 0
    for _ in range(30):
        e = ["I"] * big.n
        for q, p in zip(rng.choice(big.n, 2, replace=False), rng.choice(list("XYZ"), 2)):
            e[q] = p
        e = "".join(e)
        bits = pr.syndrome_of(big, e)
        assert pr.outcome(big, e, dec(bits)) == "corrected", e
        wins += pr.outcome(big, e, tab(bits)) != "corrected"
    assert wins > 20, wins
