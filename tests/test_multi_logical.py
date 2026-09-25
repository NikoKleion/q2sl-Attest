# codes with more than one logical qubit
import numpy as np

import syndrome_leakage as sl
from syndrome_leakage import expectations as ex
from syndrome_leakage import protection as pr
from syndrome_leakage.analyze import analytic_leak, population_leak
from syndrome_leakage.channels import amplitude_damping
from syndrome_leakage.core import pmul
from syndrome_leakage.css import _StringCode, code_distance, css_from_matrices, css_strings

FOUR_TWO_TWO = (np.array([[1, 1, 1, 1]], np.uint8), np.array([[1, 1, 1, 1]], np.uint8))


def toric(L):
    n = 2 * L * L
    H = lambda r, c: (r % L) * L + (c % L)
    V = lambda r, c: L * L + (r % L) * L + (c % L)
    Hx = [np.zeros(n, np.uint8) for _ in range(L * L)]
    Hz = [np.zeros(n, np.uint8) for _ in range(L * L)]
    i = 0
    for r in range(L):
        for c in range(L):
            for q in (H(r, c), H(r + 1, c), V(r, c), V(r, c + 1)):
                Hz[i][q] ^= 1
            for q in (H(r, c), H(r, c - 1), V(r, c), V(r - 1, c)):
                Hx[i][q] ^= 1
            i += 1
    return np.array(Hx), np.array(Hz)


def _symplectic(pairs):
    v = lambda s, sym: np.array([1 if ch == sym else 0 for ch in s], np.uint8)
    for i, (zi, _xi) in enumerate(pairs):
        for j, (_zj, xj) in enumerate(pairs):
            want = 1 if i == j else 0
            if int(np.dot(v(zi, "Z"), v(xj, "X"))) % 2 != want:
                return False
    return True


def test_k_two_codes_build_with_paired_logicals():
    for Hx, Hz in (FOUR_TWO_TWO, toric(2)):
        c = css_strings(Hx, Hz, "k=2")
        assert c.k == 2 and len(c.logical_pairs) == 2
        assert _symplectic(c.logical_pairs)
        assert code_distance(Hx, Hz) == 2


def test_projector_path_and_codewords():
    c = css_from_matrices(*FOUR_TWO_TWO, name="[[4,2,2]]")
    assert c.k == 2
    assert abs(np.vdot(c.V0, c.V0) - 1) < 1e-12 and abs(np.vdot(c.V0, c.V1)) < 1e-12
    c.verify_projectors()
    for i in range(c.k):
        ci = c.with_logical(i)
        assert ci.zl_str, ci.xl_str
        assert abs(np.vdot(ci.V0, ci.V1)) < 1e-12


def test_both_engines_agree_for_k_two():
    c = css_from_matrices(*toric(2), name="toric L=2")
    for i in range(c.k):
        ci = c.with_logical(i)
        a = population_leak(ci, amplitude_damping(0.2))[0]
        b = ex.population_leak(ci, amplitude_damping(0.2))[0]
        assert abs(a - b) < 1e-12, (i, a, b)
        assert analytic_leak(ci)[1] == 2
        assert pr.distance(ci) == 2


def test_k_one_codes_are_untouched():
    for name in ("repetition", "five_qubit", "steane", "code_4_1_2", "hamming_7"):
        c = sl.codes.STANDARD[name]()
        assert c.k == 1 and c.logical_pairs == [(c.zl_str, c.xl_str)]


def test_the_leak_belongs_to_the_logical_basis():
    c = css_strings(*toric(2), name="toric L=2")
    (z0, x0), (z1, x1) = c.logical_pairs
    alt = [(pmul(z0, z1), x0), (z1, pmul(x0, x1))]
    other = _StringCode("other basis", c.n, c.stab_strings, alt[0][0], alt[0][1], alt)
    first = [ex.population_leak(c.with_logical(i), amplitude_damping(0.2))[0] for i in range(2)]
    second = [ex.population_leak(other.with_logical(i), amplitude_damping(0.2))[0] for i in range(2)]
    assert abs(first[0] - second[0]) < 1e-12
    assert abs(first[1] - second[1]) > 1e-3, (first, second)


def _bposd_compute_lz(hx, hz):
    # compute_lz of bposd css.py: the pivot rows of [hz; ker(hx)] past rank(hz)
    from syndrome_leakage.css import _kernel, _rank
    stack = np.vstack([hz, _kernel(hx)])
    rows, picked = [], []
    for i, r in enumerate(stack):
        if _rank(np.array(rows + [r])) > len(rows):
            rows.append(r)
            picked.append(i)
    return stack[picked[_rank(hz):]]


def test_supplied_z_logicals_are_kept_and_x_is_repaired():
    from syndrome_leakage import load

    class Obj:
        pass

    hx, hz = toric(2)
    o = Obj()
    o.hx, o.hz = hx, hz
    o.lz = _bposd_compute_lz(hx, hz)
    lx = _bposd_compute_lz(hz, hx)
    o.lx = np.array([(lx[0] + lx[1]) % 2, lx[1]], np.uint8)          # valid, and not paired with lz
    assert not np.array_equal((o.lz @ o.lx.T) % 2, np.eye(2, dtype=np.uint8))
    c = load.from_css_code(o, "toric, supplied logicals")
    kept = ["".join("1" if ch == "Z" else "0" for ch in z) for z in c.zl_strs]
    assert kept == ["".join(map(str, r)) for r in o.lz]
    assert _symplectic(c.logical_pairs)
    plain = load.from_css_code(o, "toric, computed basis", use_logicals=False)
    assert plain.k == 2 and _symplectic(plain.logical_pairs)


def test_supplied_logicals_that_are_not_logicals_are_refused():
    hx, hz = toric(2)
    for bad in (np.zeros((2, 8), np.uint8), np.array([hz[0], hz[1]], np.uint8)):
        try:
            css_strings(hx, hz, "bad", lz=bad)
        except AssertionError:
            continue
        raise AssertionError("stabilizer rows were accepted as logical operators")
