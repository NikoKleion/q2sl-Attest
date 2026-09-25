# codes built by qLDPC read into this package: same leak as the package's own constructions, and
# qLDPC's logical basis kept as it is
from _harness import needs
from syndrome_leakage import expectations as ex
from syndrome_leakage import load
from syndrome_leakage.analyze import analytic_leak
from syndrome_leakage.channels import amplitude_damping


def _codes():
    needs("qldpc")
    from qldpc import codes
    return codes


def test_surface_code_matches_the_package_surface_code():
    codes = _codes()
    code = load.from_qldpc(codes.SurfaceCode(3))
    assert (code.n, code.k) == (9, 1)
    assert analytic_leak(code)[1] == 3
    leak = ex.population_leak(code, amplitude_damping(0.2))[0]
    ours = ex.population_leak(ex.surface_code_3(), amplitude_damping(0.2))[0]
    assert abs(leak - ours) < 1e-12, (leak, ours)


def test_steane_is_protected():
    codes = _codes()
    code = load.from_qldpc(codes.SteaneCode())
    assert analytic_leak(code) == (False, None)
    assert ex.population_leak(code, amplitude_damping(0.2))[0] < 1e-12


def test_the_qldpc_logical_basis_is_kept():
    codes = _codes()
    import numpy as np
    from qldpc.objects import Pauli
    q = codes.ToricCode(4)
    code = load.from_qldpc(q)
    assert code.k == 2
    lz = np.asarray(q.get_logical_ops(Pauli.Z), dtype=int)
    for i, (z, _x) in enumerate(code.logical_pairs):
        assert [1 if c == "Z" else 0 for c in z] == lz[i].tolist()


def test_from_css_code_recognises_a_qldpc_object():
    codes = _codes()
    code = load.from_css_code(codes.SurfaceCode(3))
    assert (code.n, code.k) == (9, 1)
