# codes read from files, in each format the loader accepts
import os
import tempfile

import numpy as np

from syndrome_leakage import expectations as ex
from syndrome_leakage import load


def _files(tmp):
    Hx, Hz = ex.rotated_surface_code(3)
    dense = lambda H: "\n".join(" ".join(str(int(v)) for v in row) for row in H)
    open(os.path.join(tmp, "hx.txt"), "w").write(dense(Hx) + "\n")
    open(os.path.join(tmp, "hz.txt"), "w").write(dense(Hz) + "\n")
    open(os.path.join(tmp, "pair.txt"), "w").write(dense(Hx) + "\n\n" + dense(Hz) + "\n")
    np.savez(os.path.join(tmp, "code.npz"), hx=Hx, hz=Hz)
    open(os.path.join(tmp, "gens.txt"), "w").write("\n".join(ex.surface_code_3().stab_strings) + "\n")
    m, n = Hz.shape
    cols = [list(np.flatnonzero(Hz[:, j]) + 1) for j in range(n)]
    rows = [list(np.flatnonzero(Hz[i, :]) + 1) for i in range(m)]
    cw, rw = max(len(c) for c in cols), max(len(r) for r in rows)
    lines = [f"{n} {m}", f"{cw} {rw}", " ".join(str(len(c)) for c in cols), " ".join(str(len(r)) for r in rows)]
    lines += [" ".join(map(str, c + [0] * (cw - len(c)))) for c in cols]
    lines += [" ".join(map(str, r + [0] * (rw - len(r)))) for r in rows]
    open(os.path.join(tmp, "hz.alist"), "w").write("\n".join(lines) + "\n")
    ij = np.argwhere(Hz)
    mtx = ["%%MatrixMarket matrix coordinate integer general", f"{m} {n} {len(ij)}"]
    mtx += [f"{i + 1} {j + 1} 1" for i, j in ij]
    open(os.path.join(tmp, "hz.mtx"), "w").write("\n".join(mtx) + "\n")
    return Hx, Hz


def test_every_format_gives_the_same_code():
    ref = sorted(ex.surface_code_3().stab_strings)
    with tempfile.TemporaryDirectory() as tmp:
        _files(tmp)
        j = lambda f: os.path.join(tmp, f)
        codes = [load.load_code(j("hx.txt"), j("hz.txt")),
                 load.load_code(j("pair.txt")),
                 load.load_code(j("code.npz")),
                 load.load_code(j("gens.txt")),
                 load.load_code(j("hx.txt"), j("hz.alist")),
                 load.load_code(j("hx.txt"), j("hz.mtx"))]
        for c in codes:
            assert c.n == 9 and len(c.stab_strings) == 8
            assert sorted(c.stab_strings) == ref


def test_matrix_readers_round_trip():
    with tempfile.TemporaryDirectory() as tmp:
        _Hx, Hz = _files(tmp)
        for name in ("hz.txt", "hz.alist", "hz.mtx"):
            H = load.read_matrix(os.path.join(tmp, name))
            assert H.shape == Hz.shape, name
            assert np.array_equal(H, Hz), name


def test_from_css_code_object():
    class Obj:
        hx, hz = ex.rotated_surface_code(3)

    c = load.from_css_code(Obj(), "bposd-like")
    assert sorted(c.stab_strings) == sorted(ex.surface_code_3().stab_strings)


def test_a_file_that_is_not_a_code_is_refused():
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "junk.txt")
        open(p, "w").write("1 0 1\n")
        try:
            load.load_code(p)
        except (AssertionError, ValueError, IndexError):
            return
        raise AssertionError("a single block should not load as a CSS pair")
