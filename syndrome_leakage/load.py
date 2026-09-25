# codes from files and from other packages' objects
import numpy as np

from .css import css_matrices, css_strings

PAULI = set("IXYZ")


def read_dense(path):
    """A 0/1 matrix from whitespace-separated text, blank lines and # comments ignored."""
    rows = []
    for line in open(path, encoding="utf-8"):
        line = line.split("#")[0].strip()
        if line:
            rows.append([int(x) for x in line.replace(",", " ").split()])
    return np.array(rows, np.uint8)


def read_alist(path):
    """A parity-check matrix from alist, the format of MacKay's LDPC tables."""
    nums = open(path, encoding="utf-8").read().split()
    pos = 0

    def take(k):
        nonlocal pos
        out = [int(x) for x in nums[pos:pos + k]]
        pos += k
        return out

    n, m = take(2)
    take(2)
    col_w = take(n)
    take(m)
    H = np.zeros((m, n), np.uint8)
    width = max(col_w) if col_w else 0
    for j in range(n):
        for r in take(width):
            if r:
                H[r - 1, j] = 1
    return H


def read_mtx(path):
    """A 0/1 matrix from MatrixMarket coordinate text, the format QDistRnd writes."""
    entries, shape = [], None
    for line in open(path, encoding="utf-8"):
        if line.startswith("%"):
            continue
        parts = line.split()
        if not parts:
            continue
        if shape is None:
            shape = (int(parts[0]), int(parts[1]))
            continue
        i, j = int(parts[0]) - 1, int(parts[1]) - 1
        v = float(parts[2]) if len(parts) > 2 else 1.0
        if int(round(v)) % 2:
            entries.append((i, j))
    H = np.zeros(shape, np.uint8)
    for i, j in entries:
        H[i, j] = 1
    return H


def read_matrix(path):
    """A check matrix from .alist, .mtx, .npy or 0/1 text, by file name."""
    p = str(path).lower()
    if p.endswith(".alist"):
        return read_alist(path)
    if p.endswith(".mtx"):
        return read_mtx(path)
    if p.endswith(".npy"):
        return np.asarray(np.load(path), np.uint8)
    return read_dense(path)


def read_strings(path):
    """Pauli strings, one generator per line."""
    out = []
    for line in open(path, encoding="utf-8"):
        line = line.split("#")[0].strip().upper()
        if line and set(line) <= PAULI:
            out.append(line)
    return out


def code_from_strings(strings, name="loaded"):
    """A strings-only code from CSS generators, each of them all X or all Z."""
    return css_strings(*css_matrices(list(strings)), name=name)


def load_code(path, hz_path=None, name=None):
    """A code from one file or from a pair of check-matrix files."""
    label = name or str(path).replace("\\", "/").split("/")[-1].split(".")[0]
    if hz_path is not None:
        return css_strings(read_matrix(path), read_matrix(hz_path), label)
    if str(path).lower().endswith(".npz"):
        z = np.load(path)
        keys = {k.lower(): k for k in z.files}
        assert "hx" in keys and "hz" in keys, f"expected keys hx and hz, got {z.files}"
        return css_strings(np.asarray(z[keys["hx"]], np.uint8), np.asarray(z[keys["hz"]], np.uint8), label)
    text = open(path, encoding="utf-8").read()
    if any(set(line.strip().upper()) <= PAULI and line.strip() for line in text.splitlines()):
        return code_from_strings(read_strings(path), label)
    blocks = [b for b in text.split("\n\n") if b.strip()]
    assert len(blocks) == 2, f"expected two blocks separated by a blank line, got {len(blocks)}"
    mats = []
    for b in blocks:
        mats.append(np.array([[int(x) for x in line.replace(",", " ").split()]
                              for line in b.strip().splitlines() if line.strip()], np.uint8))
    return css_strings(mats[0], mats[1], label)


def _dense(m):
    return np.asarray(m.toarray() if hasattr(m, "toarray") else m, np.uint8) % 2


def from_qldpc(code, name=None, use_logicals=True):
    """A code from a qLDPC CSSCode."""
    field = getattr(code, "field", None)
    if field is not None and getattr(field, "order", 2) != 2:
        raise ValueError(f"{type(code).__name__} is over GF({field.order}); only qubit codes are read")
    if not hasattr(code, "matrix_x"):
        raise ValueError(f"{type(code).__name__} is not a CSSCode; only CSS codes are read")
    Hx = np.asarray(code.matrix_x, dtype=np.uint8) % 2
    Hz = np.asarray(code.matrix_z, dtype=np.uint8) % 2
    lx = lz = None
    if use_logicals:
        from qldpc.objects import Pauli
        lx = np.asarray(code.get_logical_ops(Pauli.X), dtype=np.uint8) % 2
        lz = np.asarray(code.get_logical_ops(Pauli.Z), dtype=np.uint8) % 2
    return css_strings(Hx, Hz, name or type(code).__name__, lx=lx, lz=lz)


def from_css_code(obj, name=None, use_logicals=True):
    """A code from an object carrying hx and hz, as bposd's css_code does, or a qLDPC CSSCode."""
    if hasattr(obj, "matrix_x") and not hasattr(obj, "hx"):
        return from_qldpc(obj, name, use_logicals)
    lz = getattr(obj, "lz", None) if use_logicals else None
    lx = getattr(obj, "lx", None) if use_logicals else None
    return css_strings(_dense(obj.hx), _dense(obj.hz), name or type(obj).__name__,
                       lx=None if lx is None else _dense(lx), lz=None if lz is None else _dense(lz))


# qecdb.org: a catalogue of published stabilizer codes, served as Pauli strings

QECDB = "https://qecdb.org/api"


def _qecdb_get(path, timeout=30):
    import json
    import urllib.error
    import urllib.request
    req = urllib.request.Request(f"{QECDB}/{path}", headers={"User-Agent": "q2sl"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.load(r)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise LookupError(f"qecdb.org has no entry at {path}") from e
        raise RuntimeError(f"qecdb.org returned HTTP {e.code} for {path}") from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"could not reach qecdb.org: {e.reason}") from e


def _pauli_rows(strings, kind, n):
    if not strings:
        return np.zeros((0, n), np.uint8)
    return np.array([[1 if c == kind else 0 for c in s] for s in strings], np.uint8)


def code_from_qecdb_record(rec, name=None):
    """A code from one qecdb.org record, whose H field holds the stabilizers and L the logicals."""
    stabs = [s.lstrip("+-") for s in rec["H"].split()]
    n = len(stabs[0])
    mixed = [s for s in stabs if "Y" in s or ("X" in s and "Z" in s)]
    if mixed:
        raise ValueError(f"qecdb code {rec.get('_id', '?')} is not CSS; generator {mixed[0]} mixes X and Z")
    Hx = _pauli_rows([s for s in stabs if "X" in s], "X", n)
    Hz = _pauli_rows([s for s in stabs if "Z" in s], "Z", n)
    logs = [s.lstrip("+-") for s in rec.get("L", "").split()]
    lx = [s for s in logs if "X" in s and set(s) <= {"I", "X"}]
    lz = [s for s in logs if "Z" in s and set(s) <= {"I", "Z"}]
    k = int(rec.get("k", len(lz)))
    usable = len(lx) == len(lz) == k and len(lx) + len(lz) == len(logs)
    label = name or " ".join(p for p in (rec.get("name", ""), rec.get("desc", ""),
                                         f"(qecdb {rec.get('_id', '?')})") if p)
    return css_strings(Hx, Hz, label, lx=_pauli_rows(lx, "X", n) if usable else None,
                       lz=_pauli_rows(lz, "Z", n) if usable else None)


def from_qecdb(code_id, timeout=30):
    """The code with this qecdb.org id, fetched and built as a strings-only CSS code."""
    return code_from_qecdb_record(_qecdb_get(f"codes/{code_id}", timeout), None)


def qecdb_search(n=None, k=None, d=None, family=None, limit=25, page=0, timeout=30):
    """Index records of qecdb.org CSS codes."""
    from urllib.parse import urlencode
    params = {"sort": "n", "dir": "asc", "page": page, "pageSize": limit}
    for key, value in (("n", n), ("k", k), ("d", d)):
        if value is not None:
            lo, hi = value if isinstance(value, (tuple, list)) else (value, value)
            params[f"{key}Min"], params[f"{key}Max"] = lo, hi
    if family:
        params["family"] = family
    rows = _qecdb_get("codes?" + urlencode(params), timeout)
    return [{f: r.get(f) for f in ("_id", "name", "n", "k", "d", "desc")} for r in rows.get("items", [])]
