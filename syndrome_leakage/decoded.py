# the leak through a decoder's output: the correction it returns and coarser views of that correction
import itertools

import numpy as np

from .core import pmul, tvd
from .expectations import population_leak, walsh_hadamard
from .protection import _anticommutes, syndrome_of

VIEWS = ("correction", "support", "weight", "acted", "frame_z", "frame_x")


def _index(bits):
    return sum(int(b) << j for j, b in enumerate(bits))


def min_weight_decoder(code):
    """Syndrome bits -> a Pauli of least weight with that syndrome, the first found in qubit and X, Y, Z order."""
    n, m = code.n, len(code.stab_strings)
    single = {(q, p): _index(syndrome_of(code, "I" * q + p + "I" * (n - q - 1))) for q in range(n) for p in "XYZ"}
    table = {0: ()}
    for w in range(1, n + 1):
        if len(table) == 2 ** m:
            break
        for qs in itertools.combinations(range(n), w):
            for ps in itertools.product("XYZ", repeat=w):
                s = 0
                for q, p in zip(qs, ps):
                    s ^= single[(q, p)]
                if s not in table:
                    table[s] = tuple(zip(qs, ps))

    def decode(bits):
        e = ["I"] * n
        for q, p in table.get(_index(bits), ()):
            e[q] = p
        return "".join(e)

    return decode


def view_labels(code, decoder, checks=None):
    """One integer per record for each view of the correction; a record holds the selected generators, the rest zero."""
    m = len(code.stab_strings)
    idx = list(range(m)) if checks is None else list(checks)
    out, seen = {v: [] for v in VIEWS}, {}
    for r in range(2 ** len(idx)):
        bits = [0] * m
        for j, g in enumerate(idx):
            bits[g] = (r >> j) & 1
        c = decoder(tuple(bits))
        out["correction"].append(seen.setdefault(c, len(seen)))
        out["support"].append(sum(1 << q for q, ch in enumerate(c) if ch != "I"))
        out["weight"].append(sum(ch != "I" for ch in c))
        out["acted"].append(int(c != "I" * code.n))
        out["frame_z"].append(int(_anticommutes(c, code.zl_str)))
        out["frame_x"].append(int(_anticommutes(c, code.xl_str)))
    return {v: np.array(x, np.int64) for v, x in out.items()}


def push(d, labels):
    """The distribution of a view: the weights of d summed over the records that share a label."""
    return np.bincount(labels, weights=np.asarray(d, float), minlength=int(labels.max()) + 1)


def transform(labels):
    """A view as the transform argument of estimate.py."""
    return lambda d: push(d, labels)


def view_leaks(d0, d1, labels):
    """Distance between two record distributions, and between their images under each view."""
    out = {"syndrome": tvd(d0, d1)}
    out.update({v: tvd(push(d0, lab), push(d1, lab)) for v, lab in labels.items()})
    return out


def decoded_leak(code, kraus, decoder=None, readout=None):
    """view_leaks for |0_L> against |1_L> after one channel on every qubit."""
    _, d0, d1 = population_leak(code, kraus, readout=readout)
    return view_leaks(d0, d1, view_labels(code, decoder or min_weight_decoder(code)))


def frame_spectrum(d0, d1, frame):
    """Leak of the frame bit for every representative of the logical: entry a is the logical times S^a."""
    g = (np.asarray(d0, float) - np.asarray(d1, float)) * (1.0 - 2.0 * np.asarray(frame))
    return 0.5 * np.abs(np.real(walsh_hadamard(g)))


def representative(code, logical, a):
    """The logical string multiplied by the generators in the mask a."""
    out = logical
    for j, g in enumerate(code.stab_strings):
        if (a >> j) & 1:
            out = pmul(out, g)
    return out


def shor_closed_form(gammas, blocks=((0, 1, 2), (3, 4, 5), (6, 7, 8))):
    """Leak of the Shor code between |0_L> and |1_L> under amplitude damping of strength gammas[q] on qubit q."""
    out = 1.0
    for b in blocks:
        out *= 1.0 - np.prod([1.0 - gammas[q] for q in b]) - np.prod([gammas[q] for q in b])
    return float(out)


def type_masks(code, letter):
    """Masks over the generators written with I and one Pauli letter only."""
    keep = [j for j, g in enumerate(code.stab_strings) if set(g) <= {"I", letter}]
    return [sum(1 << keep[i] for i in range(len(keep)) if (b >> i) & 1) for b in range(2 ** len(keep))]
