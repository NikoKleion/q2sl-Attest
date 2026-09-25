# logical error after one round of syndrome measurement and single-qubit recovery
import itertools
import math

import numpy as np

from .core import per_qubit_channels
from .expectations import apply_pauli_string, logical_state_vector


def _anticommutes(a, b):
    return sum(1 for x, y in zip(a, b) if x != "I" and y != "I" and x != y) % 2 == 1


def syndrome_of(code, pauli):
    """Syndrome bits of a Pauli string, bit j for generator j."""
    return tuple(int(_anticommutes(pauli, g)) for g in code.stab_strings)


def distance(code, max_weight=4, max_coset=2 ** 22, max_search=2_000_000):
    """Smallest weight of a Pauli that commutes with every generator and anticommutes with Z_L or X_L."""
    from .css import code_distance, css_matrices
    try:
        Hx, Hz = css_matrices(code.stab_strings)
    except ValueError:
        Hx = None
    if Hx is not None:
        from .css import _rank
        if 2 ** max(_rank(Hx), _rank(Hz)) <= max_coset:
            return code_distance(Hx, Hz)
    cost = sum(math.comb(code.n, w) * 3 ** w for w in range(1, max_weight + 1))
    if cost > max_search:
        return None
    for w in range(1, max_weight + 1):
        for qs in itertools.combinations(range(code.n), w):
            for ps in itertools.product("XYZ", repeat=w):
                e = ["I"] * code.n
                for q, p in zip(qs, ps):
                    e[q] = p
                e = "".join(e)
                if not any(syndrome_of(code, e)) and (_anticommutes(e, code.zl_str) or _anticommutes(e, code.xl_str)):
                    return w
    return None


def recovery_strings(code):
    """Syndrome bits -> recovery Pauli string, for the trivial syndrome and every single-qubit error."""
    table = {tuple(0 for _ in code.stab_strings): "I" * code.n}
    for q in range(code.n):
        for pauli in "XYZ":
            e = "I" * q + pauli + "I" * (code.n - q - 1)
            table.setdefault(syndrome_of(code, e), e)
    return table


def evolve(psi, kraus, n):
    """E(|psi><psi|) for one single-qubit channel on every qubit, or one channel per qubit."""
    chans = per_qubit_channels(kraus, n)
    T = np.outer(psi, psi.conj()).reshape([2] * n + [2] * n)
    for q in range(n):
        out = np.zeros_like(T)
        for K in chans[q]:
            A = np.moveaxis(np.tensordot(K, T, axes=([1], [q])), 0, q)
            out += np.moveaxis(np.tensordot(K.conj(), A, axes=([1], [n + q])), 0, n + q)
        T = out
    return T.reshape(2 ** n, 2 ** n)


def table_decoder(code):
    """Syndrome bits -> the first single-qubit Pauli with that syndrome, identity when none has it."""
    table = recovery_strings(code)
    blank = "I" * code.n
    return lambda bits: table.get(tuple(bits), blank)


def css_decoder(code, error_rate=0.05, max_iter=48, osd_order=7):
    """Syndrome bits -> a correction from ldpc's BP+OSD, X errors from the Z checks and Z from the X."""
    from ldpc import bposd_decoder
    from .css import css_matrices
    Hx, Hz = css_matrices(code.stab_strings)
    xj = [j for j, g in enumerate(code.stab_strings) if set(g) <= {"I", "X"} and set(g) != {"I"}]
    zj = [j for j, g in enumerate(code.stab_strings) if set(g) <= {"I", "Z"} and set(g) != {"I"}]
    from .css import _rank

    def mk(H):
        # ldpc caps the OSD order at the dimension of the classical code
        order = max(0, min(osd_order, code.n - _rank(H)))
        return bposd_decoder(H, channel_probs=[error_rate] * code.n, max_iter=max_iter,
                             bp_method="ms", osd_method="osd_cs", osd_order=order)
    dec_x, dec_z = (mk(Hz) if len(Hz) else None), (mk(Hx) if len(Hx) else None)

    def decode(bits):
        bits = tuple(bits)
        xs = np.zeros(code.n, int) if dec_x is None else np.asarray(
            dec_x.decode(np.array([bits[j] for j in zj], int)), int)
        zs = np.zeros(code.n, int) if dec_z is None else np.asarray(
            dec_z.decode(np.array([bits[j] for j in xj], int)), int)
        return "".join("Y" if x and z else "X" if x else "Z" if z else "I" for x, z in zip(xs, zs))

    return decode


def outcome(code, error, correction):
    """What a correction does to an error: corrected, a logical, or a syndrome it does not explain."""
    from .core import pmul
    resid = pmul(error, correction)
    if any(syndrome_of(code, resid)):
        return "detected"
    if _anticommutes(resid, code.zl_str) or _anticommutes(resid, code.xl_str):
        return "logical"
    return "corrected"


def recovered_fidelity(code, rho, psi, decoder=None):
    """<psi| rho' |psi> for a state rho before syndrome measurement."""
    m = len(code.stab_strings)
    dec = decoder if decoder is not None else table_decoder(code)
    f = 0.0
    for s in range(2 ** m):
        bits = tuple((s >> j) & 1 for j in range(m))
        r = dec(bits)
        if tuple(syndrome_of(code, r)) != bits:
            continue
        v = apply_pauli_string(psi, r)
        f += float(np.real(np.vdot(v, rho @ v)))
    return f


def logical_error(code, kraus, t=0.0, phi=0.0, decoder=None):
    """1 - <psi| rho' |psi> for cos(t/2)|0_L> + e^{i phi} sin(t/2)|1_L> under the channel on every qubit."""
    psi = logical_state_vector(code, t, phi)
    return 1.0 - recovered_fidelity(code, evolve(psi, kraus, code.n), psi, decoder)


STATES = (("0_L", 0.0, 0.0), ("1_L", math.pi, 0.0), ("+_L", math.pi / 2, 0.0))


def protection_row(code, kraus, decoder=None):
    """Logical error after recovery for |0_L>, |1_L> and |+_L>."""
    return {name: logical_error(code, kraus, t, phi, decoder) for name, t, phi in STATES}
