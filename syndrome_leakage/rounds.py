# repeated extraction with recovery, as one channel on the logical qubit
import numpy as np

from .core import op, per_qubit_channels, pmul
from .expectations import logical_state_vector
from .protection import evolve, syndrome_of


def coset_decoder(code):
    """Syndrome bits -> a Pauli with that syndrome, for every syndrome, by GF(2) elimination."""
    m = len(code.stab_strings)
    basis, seen = [], {tuple([0] * m): "I" * code.n}
    for q in range(code.n):
        for p in "XYZ":
            e = "I" * q + p + "I" * (code.n - q - 1)
            s = tuple(syndrome_of(code, e))
            if s not in seen:
                for prev_s, prev_e in list(seen.items()):
                    combo = tuple(a ^ b for a, b in zip(s, prev_s))
                    if combo not in seen:
                        seen[combo] = pmul(e, prev_e)
                seen[s] = e
                basis.append((s, e))
    return lambda bits: seen.get(tuple(bits), "I" * code.n)


def is_complete(code, decoder):
    """Does the decoder return a correction carrying the measured syndrome, for every syndrome?"""
    m = len(code.stab_strings)
    for s in range(2 ** m):
        bits = tuple((s >> j) & 1 for j in range(m))
        if tuple(syndrome_of(code, decoder(bits))) != bits:
            return False
    return True


def pauli_action(s):
    """(permutation, phase) with P|j> = phase[j] |perm[j]>, for a Pauli string."""
    n = len(s)
    j = np.arange(2 ** n, dtype=np.int64)
    xmask = sum(1 << (n - 1 - q) for q, ch in enumerate(s) if ch in "XY")
    zmask = sum(1 << (n - 1 - q) for q, ch in enumerate(s) if ch in "ZY")
    ny = sum(1 for ch in s if ch == "Y")
    bits = np.zeros_like(j)
    m = zmask
    while m:
        low = m & -m
        bits ^= (j & low) // low
        m ^= low
    phase = np.where(bits & 1, -1.0, 1.0).astype(complex) * (1j ** ny)
    return j ^ xmask, phase


def left(rho, act):
    # P rho
    perm, phase = act
    return phase[perm][:, None] * rho[perm, :]


def right(rho, act):
    # rho P
    perm, phase = act
    return rho[:, perm] * phase[None, :]


def _logical_ops(code):
    XL, ZL = op(code.xl_str), op(code.zl_str)
    return XL, 1j * XL @ ZL, ZL


def code_projector(code):
    """The projector onto the code space, stabilizers and every logical Z at +1 left free."""
    I = np.eye(2 ** code.n, dtype=complex)
    P = I
    for g in code.stab_strings:
        P = P @ ((I + op(g)) / 2)
    for z in getattr(code, "zl_strs", [code.zl_str]):
        if z != code.zl_str:
            P = P @ ((I + op(z)) / 2)
    return P


def round_map(code, kraus, decoder=None):
    """(transfer matrix, syndrome map) for one round: noise, ideal measurement, then recovery."""
    from .protection import table_decoder
    dec = decoder if decoder is not None else coset_decoder(code)
    assert is_complete(code, dec), "this path needs a decoder that explains every syndrome"
    n, m = code.n, len(code.stab_strings)
    XL, YL, ZL = _logical_ops(code)
    P0 = code_projector(code)
    basis = [P0 / 2, XL @ P0 / 2, YL @ P0 / 2, ZL @ P0 / 2]
    corrections = [pauli_action(dec(tuple((s >> j) & 1 for j in range(m)))) for s in range(2 ** m)]
    factors = [pauli_action(g) for g in code.stab_strings]
    reads = [pauli_action(code.xl_str), pauli_action(_y_string(code)), pauli_action(code.zl_str)]
    yphase = _y_phase(code)
    T = np.zeros((4, 4))
    A = np.zeros((2 ** m, 4))
    chans = per_qubit_channels(kraus, n)
    for col, rho in enumerate(basis):
        out = _apply_channel(rho, chans, n)
        for s in range(2 ** m):
            branch = out
            for j, act in enumerate(factors):
                sign = -1.0 if (s >> j) & 1 else 1.0
                branch = (branch + sign * (left(branch, act) + right(branch, act))
                          + left(right(branch, act), act)) / 4
            A[s, col] = float(np.real(np.trace(branch)))
            R = corrections[s]
            fixed = left(right(branch, R), R)
            T[0, col] += float(np.real(np.trace(fixed)))
            for row, act in enumerate(reads, start=1):
                v = _trace_with(fixed, act) * (yphase if row == 2 else 1.0)
                T[row, col] += float(np.real(v))
    return T, A


def _trace_with(rho, act):
    # Tr(P rho) for a Pauli action
    perm, phase = act
    return complex(np.sum(phase[perm] * rho[perm, np.arange(rho.shape[0])]))


def _y_string(code):
    # the Pauli string of X_L Z_L, whose phase i is carried separately
    return pmul(code.xl_str, code.zl_str)


def _y_phase(code):
    # i X_L Z_L = phase * (the string X_L Z_L)
    from .expectations import pauli_product
    ph, _s = pauli_product(code.xl_str, code.zl_str)
    return 1j * ph



def _apply_channel(rho, chans, n):
    T = rho.reshape([2] * n + [2] * n)
    for q in range(n):
        acc = np.zeros_like(T)
        for K in chans[q]:
            B = np.moveaxis(np.tensordot(K, T, axes=([1], [q])), 0, q)
            acc += np.moveaxis(np.tensordot(K.conj(), B, axes=([1], [n + q])), 0, n + q)
        T = acc
    return T.reshape(2 ** n, 2 ** n)


def iterate(code, kraus, rounds=100, decoder=None, states=((0.0, 0.0), (np.pi, 0.0))):
    """Per-round syndrome distribution for two logical states, held through recovery each round."""
    from .core import tvd
    from .eavesdrop import chernoff_exponent
    T, A = round_map(code, kraus, decoder)
    XL, YL, ZL = _logical_ops(code)
    P0 = code_projector(code)
    vecs = []
    for t, phi in states:
        psi = logical_state_vector(code, t, phi)
        rho = np.outer(psi, psi.conj())
        vecs.append(np.array([1.0] + [float(np.real(np.trace(L @ rho))) for L in (XL, YL, ZL)]))
    out, total = [], 0.0
    for r in range(rounds):
        dists = []
        for i, v in enumerate(vecs):
            d = A @ v
            d = np.clip(d, 0.0, None)
            dists.append(d / max(d.sum(), 1e-300))
            vecs[i] = T @ v
        C, _s = chernoff_exponent(dists[0], dists[1])
        total += C
        out.append({"round": r + 1, "tvd": float(tvd(dists[0], dists[1])), "chernoff": C,
                    "cumulative": total})
    return out
