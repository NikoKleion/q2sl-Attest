# syndrome_leakage.approx_qec: Knill-Laflamme matrices of a code under a single-qubit channel
import numpy as np

MAX_BRANCHES = 1 << 16


def _apply_1q_vec(v, K, q, n):
    # 2x2 operator on qubit q of an n-qubit vector; qubit 0 is the most significant axis
    t = np.tensordot(K, v.reshape([2] * n), axes=([1], [q]))
    return np.moveaxis(t, 0, q).reshape(-1)


def branch_matrices(code, kraus):
    # for every Kraus product A_e, the 2x2 matrix <V_i| A_e^dagger A_e |V_j> on the codeword basis
    n = code.n
    if len(kraus) ** n > MAX_BRANCHES:
        raise ValueError(f"{len(kraus)}^{n} Kraus products exceeds {MAX_BRANCHES}")
    branches = [((), code.V0, code.V1)]
    for q in range(n):
        nxt = []
        for e, w0, w1 in branches:
            for k, K in enumerate(kraus):
                nxt.append((e + (k,), _apply_1q_vec(w0, K, q, n), _apply_1q_vec(w1, K, q, n)))
        branches = nxt
    out = {}
    for e, w0, w1 in branches:
        W = np.stack([w0, w1], axis=1)
        out[e] = W.conj().T @ W
    return out


def branch_spread(code, kraus):
    # largest eigenvalue gap of P_C A_e^dagger A_e P_C over branches, p_n (1 - lambda_n) in Leung et al.
    best, arg = 0.0, None
    for e, M in branch_matrices(code, kraus).items():
        ev = np.linalg.eigvalsh(M)
        gap = float(ev[-1] - ev[0])
        if gap > best:
            best, arg = gap, e
    return best, arg


def syndrome_blocks(code, kraus):
    # per syndrome s, B_s[i, j] = sum_e <V_i| A_e^dagger Pi_s A_e |V_j> = Tr(Pi_s E(|V_j><V_i|))
    V = (code.V0, code.V1)
    evolved = {(i, j): code.apply(np.outer(V[j], V[i].conj()), kraus) for i in (0, 1) for j in (0, 1)}
    blocks = {}
    for s, P in code.PROJ.items():
        B = np.empty((2, 2), dtype=complex)
        for (i, j), R in evolved.items():
            B[i, j] = np.trace(P @ R)
        blocks[s] = B
    return blocks


def population_split(blocks):
    # total variation distance between the two diagonal entries across syndromes
    return 0.5 * float(sum(abs(B[0, 0] - B[1, 1]) for B in blocks.values()))


def coherence_by_syndrome(blocks):
    # |B_s[0, 1]| for each syndrome
    return {s: float(abs(B[0, 1])) for s, B in blocks.items()}


def fitted_order(fn, gammas=(0.005, 0.01, 0.02)):
    # log-log slope of fn(gamma) over small gamma; None when fn stays at machine zero
    ys = [fn(g) for g in gammas]
    if max(ys) < 1e-12:
        return None
    return float(np.polyfit(np.log(gammas), np.log(np.maximum(ys, 1e-300)), 1)[0])
