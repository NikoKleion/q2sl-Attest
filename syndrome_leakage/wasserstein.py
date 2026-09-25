# Wasserstein-1 syndrome distance under the Hamming ground metric
import math
import numpy as np

from .analyze import population_leak
from .core import tvd
from . import channels as _ch


def _marginal_lower_bound(keys, d0, d1):
    # sum of per-bit marginal gaps, a lower bound on the joint Hamming-W1
    K = np.array(keys)
    p0 = d0 @ K
    p1 = d1 @ K
    return float(np.abs(p0 - p1).sum())


def _exact_w1(keys, d0, d1):
    # exact Wasserstein-1 with Hamming ground metric via the transportation LP
    from scipy.optimize import linprog
    from scipy.sparse import lil_matrix
    K = np.array(keys); m = len(keys)
    C = (K[:, None, :] != K[None, :, :]).sum(axis=2).astype(float).ravel()
    A = lil_matrix((2 * m, m * m))
    for i in range(m):
        A[i, i * m:(i + 1) * m] = 1.0
        A[m + i, i::m] = 1.0
    b = np.concatenate([d0, d1])
    res = linprog(C, A_eq=A.tocsr(), b_eq=b, bounds=(0, None), method="highs")
    if not res.success:
        return None
    return float(res.fun)


def w1_leakage(code, channel="amplitude_damping", gamma=0.2, exact=True, max_exact_dim=512):
    # W1 leakage distance between the |0_L> and |1_L> syndrome distributions
    kraus = _ch.LIBRARY[channel](gamma) if isinstance(channel, str) else channel
    _, d0, d1 = population_leak(code, kraus)
    keys = list(code.PROJ.keys())
    lb = _marginal_lower_bound(keys, d0, d1)
    w1, is_exact = lb, False
    if exact and len(keys) <= max_exact_dim:
        try:
            v = _exact_w1(keys, d0, d1)
            if v is not None:
                w1, is_exact = v, True
        except Exception:
            pass
    return {"w1": w1, "exact": is_exact, "w1_lower_bound": lb,
            "tvd": tvd(d0, d1), "n_syndromes": len(keys),
            "channel": channel if isinstance(channel, str) else "custom", "gamma": gamma}


def main():
    from .codes import repetition, code_4_1_2, hamming_7
    print("# Wasserstein-1 syndrome distance (Hamming ground metric) and TVD")
    print()
    print(f"#   {'code':>16}  {'TVD':>10}  {'W1 (exact)':>12}  {'W1/TVD':>8}  {'syndromes':>9}")
    for maker in (repetition, code_4_1_2, hamming_7):
        c = maker()
        r = w1_leakage(c, "amplitude_damping", 0.2)
        ratio = r["w1"] / r["tvd"] if r["tvd"] > 1e-12 else float("nan")
        tag = "" if r["exact"] else "  (lower bound; no scipy)"
        print(f"#   {c.name:>16}  {r['tvd']:>10.3e}  {r['w1']:>12.3e}  {ratio:>8.2f}  {r['n_syndromes']:>9}{tag}")


if __name__ == "__main__":
    main()
