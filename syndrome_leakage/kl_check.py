# analytic leak order against the amplitude-damping population distance
import itertools

import numpy as np

from .analyze import analytic_leak


def ad_population_distance(code, tol=1e-9):
    # min |A| with <0_L| n_A |0_L> != <1_L| n_A |1_L>; None if no A separates the codewords.
    n = code.n
    xs = np.arange(code.dim)
    p0 = np.abs(code.V0) ** 2
    p1 = np.abs(code.V1) ** 2
    for w in range(1, n + 1):
        for A in itertools.combinations(range(n), w):
            mask = np.ones(code.dim, bool)
            for q in A:
                mask &= ((xs >> (n - 1 - q)) & 1).astype(bool)
            if abs(float(p0[mask].sum() - p1[mask].sum())) > tol:
                return w
    return None


def correspondence_report():
    from .codes import STANDARD
    L = ["Analytic leak order and amplitude-damping population distance.",
         "The AD distance is the smallest weight of a number operator that separates the codewords.",
         "",
         "  code                 analytic order   AD population distance   match",
         "  " + "-" * 72]
    all_match = True
    for name, maker in STANDARD.items():
        code = maker()
        leaks, order = analytic_leak(code)
        a = order if leaks else None
        r_kl = ad_population_distance(code)
        match = (a == r_kl)
        all_match = all_match and match
        L.append(f"  {code.name:<18} {str(a):>15}   {str(r_kl):>22}   {'yes' if match else 'NO'}")
    L += ["", f"  analytic order equals the AD population distance on every code: {all_match}"]
    return "\n".join(L)


if __name__ == "__main__":
    print(correspondence_report())
