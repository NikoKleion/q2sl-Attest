# syndrome probabilities of a CSS code under a transversal Z-rotation (Hu, Liang and Calderbank)
import itertools
import math

import numpy as np


def is_css(code):
    return all(set(g) <= set("IX") or set(g) <= set("IZ") for g in code.stab_strings)


def generator_coefficients(code, theta):
    # A[(syndrome key, g)]: sum of f(z) over the coset of C1-perp with that X-syndrome and parity g against X_L
    if not is_css(code):
        raise ValueError(f"{code.name} is not a CSS code")
    n = code.n
    xl = [q for q, ch in enumerate(code.xl_str) if ch in "XY"]
    xsup = [[q for q, ch in enumerate(g) if ch == "X"] if set(g) <= set("IX") else None
            for g in code.stab_strings]
    c, s = math.cos(theta / 2), -1j * math.sin(theta / 2)
    A = {}
    for z in itertools.product((0, 1), repeat=n):
        w = sum(z)
        key = tuple(0 if sup is None else sum(z[q] for q in sup) % 2 for sup in xsup)
        g = sum(z[q] for q in xl) % 2
        A[(key, g)] = A.get((key, g), 0) + c ** (n - w) * s ** w
    return A


def syndrome_probabilities(code, theta, z_expectation):
    # eq. 91 with k = 1: p_mu = |A_mu,0|^2 + |A_mu,1|^2 + 2 Re(A_mu,0 conj(A_mu,1)) <Z_L>
    A = generator_coefficients(code, theta)
    out = {}
    for key in {k for k, _g in A}:
        a0, a1 = A.get((key, 0), 0), A.get((key, 1), 0)
        out[key] = float(abs(a0) ** 2 + abs(a1) ** 2 + 2 * (a0 * np.conj(a1)).real * z_expectation)
    return out
