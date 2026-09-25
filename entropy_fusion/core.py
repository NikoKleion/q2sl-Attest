# entropy_fusion.core: fusion, scoring, evaluation and Shapley attribution
import math

import numpy as np


class Payload:
    # synthetic structured value; subclasses set n, A and implement generate
    name = "payload"
    n = 0
    A = 10

    def generate(self, rng):
        # return an int array of length n with entries in [0, A)
        raise NotImplementedError

    @property
    def field_alphabets(self):
        # per-position allowed values; default is the full alphabet
        return [np.arange(self.A) for _ in range(self.n)]

    checksum = None

    def blind_bits(self):
        # log2 of the number of structurally valid instances
        return float(sum(np.log2(len(v)) for v in self.field_alphabets))


class Expert:
    name = "expert"

    def likelihood(self, payload, instance, known):
        # return a (payload.n, payload.A) non-negative likelihood
        raise NotImplementedError


def _constrain_checksum_ref(post, F, mod):
    # reference sum-product, checked against the vectorized version
    n, A = post.shape
    P = np.clip(post, 1e-300, None)
    alpha = np.zeros((n + 1, mod)); alpha[0, 0] = 1.0
    for i in range(n):
        for d in range(A):
            alpha[i + 1] += np.roll(alpha[i], int(F[i, d])) * P[i, d]
        t = alpha[i + 1].sum()
        if t > 0:
            alpha[i + 1] /= t
    beta = np.zeros((n + 1, mod)); beta[n, 0] = 1.0
    for i in range(n - 1, -1, -1):
        for d in range(A):
            beta[i] += np.roll(beta[i + 1], -int(F[i, d])) * P[i, d]
        t = beta[i].sum()
        if t > 0:
            beta[i] /= t
    out = np.zeros_like(post)
    for i in range(n):
        for d in range(A):
            out[i, d] = P[i, d] * float(alpha[i] @ np.roll(beta[i + 1], -int(F[i, d])))
    s = out.sum(axis=1, keepdims=True)
    return np.where(s > 0, out / np.maximum(s, 1e-300), post)


def _constrain_checksum(post, F, mod):
    # exact sum-product for sum_i F[i, x_i] == 0 (mod), a circular convolution over Z_mod
    n, A = post.shape
    P = np.clip(post, 1e-300, None)
    R = np.arange(mod)
    fwd = (R[:, None] - R[None, :]) % mod
    bwd = (R[:, None] + R[None, :]) % mod
    C = np.stack([np.bincount(F[i], weights=P[i], minlength=mod) for i in range(n)])
    alpha = np.zeros((n + 1, mod)); alpha[0, 0] = 1.0
    for i in range(n):
        a = (alpha[i][fwd] * C[i][None, :]).sum(axis=1)
        t = a.sum(); alpha[i + 1] = a / t if t > 0 else a
    beta = np.zeros((n + 1, mod)); beta[n, 0] = 1.0
    for i in range(n - 1, -1, -1):
        b = (beta[i + 1][bwd] * C[i][None, :]).sum(axis=1)
        t = b.sum(); beta[i] = b / t if t > 0 else b
    out = np.zeros_like(post)
    for i in range(n):
        g = (alpha[i][None, :] * beta[i + 1][bwd]).sum(axis=1)
        out[i] = P[i] * g[F[i]]
    s = out.sum(axis=1, keepdims=True)
    return np.where(s > 0, out / np.maximum(s, 1e-300), post)


def _pin_known(post, instance, known):
    for i in range(len(known)):
        if known[i]:
            post[i] = 0.0; post[i, int(instance[i])] = 1.0
    return post


def fuse(payload, instance, known, experts):
    post = np.ones((payload.n, payload.A))
    for e in experts:
        post = post * e.likelihood(payload, instance, known)
    s = post.sum(axis=1, keepdims=True)
    post = np.where(s > 0, post / np.maximum(s, 1e-300), 1.0 / payload.A)
    if np.any(known):
        post = _pin_known(post, instance, known)
    if payload.checksum is not None:
        F, mod = payload.checksum
        post = _constrain_checksum(post, F, mod)
        if np.any(known):
            post = _pin_known(post, instance, known)
    return post


def residual_bits(post):
    # sum over positions of -log2(max posterior)
    return float(np.sum(-np.log2(np.clip(post.max(axis=1), 1e-12, 1.0))))


def evaluate(payload, experts, reveals=0, samples=300, seed=0):
    # Monte Carlo mean residual bits and accuracy, paired instances across expert subsets
    bits, acc = [], []
    for s in range(samples):
        rng = np.random.default_rng(seed * 100003 + s)
        inst = payload.generate(rng)
        known = np.zeros(payload.n, bool)
        if reveals:
            known[rng.permutation(payload.n)[:reveals]] = True
        post = fuse(payload, inst, known, experts)
        bits.append(residual_bits(post))
        acc.append(float(np.mean(np.argmax(post, axis=1) == inst)))
    return float(np.mean(bits)), float(np.mean(acc))


def shapley(payload, experts, reveals=0, samples=200, seed=1):
    # exact Shapley attribution over all 2^k subsets
    import itertools
    k = len(experts)
    blind, _ = evaluate(payload, [], reveals, samples, seed)
    val = {}
    for r in range(k + 1):
        for sub in itertools.combinations(range(k), r):
            val[sub] = blind - evaluate(payload, [experts[i] for i in sub], reveals, samples, seed)[0]
    phi = np.zeros(k)
    for i in range(k):
        for r in range(k):
            for sub in itertools.combinations([j for j in range(k) if j != i], r):
                w = math.factorial(r) * math.factorial(k - r - 1) / math.factorial(k)
                phi[i] += w * (val[tuple(sorted(sub + (i,)))] - val[sub])
    return phi

