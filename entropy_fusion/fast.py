# entropy_fusion.fast: closed-form evaluate and shapley
import numpy as np

from .experts import HardwareExpert, PatternExpert, FieldFormatExpert, RevealExpert


def generating_marginals(payload, samples=200000, seed=12345):
    # per-position marginal g_i[d] = P(instance[i]==d) from marginals() or a large sample
    m = getattr(payload, "marginals", None)
    if callable(m):
        g = np.asarray(m(), float)
        return g / g.sum(axis=1, keepdims=True)
    rng = np.random.default_rng(seed)
    counts = np.full((payload.n, payload.A), 1e-9)
    for _ in range(samples):
        x = payload.generate(rng)
        counts[np.arange(payload.n), x.astype(int)] += 1.0
    return counts / counts.sum(axis=1, keepdims=True)


def _classify(payload, experts):
    # base product B and at most one HardwareExpert; None if unsupported
    B = np.ones((payload.n, payload.A))
    hw = None
    dummy = np.zeros(payload.n, int); noknown = np.zeros(payload.n, bool)
    for e in experts:
        if isinstance(e, HardwareExpert):
            if hw is not None:
                return None
            hw = e
        elif isinstance(e, (PatternExpert, FieldFormatExpert, RevealExpert)):
            B = B * e.likelihood(payload, dummy, noknown)
        else:
            return None
    return B, hw


def _norm(v):
    s = v.sum()
    return v / s if s > 0 else np.full_like(v, 1.0 / len(v))


def evaluate_fast(payload, experts, reveals=0, marg=None):
    # exact (residual_bits, accuracy); None if the closed form does not apply
    if getattr(payload, "checksum", None) is not None:
        return None
    cl = _classify(payload, experts)
    if cl is None:
        return None
    B, hw = cl
    g = generating_marginals(payload) if marg is None else marg
    n, A = payload.n, payload.A
    p_rev = reveals / n
    reads = hw._reads(n) if hw is not None else set()
    p_read = hw.p if hw is not None else None

    tot_bits = 0.0; acc_terms = []
    for i in range(n):
        gi = g[i]
        if hw is not None and i in reads:
            bits_obs = np.empty(A); a_obs = np.empty(A, int)
            for obs in range(A):
                h = np.full(A, (1.0 - p_read) / (A - 1)); h[obs] = p_read
                post = _norm(B[i] * h)
                bits_obs[obs] = -np.log2(max(post.max(), 1e-12))
                a_obs[obs] = int(np.argmax(post))
            p_obs = np.zeros(A); e_acc = 0.0
            for t in range(A):
                pch = np.full(A, (1.0 - p_read) / (A - 1)); pch[t] = p_read
                p_obs += gi[t] * pch
                e_acc += gi[t] * float(np.sum(pch * (a_obs == t)))
            e_bits = float(np.sum(p_obs * bits_obs))
        else:
            post = _norm(B[i])
            a = int(np.argmax(post))
            e_bits = -np.log2(max(post.max(), 1e-12))
            e_acc = float(gi[a])
        tot_bits += (1.0 - p_rev) * e_bits
        acc_terms.append(p_rev * 1.0 + (1.0 - p_rev) * e_acc)
    return float(tot_bits), float(np.mean(acc_terms))


def shapley_fast(payload, experts, reveals=0, marg=None):
    # exact Shapley attribution reusing evaluate_fast; None if the closed form does not apply
    import itertools, math
    if getattr(payload, "checksum", None) is not None:
        return None
    g = generating_marginals(payload) if marg is None else marg
    k = len(experts)
    base = evaluate_fast(payload, [], reveals, g)
    if base is None:
        return None
    blind = base[0]
    val = {}
    for r in range(k + 1):
        for sub in itertools.combinations(range(k), r):
            ev = evaluate_fast(payload, [experts[i] for i in sub], reveals, g)
            if ev is None:
                return None
            val[sub] = blind - ev[0]
    phi = np.zeros(k)
    for i in range(k):
        for r in range(k):
            for sub in itertools.combinations([j for j in range(k) if j != i], r):
                w = math.factorial(r) * math.factorial(k - r - 1) / math.factorial(k)
                phi[i] += w * (val[tuple(sorted(sub + (i,)))] - val[sub])
    return phi


def main():
    import time
    from .payloads import PhoneNumber, DeviceErrorString
    from .experts import PatternExpert, HardwareExpert, FieldFormatExpert, RevealExpert
    from .core import evaluate

    for pay in (PhoneNumber(), DeviceErrorString()):
        experts = [FieldFormatExpert(), PatternExpert(pay), HardwareExpert(reliability=0.85), RevealExpert()]
        print(f"# {pay.name}")
        print(f"#   {'reveals':>7}  {'MC bits':>9}  {'fast bits':>9}  {'MC acc':>7}  {'fast acc':>8}  {'|dbits|':>8}")
        for r in (0, 2, 5):
            t0 = time.time(); mb, ma = evaluate(pay, experts, reveals=r, samples=4000); t_mc = time.time() - t0
            t0 = time.time(); fb, fa = evaluate_fast(pay, experts, reveals=r); t_fast = time.time() - t0
            print(f"#   {r:>7}  {mb:>9.3f}  {fb:>9.3f}  {ma:>7.3f}  {fa:>8.3f}  {abs(mb-fb):>8.3f}"
                  f"   (MC {t_mc*1e3:.0f} ms vs fast {t_fast*1e3:.1f} ms)")
        print()


if __name__ == "__main__":
    main()
