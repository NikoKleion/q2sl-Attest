# entropy_fusion.adaptive: fusion with the classical prior weighted by an estimated trust
import numpy as np

from .core import residual_bits


def _normrows(M):
    s = M.sum(axis=1, keepdims=True)
    return np.where(s > 0, M / np.maximum(s, 1e-300), 1.0 / M.shape[1])


def trust_weight(payload, prior_L, quantum_L, instance, known, smooth=True):
    # trust in [0, 1] for the classical prior
    n, A = prior_L.shape
    signals, weights = [], []
    rev = [i for i in range(n) if known[i]]
    if rev:
        cal = np.mean([prior_L[i, int(instance[i])] for i in rev])
        signals.append(float(np.clip((cal - 1.0 / A) / (1.0 - 1.0 / A), 0.0, 1.0)))
        weights.append(2.0)
    unrev = [i for i in range(n) if not known[i]]
    if unrev and smooth:
        Pc, Pq = _normrows(prior_L), _normrows(quantum_L)
        bc = np.sqrt(np.clip(Pc, 0, None) * np.clip(Pq, 0, None)).sum(axis=1)
        conf = (Pc.max(axis=1) - 1.0 / A) * (Pq.max(axis=1) - 1.0 / A)
        w = np.array([max(0.0, conf[i]) for i in unrev]); bcu = np.array([bc[i] for i in unrev])
        if w.sum() > 1e-9:
            f = 1.0 / np.sqrt(A)
            bc_bar = float(np.average(bcu, weights=w))
            signals.append(float(np.clip((bc_bar - f) / (1.0 - f), 0.0, 1.0)))
            weights.append(1.0)
    elif unrev and not smooth:
        conf = quantum_L.max(axis=1) > 1.5 / A
        idx = [i for i in unrev if conf[i]]
        if idx:
            agree = np.mean([np.argmax(prior_L[i]) == np.argmax(quantum_L[i]) for i in idx])
            signals.append(float(np.clip((agree - 1.0 / A) / (1.0 - 1.0 / A), 0.0, 1.0)))
            weights.append(1.0)
    if not signals:
        return 0.5
    return float(np.average(signals, weights=weights))


def adaptive_fuse(payload, instance, known, prior_expert, quantum_expert, other_experts=()):
    Pc = prior_expert.likelihood(payload, instance, known)
    Pq = quantum_expert.likelihood(payload, instance, known)
    w = trust_weight(payload, Pc, Pq, instance, known)
    Pc_t = np.clip(Pc, 1e-12, None) ** w
    Pc_t = Pc_t / Pc_t.sum(axis=1, keepdims=True)
    post = Pc_t * Pq
    for e in other_experts:
        post = post * e.likelihood(payload, instance, known)
    for i in range(payload.n):
        if known[i]:
            post[i] = 0.0; post[i, int(instance[i])] = 1.0
    s = post.sum(axis=1, keepdims=True)
    return np.where(s > 0, post / np.maximum(s, 1e-300), 1.0 / payload.A), w


def evaluate_adaptive(payload, prior_expert, quantum_expert, other_experts=(), reveals=0, samples=400, seed=0):
    bits, acc, ws = [], [], []
    for s in range(samples):
        rng = np.random.default_rng(seed * 100003 + s)
        inst = payload.generate(rng)
        known = np.zeros(payload.n, bool)
        if reveals:
            known[rng.permutation(payload.n)[:reveals]] = True
        post, w = adaptive_fuse(payload, inst, known, prior_expert, quantum_expert, other_experts)
        bits.append(residual_bits(post)); acc.append(float(np.mean(np.argmax(post, axis=1) == inst))); ws.append(w)
    return float(np.mean(bits)), float(np.mean(acc)), float(np.mean(ws))


def main():
    import os, sys
    from .registry import SpecPayload
    from .experts import PatternExpert, HardwareExpert, RevealExpert
    from .core import evaluate

    def structured_bits(n, c, seed):
        rng = np.random.default_rng(seed); q = 1.0 / (1.0 + np.exp(-c * rng.standard_normal(n)))
        return SpecPayload(f"bits(c={c})", n, 2, lambda r: (r.random(n) < q).astype(int))

    n = 16
    pay = structured_bits(n, 2.5, seed=0)
    hw = HardwareExpert(reliability=0.9); rev = RevealExpert()
    prior_good = PatternExpert(pay, train_samples=40000)
    prior_bad = PatternExpert(structured_bits(n, 2.5, seed=999), train_samples=40000, temperature=0.5)

    print("# accuracy: hardware only, static fusion with the prior, adaptive fusion, mean prior weight w\n")
    for label, prior in (("calibrated prior", prior_good), ("miscalibrated prior", prior_bad)):
        print(f"# {label}")
        print(f"#   {'reveals':>7}  {'hardware only':>13}  {'static prior+hw':>15}  {'adaptive':>9}  {'mean w':>7}")
        for r in (0, 2, 4, 8):
            _, a_hw = evaluate(pay, [hw, rev], reveals=r)
            _, a_st = evaluate(pay, [prior, hw, rev], reveals=r)
            _, a_ad, w = evaluate_adaptive(pay, prior, hw, (rev,), reveals=r)
            print(f"#   {r:>7}  {a_hw:>13.3f}  {a_st:>15.3f}  {a_ad:>9.3f}  {w:>7.2f}")
        print()


if __name__ == "__main__":
    main()
