# covert_channel: a classical bit carried by data-dependent noise, read from syndromes
import os, math
import numpy as np

L = 4
N = 2 * L * L
BASE_P = 0.03


def toric_Hz(L):
    n = 2 * L * L
    Hh = lambda r, c: (r % L) * L + (c % L)
    Vt = lambda r, c: L * L + (r % L) * L + (c % L)
    rows = []
    for r in range(L):
        for c in range(L):
            row = np.zeros(n, dtype=int)
            for idx in (Hh(r, c), Hh(r + 1, c), Vt(r, c), Vt(r, c + 1)):
                row[idx] ^= 1
            rows.append(row)
    return np.array(rows, dtype=int)


H = toric_Hz(L)
M = H.shape[0]
SUPP = [np.nonzero(H[j])[0] for j in range(M)]
Q = list(range(L))


def device_rates(seed=7, hot_frac=0.25, hot_mult=8.0, eta_hot=1.0):
    # per-qubit X error rates
    rng = np.random.default_rng(seed)
    hot = np.ones(N); hot[rng.random(N) < hot_frac] = hot_mult
    return np.clip(BASE_P * hot, 1e-4, 0.45)


def rates_for(base, b, alpha):
    # b=1 raises the noise on loop Q by (1 + alpha); alpha=0 is no data dependence
    r = base.copy()
    if b == 1:
        r[Q] = np.clip(r[Q] * (1.0 + alpha), 1e-4, 0.45)
    return r


def fire_prob(rates, obs):
    # P(check j fires) = (1 - prod_{i in support}(1 - 2 r_i)) / 2, for each observed check
    p = np.empty(len(obs))
    for k, j in enumerate(obs):
        p[k] = 0.5 * (1.0 - np.prod(1.0 - 2.0 * rates[SUPP[j]]))
    return np.clip(p, 1e-6, 1 - 1e-6)


def attack_accuracy(base, alpha, K, obs, trials, rng, trained=True):
    # payload bit from K syndrome shots by a per-check log-likelihood ratio
    r0 = rates_for(base, 0, alpha); r1 = rates_for(base, 1, alpha)
    p0 = fire_prob(r0, obs); p1 = fire_prob(r1, obs)
    if not trained:
        r_uni = np.full(N, base.mean()); p0 = p1 = fire_prob(r_uni, obs)
    w1 = np.log(p1 / p0); w0 = np.log((1 - p1) / (1 - p0))
    hits = 0
    for _ in range(trials):
        b = int(rng.integers(2))
        e = (rng.random((K, N)) < (r1 if b else r0)).astype(np.int8)
        s = (e @ H[obs].T) % 2
        llr = float((s * w1 + (1 - s) * w0).sum())
        b_hat = 1 if llr > 0 else (0 if llr < 0 else int(rng.integers(2)))
        hits += int(b_hat == b)
    return hits / trials


def bits_leaked(acc):
    # mutual information of a binary symmetric channel with error p = 1-acc
    p = min(max(1 - acc, 1e-9), 1 - 1e-9)
    return 1.0 + p * math.log2(p) + (1 - p) * math.log2(1 - p)


def main():
    rng = np.random.default_rng(20260816)
    base = device_rates()
    allrows = np.arange(M)
    trials = 400
    print(f"# covert channel on toric L={L}: one bit carried by data-dependent noise, read from K syndrome shots")
    print(f"# payload footprint Q: logical loop of {len(Q)} qubits\n")

    acc_null = attack_accuracy(base, 0.0, 200, allrows, trials, np.random.default_rng(1))
    acc_open = attack_accuracy(base, 0.6, 200, allrows, trials, np.random.default_rng(1))
    print("# alpha 0 and alpha 0.6:")
    print(f"#   alpha=0.0 (data-independent noise), full access, 200 shots : accuracy {acc_null:.3f}")
    print(f"#   alpha=0.6 (data-dependent noise)   , full access, 200 shots : accuracy {acc_open:.3f}")

    acc_agn = attack_accuracy(base, 0.6, 200, allrows, trials, np.random.default_rng(1), trained=False)
    print(f"#\n# decoder with uniform rates against decoder with the device rates:")
    print(f"#   agnostic (uniform rates), alpha=0.6, 200 shots : accuracy {acc_agn:.3f}")
    print(f"#   trained  (device rates)                        : accuracy {acc_open:.3f}")

    print(f"#\n# shots sweep (alpha=0.5, full access): accuracy vs number of syndrome shots observed")
    print(f"#   {'shots':>6}  {'accuracy':>8}  {'bits':>6}")
    shot_rows = []
    for K in (5, 15, 50, 150, 400):
        a = attack_accuracy(base, 0.5, K, allrows, trials, np.random.default_rng(2))
        shot_rows.append((K, a)); print(f"#   {K:>6}  {a:>8.3f}  {bits_leaked(a):>6.2f}")

    print(f"#\n# access sweep (alpha=0.5, 150 shots): accuracy vs fraction of syndrome observed")
    print(f"#   {'access':>6}  {'accuracy':>8}")
    acc_rows = []
    for f in (0.15, 0.3, 0.5, 0.75, 1.0):
        k = max(1, int(round(f * M)))
        accs = [attack_accuracy(base, 0.5, 150, np.sort(rng.permutation(M)[:k]), trials // 2,
                                np.random.default_rng(300 + i)) for i in range(3)]
        a = float(np.mean(accs)); acc_rows.append((f, a)); print(f"#   {f:>6.2f}  {a:>8.3f}")

    print(f"#\n# strength sweep (150 shots, full access): accuracy vs data-dependence strength alpha")
    print(f"#   {'alpha':>6}  {'accuracy':>8}  {'bits':>6}")
    a_rows = []
    for al in (0.0, 0.1, 0.25, 0.5, 1.0, 2.0):
        a = attack_accuracy(base, al, 150, allrows, trials, np.random.default_rng(4))
        a_rows.append((al, a)); print(f"#   {al:>6.2f}  {a:>8.3f}  {bits_leaked(a):>6.2f}")


    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15.2, 4.6))
        ks = [k for k, _ in shot_rows]
        ax[0].plot(ks, [a for _, a in shot_rows], "o-", color="crimson", lw=2.2)
        ax[0].axhline(0.5, color="0.5", ls=":", lw=1.2, label="0.5")
        ax[0].set_xlabel("syndrome shots observed"); ax[0].set_ylabel("accuracy on the bit")
        ax[0].set_title("Accuracy against shots\n(alpha=0.5, full access)"); ax[0].set_ylim(0.45, 1.02)
        ax[0].legend(fontsize=8)
        fs = [f for f, _ in acc_rows]
        ax[1].plot(fs, [a for _, a in acc_rows], "s-", color="#1f77b4", lw=2.2)
        ax[1].axhline(0.5, color="0.5", ls=":", lw=1.2)
        ax[1].set_xlabel("fraction of syndrome observed (partial access)"); ax[1].set_ylabel("accuracy")
        ax[1].set_title("Accuracy against syndrome access\n(alpha=0.5, 150 shots)"); ax[1].set_ylim(0.45, 1.02)
        als = [a for a, _ in a_rows]
        ax[2].plot(als, [a for _, a in a_rows], "^-", color="seagreen", lw=2.2)
        ax[2].axhline(0.5, color="0.5", ls=":", lw=1.2, label="0.5")
        ax[2].set_xlabel("data-dependence strength alpha"); ax[2].set_ylabel("accuracy")
        ax[2].set_title("Accuracy against alpha\n(150 shots, full access)"); ax[2].set_ylim(0.45, 1.02)
        ax[2].legend(fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "covert_channel.png"), dpi=130)
        print("#\n# wrote covert_channel.png")
    except Exception as ex:
        print(f"# plot skipped, {ex}")


if __name__ == "__main__":
    main()
