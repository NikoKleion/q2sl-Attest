# account seed: a 9 digit account holds min(29.9, b) bits from a b bit seed
import os
import numpy as np

BITS_PER_DIGIT = np.log2(10)


def account_from_seed(seed, prng="strong"):
    if prng == "strong":
        return np.random.default_rng(seed).integers(0, 10, 9)
    if prng == "lcg":
        x = (seed * 1103515245 + 12345) & 0x7FFFFFFF
        out = []
        for _ in range(9):
            x = (x * 1103515245 + 12345) & 0x7FFFFFFF
            out.append((x >> 16) % 10)
        return np.array(out)
    if prng == "structured":
        d = [seed % 10, (seed // 10) % 10]
        for i in range(2, 9):
            d.append((d[-1] + d[-2]) % 10)
        return np.array(d[:9])
    raise ValueError(prng)


def brute_force_residual(b, k, prng="strong", trials=40, seed0=12345):
    # reveal k digits, count surviving seeds; residual bits = log2(survivors), mean over trials
    N = 1 << b
    accts = np.stack([account_from_seed(s, prng) for s in range(N)])
    rng = np.random.default_rng(seed0)
    res = []
    for _ in range(trials):
        s_star = rng.integers(N)
        target = accts[s_star]
        pos = rng.permutation(9)[:k]
        mask = np.all(accts[:, pos] == target[pos], axis=1)
        res.append(np.log2(max(1, int(mask.sum()))))
    return float(np.mean(res))


def analytic_residual(b, k):
    return max(0.0, min(9 * BITS_PER_DIGIT, b) - BITS_PER_DIGIT * k)


def trained_prediction_accuracy(steps=2000, B=256):
    # pattern model on uniform against LCG accounts, accuracy at half reveal
    import torch
    from .pattern_net import Net, encode, luhn_check
    BINS = [list(map(int, s)) for s in ["453201", "510510", "401288", "601100"]]
    rng = np.random.default_rng(0)

    def gen(kind, n):
        rows = []
        for _ in range(n):
            b = BINS[rng.integers(len(BINS))]
            if kind == "uniform":
                acct = rng.integers(0, 10, 9)
            elif kind == "lcg":
                acct = account_from_seed(int(rng.integers(1 << 16)), "lcg")
            else:
                acct = account_from_seed(int(rng.integers(100)), "structured")
            p15 = b + [int(a) for a in acct]
            rows.append(np.array(p15 + [luhn_check(p15)], dtype=int))
        return np.stack(rows)

    def train_eval(kind):
        net = Net(); opt = torch.optim.Adam(net.parameters(), 1e-3); lf = torch.nn.CrossEntropyLoss()
        for t in range(steps):
            batch = gen(kind, B)
            known = rng.random((B, 16)) < rng.uniform(0.1, 0.9)
            x, digits, kn = encode(batch, known)
            unk = ~kn
            if unk.sum() == 0:
                continue
            loss = lf(net(x)[unk], digits[unk])
            opt.zero_grad(); loss.backward(); opt.step()
        test = gen(kind, 400); acc = 0.0; cnt = 0
        for row in test:
            known = rng.random(16) < 0.5
            x, _, _ = encode(row[None, :], known[None, :])
            with torch.no_grad():
                p = torch.softmax(net(x)[0], dim=1)
            pr = np.array(p.tolist())
            for i in range(6, 15):
                if not known[i]:
                    acc += int(np.argmax(pr[i]) == row[i]); cnt += 1
        return acc / max(cnt, 1)

    return {"uniform": train_eval("uniform"), "lcg": train_eval("lcg"), "structured": train_eval("structured")}


def main():
    print("# residual account bits against seed entropy")
    b_val = 16
    print(f"# (1) brute-force validation at seed space b={b_val} bits (strong PRNG), residual account bits:")
    print(f"#   {'revealed k':>10}  {'brute-force':>12}  {'analytic b-3.32k':>18}")
    bf_pts = []
    for k in (0, 1, 2, 3, 4):
        bf = brute_force_residual(b_val, k, "strong")
        an = analytic_residual(b_val, k)
        bf_pts.append((k, bf))
        print(f"#   {k:>10}  {bf:>12.2f}  {an:>18.2f}")
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.7))
        bs = np.arange(0, 31)
        for k in (0, 2, 4, 6):
            ax[0].plot(bs, [analytic_residual(b, k) for b in bs], label=f"{k} account digits revealed")
        ax[0].axhline(9 * BITS_PER_DIGIT, ls="--", color="crimson", label="uniform account (29.9 bits)")
        ax[0].set_xlabel("RNG seed entropy b (bits)"); ax[0].set_ylabel("residual account bits")
        ax[0].set_title("Residual account bits against seed entropy\nmin(29.9, b) - 3.32 x revealed")
        ax[0].legend(fontsize=8, loc="upper left")
        ax[1].plot([k for k, _ in bf_pts], [analytic_residual(b_val, k) for k, _ in bf_pts], "-",
                   color="#1f77b4", label="analytic")
        ax[1].plot([k for k, _ in bf_pts], [v for _, v in bf_pts], "o", color="crimson",
                   label=f"brute force (b={b_val})")
        ax[1].set_xlabel("account digits revealed"); ax[1].set_ylabel("residual seed-set bits")
        ax[1].set_title("Brute-force validation: revealed digits\nfilter the seed hypothesis set")
        ax[1].legend(fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(os.getcwd(), "seed_collapse.png"), dpi=130)
        print("# wrote seed_collapse.png in the working directory")
    except Exception as ex:
        print(f"# plot skipped, {ex}")
    print("# (2) trained model account-digit prediction accuracy at half reveal:")
    try:
        r = trained_prediction_accuracy()
        print(f"#   uniform accounts (true RNG)      : {r['uniform']:.2f}")
        print(f"#   weak-LCG high bits               : {r['lcg']:.2f}")
        print(f"#   internally-structured field      : {r['structured']:.2f}")
    except Exception as ex:
        print(f"#   trained-model demo skipped, {ex}")


if __name__ == "__main__":
    main()
