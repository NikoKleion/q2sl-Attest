# reconstruction scenario: pattern, device and seed likelihoods fused on one posterior
import os
import numpy as np

try:
    import torch
    _HAVE_TORCH = True
except Exception:
    _HAVE_TORCH = False
try:
    from ldpc import bposd_decoder
    _HAVE_LDPC = True
except Exception:
    try:
        from ldpc.bposd_decoder import bposd_decoder
        _HAVE_LDPC = True
    except Exception:
        _HAVE_LDPC = False

from .pattern_net import Net, luhn_check, encode
from .device import Device
from .seed import account_from_seed
from .magic import sre_m2, _kron, _H, _T, _I

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(2026)
BINS = [list(map(int, s)) for s in ["453201", "510510", "401288", "601100", "371449", "353011"]]

DEV = Device(); H = DEV.Hz; M, N = H.shape
_Str, _Etr = DEV.simulate(40000, "x", np.random.default_rng(7))
R_HAT = np.clip(_Etr.mean(axis=0), 1e-4, 0.45)
READ_FIDELITY = 0.9
_ACCT_TABLE = {}
_RELIABILITY = {}
_NET = None


def _net():
    global _NET
    if _NET is None:
        if not _HAVE_TORCH:
            raise RuntimeError("reconstruction scenario needs torch (pip install torch)")
        _NET = Net(); _NET.load_state_dict(torch.load(os.path.join(HERE, "models", "pattern_model.pt"))); _NET.eval()
    return _NET


def acct_table(b):
    if b not in _ACCT_TABLE:
        _ACCT_TABLE[b] = np.stack([account_from_seed(s, "strong") for s in range(1 << b)]).astype(np.int8)
    return _ACCT_TABLE[b]


def measure_reliability(access, samples=150):
    # exact-recovery rate of the device-trained decoder at this syndrome access
    if access not in _RELIABILITY:
        k = max(1, int(round(access * M))); obs = np.sort(rng.permutation(M)[:k])
        dec = bposd_decoder(H[obs], channel_probs=R_HAT.tolist(), max_iter=48,
                            bp_method="ms", osd_method="osd_cs", osd_order=7)
        ok = 0
        for _ in range(samples):
            e = (rng.random(N) < DEV.rx).astype(int); s = (H @ e) % 2
            ok += int(np.array_equal(dec.decode(s[obs]), e))
        _RELIABILITY[access] = ok / samples
    return _RELIABILITY[access]


def carrier_trust(t_gates, n=4, kappa=1.0):
    # gamma in (0, 1] from the stabilizer Renyi entropy; 1 at zero magic
    plus = (_H @ np.array([1, 0], complex))
    psi = _kron([plus] * n).reshape(-1)
    for q in range(min(t_gates, n)):
        op = _kron([_T if i == q else _I for i in range(n)])
        psi = op @ psi
    m2 = sre_m2(psi, n) + max(0, t_gates - n) * 0.4150
    return 1.0 / (1.0 + kappa * m2), m2


TRUST_STABILIZER, _M2_STAB = carrier_trust(0)
TRUST_MAGIC, _M2_MAGIC = carrier_trust(4)


def gen_card(b_seed=None):
    # issuer prefix + account (uniform, or seed-derived if b_seed is set) + Luhn digit
    b = BINS[rng.integers(len(BINS))]
    if b_seed is None:
        acct = rng.integers(0, 10, 9)
    else:
        acct = account_from_seed(int(rng.integers(1 << b_seed)), "strong")
    p15 = b + [int(a) for a in acct]
    return np.array(p15 + [luhn_check(p15)], dtype=int)


# each returns a (16, 10) per-position likelihood
def pattern_expert(card, known):
    x, _, _ = encode(card[None, :], known[None, :])
    with torch.no_grad():
        p = torch.softmax(_net()(x)[0], dim=1)
    return np.array(p.tolist())


def emitter_expert(card, access, twirl_known, magic_trust=1.0):
    # per-digit decoder read at the measured reliability, scaled by gamma, zero if the twirl seed is unknown
    rel = (measure_reliability(access) if twirl_known else 0.0) * magic_trust
    L = np.full((16, 10), 0.1)
    for i in range(16):
        if rng.random() < rel:
            L[i] = (1 - READ_FIDELITY) / 9.0
            L[i, card[i]] = READ_FIDELITY
    return L


def seed_expert(card, known, b_seed):
    # seeds whose account matches the revealed digits; likelihood is the survivor histogram
    tab = acct_table(b_seed)
    acct_pos = np.arange(6, 15)
    rev = [i - 6 for i in acct_pos if known[i]]
    if rev:
        vals = np.array([card[6 + j] for j in rev])
        keep = np.all(tab[:, rev] == vals, axis=1)
        surv = tab[keep] if keep.any() else tab
    else:
        surv = tab
    L = np.full((16, 10), 0.1)
    for j in range(9):
        counts = np.bincount(surv[:, j], minlength=10).astype(float) + 0.5
        L[6 + j] = counts / counts.sum()
    return L


def luhn_clamp(post, card, known):
    # the Luhn digit is fixed once the first 15 are known
    if known[:15].all() or np.argmax(post[:15], axis=1).shape[0] == 15:
        p15 = list(np.argmax(post[:15], axis=1))
        d = luhn_check(p15)
        post[15] = 0.0; post[15, d] = 1.0
    return post


def aggregate(card, known, experts, use_luhn=True):
    post = np.ones((16, 10))
    for L in experts:
        post = post * L
    for i in range(16):
        if known[i]:
            post[i] = 0.0; post[i, card[i]] = 1.0
    post /= post.sum(axis=1, keepdims=True)
    if use_luhn:
        post = luhn_clamp(post, card, known)
        post /= post.sum(axis=1, keepdims=True)
    return post


def accuracy(post, card):
    return float(np.mean(np.argmax(post, axis=1) == card))


# which factors are active
class Posture:
    def __init__(self, pattern=True, emitter=True, seed=True, luhn=True,
                 b_seed=16, access=0.7, twirl_known=True, magic_trust=TRUST_STABILIZER, label=""):
        self.pattern, self.emitter, self.seed, self.luhn = pattern, emitter, seed, luhn
        self.b_seed, self.access, self.twirl_known, self.label = b_seed, access, twirl_known, label
        self.magic_trust = magic_trust

    def experts(self, card, known):
        E = []
        if self.pattern: E.append(pattern_expert(card, known))
        if self.emitter: E.append(emitter_expert(card, self.access, self.twirl_known, self.magic_trust))
        if self.seed and self.b_seed is not None: E.append(seed_expert(card, known, self.b_seed))
        return E


def reveal_curve(posture, reveals, samples=60):
    out = []
    for r in reveals:
        accs = []
        for _ in range(samples):
            card = gen_card(posture.b_seed if posture.seed else None)
            known = np.zeros(16, bool); known[rng.permutation(16)[:r]] = True
            post = aggregate(card, known, posture.experts(card, known), posture.luhn)
            accs.append(accuracy(post, card))
        out.append(np.mean(accs))
    return out


def heatmap_reveals_x_seed(bs, reveals, samples=40):
    acc = np.zeros((len(bs), len(reveals)))
    for bi, b in enumerate(bs):
        acct_table(b)
        for ri, r in enumerate(reveals):
            vals = []
            for _ in range(samples):
                card = gen_card(b)
                known = np.zeros(16, bool); known[rng.permutation(16)[:r]] = True
                post = aggregate(card, known,
                                 [pattern_expert(card, known), emitter_expert(card, 0.7, True),
                                  seed_expert(card, known, b)], True)
                vals.append(accuracy(post, card))
            acc[bi, ri] = np.mean(vals)
    return acc


def emitter_curve(accesses, twirl_known, reveal=6, b=16, samples=50, magic_trust=TRUST_STABILIZER):
    out = []
    for a in accesses:
        vals = []
        for _ in range(samples):
            card = gen_card(b)
            known = np.zeros(16, bool); known[rng.permutation(16)[:reveal]] = True
            post = aggregate(card, known,
                             [pattern_expert(card, known), emitter_expert(card, a, twirl_known, magic_trust),
                              seed_expert(card, known, b)], True)
            vals.append(accuracy(post, card))
        out.append(np.mean(vals))
    return out


def main():
    if not (_HAVE_TORCH and _HAVE_LDPC):
        print("# reconstruction scenario needs torch and ldpc (pip install torch ldpc)"); return
    reveals = [0, 2, 4, 6, 8, 10, 12, 14, 16]
    postures = [
        Posture(pattern=True, emitter=False, seed=False, label="pattern only (train on cards)"),
        Posture(pattern=True, emitter=True, seed=False, access=0.7, twirl_known=True,
                label="+ emitter (device-trained, access 0.7)"),
        Posture(pattern=True, emitter=True, seed=True, b_seed=16, access=0.7, twirl_known=True,
                label="+ account seed b=16 (all factors)"),
        Posture(pattern=True, emitter=True, seed=True, b_seed=16, access=0.7, twirl_known=False,
                label="all factors, twirl seed unknown"),
    ]
    print("# scenario: pattern, device decoder, account seed and twirl factors, against reveals")
    curves = []
    for p in postures:
        c = reveal_curve(p, reveals)
        curves.append(c)
        print(f"#  {p.label:>44}: " + "  ".join(f"{v:.2f}" for v in c))

    bs = [8, 12, 16]
    hm = heatmap_reveals_x_seed(bs, reveals)
    print("# heatmap accuracy (rows = seed bits, cols = reveals):")
    for bi, b in enumerate(bs):
        print(f"#   b={b:>2}: " + "  ".join(f"{v:.2f}" for v in hm[bi]))

    accs = [0.3, 0.5, 0.7, 1.0]
    e_known = emitter_curve(accs, True, magic_trust=TRUST_STABILIZER)
    e_magic = emitter_curve(accs, True, magic_trust=TRUST_MAGIC)
    e_unknown = emitter_curve(accs, False)
    print(f"# magic trust weights: stabilizer carrier gamma={TRUST_STABILIZER:.2f} (M2={_M2_STAB:.2f}), "
          f"magic-rich carrier gamma={TRUST_MAGIC:.2f} (M2={_M2_MAGIC:.2f})")
    print("# emitter access sweep (reveal=6, b=16):")
    print("#   twirl known, stabilizer carrier : " + "  ".join(f"{v:.2f}" for v in e_known))
    print("#   twirl known, magic-rich carrier : " + "  ".join(f"{v:.2f}" for v in e_magic))
    print("#   twirl unknown                   : " + "  ".join(f"{v:.2f}" for v in e_unknown))

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.9))
        cols = ["#888", "#1f77b4", "crimson", "#2ca02c"]
        for c, p, col in zip(curves, postures, cols):
            ax[0].plot(reveals, c, "o-", color=col, lw=2 if "twirl seed unknown" not in p.label else 1.5,
                       ls="--" if "unknown" in p.label else "-", label=p.label)
        ax[0].axhline(0.1, color="black", ls=":", lw=1.3,
                      label="session key + in-flight state (not reconstructed)")
        ax[0].set_xlabel("revealed digits (receiver cribs)"); ax[0].set_ylabel("PAN reconstruction accuracy")
        ax[0].set_title("Reconstruction accuracy by posture\nover revealed digits")
        ax[0].legend(fontsize=7, loc="lower right"); ax[0].set_ylim(0, 1.02)
        im = ax[1].imshow(hm, aspect="auto", origin="lower", cmap="viridis", vmin=0, vmax=1,
                          extent=[reveals[0], reveals[-1], bs[0], bs[-1]])
        ax[1].set_xlabel("revealed digits"); ax[1].set_ylabel("account RNG seed entropy b (bits)")
        ax[1].set_title("Accuracy by seed entropy\nand revealed digits")
        fig.colorbar(im, ax=ax[1], label="reconstruction accuracy")
        ax[2].plot(accs, e_known, "^-", color="crimson", lw=2, label="twirl known, stabilizer carrier")
        ax[2].plot(accs, e_magic, "D-", color="#9467bd", lw=2, label="twirl known, magic-rich carrier")
        ax[2].plot(accs, e_unknown, "s--", color="#888", label="twirl unknown")
        ax[2].set_xlabel("emitter syndrome access"); ax[2].set_ylabel("PAN reconstruction accuracy (reveal=6)")
        ax[2].set_title("Emitter accuracy vs syndrome access:\ntwirl seed and carrier magic")
        ax[2].legend(fontsize=7.5, loc="center right"); ax[2].set_ylim(0, 1.02)
        fig.tight_layout(); fig.savefig(os.path.join(os.getcwd(), "master_scenario.png"), dpi=130)
        print("# wrote master_scenario.png in the working directory")
    except Exception as ex:
        print(f"# plot skipped, {ex}")


if __name__ == "__main__":
    main()
