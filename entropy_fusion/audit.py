# entropy_fusion.audit: residual min-entropy with and without the classical prior
import numpy as np

from .core import evaluate
from .registry import SpecPayload
from .experts import PatternExpert, HardwareExpert, RevealExpert


def structured_bits(n, c, seed=0):
    # per-bit bias q_i = sigmoid(c * z_i); c=0 gives all 0.5
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(n)
    q = 1.0 / (1.0 + np.exp(-c * z))
    gen = lambda r: (r.random(n) < q).astype(int)
    return SpecPayload(f"structured-bits(n={n},c={c})", n, 2, gen), q


def hardware_reliability():
    # read reliability p from the covert channel, or a measured value
    try:
        from suite import covert_channel as cc
        base = cc.device_rates()
        acc = cc.attack_accuracy(base, alpha=0.6, K=120, obs=np.arange(cc.M), trials=300, rng=np.random.default_rng(1))
        return float(acc), "covert channel, alpha=0.6, 120 shots"
    except Exception as e:
        return 0.87, f"fallback (covert channel unavailable: {type(e).__name__})"


def rem(payload, experts, reveals, base_bits):
    b, a = evaluate(payload, experts, reveals=reveals, samples=400)
    return base_bits - b, a, b


def main():
    n = 16
    p, psrc = hardware_reliability()
    print("# residual bits with hardware, pattern, and both, by pattern concentration")
    print(f"# hardware read reliability p = {p:.3f}  ({psrc})\n")

    for c in (0.0, 1.0, 2.5):
        pay, q = structured_bits(n, c)
        pat = PatternExpert(pay, train_samples=40000)
        hw = HardwareExpert(reliability=p, read_fraction=1.0)
        rev = RevealExpert()
        print(f"# pattern concentration c={c} (mean |q-0.5| = {np.mean(np.abs(q-0.5)):.2f})")
        print(f"#   {'reveals':>7}  {'rem_hw':>7}  {'rem_pat':>7}  {'rem_both':>8}  {'complementarity':>15}  {'marginal_hw|pat':>15}")
        for r in (0, 2, 4, 8):
            base, _ = evaluate(pay, [rev], reveals=r, samples=400)
            rh, ah, _ = rem(pay, [hw, rev], r, base)
            rp, ap, _ = rem(pay, [pat, rev], r, base)
            rb, ab, _ = rem(pay, [pat, hw, rev], r, base)
            comp = rb - (rp + rh)
            marg = rb - rp
            print(f"#   {r:>7}  {rh:>7.2f}  {rp:>7.2f}  {rb:>8.2f}  {comp:>+15.2f}  {marg:>15.2f}")
        print()

    print("# accuracy with hardware alone and with a pattern trained on a different structure")
    pay, q = structured_bits(n, 2.5, seed=0)
    wrong = structured_bits(n, 2.5, seed=999)[0]
    pat_wrong = PatternExpert(wrong, train_samples=40000, temperature=0.5)
    hw = HardwareExpert(reliability=p); rev = RevealExpert()
    print(f"#   {'reveals':>7}  {'acc hardware':>12}  {'acc hw+wrongpattern':>20}  {'accuracy change':>15}")
    for r in (0, 2, 4, 8):
        _, ah = evaluate(pay, [hw, rev], reveals=r, samples=400)
        _, aw = evaluate(pay, [pat_wrong, hw, rev], reveals=r, samples=400)
        print(f"#   {r:>7}  {ah:>12.3f}  {aw:>20.3f}  {aw-ah:>+15.3f}")



if __name__ == "__main__":
    main()
