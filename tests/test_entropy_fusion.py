# Tests for entropy_fusion: closed-form exactness, checksum convolution, adaptivity.
import numpy as np

import entropy_fusion as ef
from entropy_fusion.core import fuse, residual_bits, evaluate, shapley, _constrain_checksum, _constrain_checksum_ref
from entropy_fusion.fast import evaluate_fast, shapley_fast
from entropy_fusion.experts import HardwareExpert, ObservedReadExpert
from entropy_fusion.adaptive import evaluate_adaptive
from entropy_fusion.fisher import fisher_cosine


def _phone_experts(pay):
    return [ef.make_expert("field_format"), ef.make_expert("pattern", pay),
            ef.make_expert("hardware", reliability=0.85), ef.make_expert("revealed")]


def test_fast_evaluate_matches_mc():
    pay = ef.make_payload("phone_number"); experts = _phone_experts(pay)
    for r in (0, 2, 5):
        fb, fa = evaluate_fast(pay, experts, reveals=r)
        mb, ma = evaluate(pay, experts, reveals=r, samples=6000)
        assert abs(fb - mb) < 0.05, f"reveals={r}: fast bits {fb} vs MC {mb}"
        assert abs(fa - ma) < 0.03


def test_fast_shapley_matches_mc():
    pay = ef.make_payload("phone_number")
    experts = [ef.make_expert("field_format"), ef.make_expert("pattern", pay),
               ef.make_expert("hardware", reliability=0.85)]
    phi_f = np.asarray(shapley_fast(pay, experts, reveals=2))
    phi_mc = np.asarray(shapley(pay, experts, reveals=2, samples=4000))
    assert np.max(np.abs(phi_f - phi_mc)) < 0.15


def test_fast_returns_none_for_checksum_payload():
    cc = ef.make_payload("credit_card")
    experts = [ef.make_expert("field_format"), ef.make_expert("pattern", cc)]
    assert evaluate_fast(cc, experts, reveals=2) is None
    assert shapley_fast(cc, experts, reveals=2) is None


def test_checksum_convolution_exact():
    cc = ef.make_payload("credit_card"); F, mod = cc.checksum
    rng = np.random.default_rng(0); worst = 0.0
    for _ in range(50):
        post = rng.random((16, 10)); post /= post.sum(axis=1, keepdims=True)
        a = _constrain_checksum(post.copy(), F, mod)
        b = _constrain_checksum_ref(post.copy(), F, mod)
        worst = max(worst, float(np.abs(a - b).max()))
    assert worst < 1e-12, f"vectorized checksum deviates from reference by {worst}"


def test_regression_hardware_channel_is_clean():
    # low-entropy instance seeding inflated P(read correct)
    for payname in ("phone_number", "device_error"):
        pay = ef.make_payload(payname); hw = HardwareExpert(reliability=0.85)
        rng = np.random.default_rng(1); ok = tot = 0
        for _ in range(4000):
            inst = pay.generate(rng); L = hw.likelihood(pay, inst, np.zeros(pay.n, bool))
            ok += int((np.argmax(L, axis=1) == inst).sum()); tot += pay.n
        assert abs(ok / tot - 0.85) < 0.015, f"{payname}: channel P(correct)={ok/tot} not ~0.85"


def test_regression_reveals_always_apply():
    # fuse() must apply reveals even without a RevealExpert
    pay = ef.make_payload("phone_number")
    experts = [ef.make_expert("field_format"), ef.make_expert("pattern", pay)]
    b0, _ = evaluate_fast(pay, experts, reveals=0)
    b5, _ = evaluate_fast(pay, experts, reveals=5)
    assert b5 < b0 - 5.0, "revealing 5 of 10 digits must drop residual entropy substantially"


def test_adaptive_beats_static_under_wrong_prior():
    from entropy_fusion.registry import SpecPayload
    import numpy as np
    q = 1.0 / (1.0 + np.exp(-2.5 * np.random.default_rng(0).standard_normal(16)))
    good = SpecPayload("tok", 16, 2, lambda r: (r.random(16) < q).astype(int)); good.marginals = lambda: np.stack([1 - q, q], 1)
    qb = 1.0 / (1.0 + np.exp(-2.5 * np.random.default_rng(999).standard_normal(16)))
    bad_pay = SpecPayload("tok-wrong", 16, 2, lambda r: (r.random(16) < qb).astype(int)); bad_pay.marginals = lambda: np.stack([1 - qb, qb], 1)
    hw = ef.make_expert("hardware", reliability=0.85)
    prior_bad = ef.make_expert("pattern", bad_pay, temperature=0.5)
    _, a_hw = evaluate_fast(good, [hw, ef.make_expert("revealed")], reveals=0)
    _, a_st = evaluate_fast(good, [prior_bad, hw, ef.make_expert("revealed")], reveals=0)
    _, a_ad, w = evaluate_adaptive(good, prior_bad, hw, (ef.make_expert("revealed"),), reveals=0, samples=200)
    assert a_st < a_hw - 0.2, "static fusion should collapse under a confident-wrong prior"
    assert a_ad > a_st + 0.2, "adaptive should recover well above the collapsed static fusion"
    assert w < 0.3, "trust weight should collapse for a wrong prior"


def test_observed_read_expert_reconstructs():
    pay = ef.make_payload("phone_number")
    rng = np.random.default_rng(1); true = pay.generate(rng)
    observed = true.copy(); flip = rng.random(pay.n) < 0.15
    observed[flip] = (observed[flip] + rng.integers(1, 10, flip.sum())) % 10
    experts = [ef.make_expert("field_format"), ef.make_expert("pattern", pay), ObservedReadExpert(observed, 0.85)]
    post = fuse(pay, true, np.zeros(pay.n, bool), experts)
    acc = float(np.mean(np.argmax(post, axis=1) == true))
    assert acc > 0.7, f"real observed reads should reconstruct most digits, got {acc}"


def test_fisher_cosine_bounds():
    pay = ef.make_payload("phone_number")
    c = fisher_cosine(pay, ef.make_expert("pattern", pay), ef.make_expert("hardware", reliability=0.9))
    assert -1.0 <= c <= 1.0


class _FixedExpert:
    # a likelihood fixed in advance, so the fused posterior can be enumerated by hand
    name = "fixed"

    def __init__(self, L):
        self.L = L

    def likelihood(self, payload, instance, known):
        return self.L


def test_fused_posterior_equals_enumeration():
    # every instance of a five-digit Luhn number, weighted by the product of two experts' likelihoods
    import itertools
    from entropy_fusion.core import Payload
    cc = ef.make_payload("credit_card")
    F, mod = cc.checksum

    class FiveDigits(Payload):
        n, A = 5, 10
        checksum = (F[-5:], mod)

    pay = FiveDigits()
    rng = np.random.default_rng(3)
    experts = [_FixedExpert(rng.random((5, 10))), _FixedExpert(rng.random((5, 10)))]
    L = experts[0].L * experts[1].L
    exact = np.zeros((5, 10))
    for x in itertools.product(range(10), repeat=5):
        if sum(int(F[-5:][i, d]) for i, d in enumerate(x)) % mod == 0:
            w = np.prod([L[i, d] for i, d in enumerate(x)])
            for i, d in enumerate(x):
                exact[i, d] += w
    exact /= exact.sum(axis=1, keepdims=True)
    post = fuse(pay, np.zeros(5, int), np.zeros(5, bool), experts)
    assert np.abs(post - exact).max() < 1e-12


def test_shapley_values_sum_to_the_total_reduction():
    # efficiency: the attributions add up to the drop from no expert to all of them
    pay = ef.make_payload("phone_number")
    experts = [ef.make_expert("field_format"), ef.make_expert("pattern", pay),
               ef.make_expert("hardware", reliability=0.85)]
    phi = np.asarray(shapley_fast(pay, experts, reveals=2))
    blind, _ = evaluate_fast(pay, [], reveals=2)
    full, _ = evaluate_fast(pay, experts, reveals=2)
    assert abs(phi.sum() - (blind - full)) < 1e-9
