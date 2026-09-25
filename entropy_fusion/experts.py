# entropy_fusion.experts: expert modules, each giving a per-position likelihood
import hashlib
import numpy as np
from .core import Expert
from .registry import register_expert


@register_expert("revealed")
class RevealExpert(Expert):
    # hard evidence for positions already known
    name = "revealed"

    def likelihood(self, payload, instance, known):
        L = np.ones((payload.n, payload.A))
        for i in range(payload.n):
            if known[i]:
                L[i] = 0.0; L[i, instance[i]] = 1.0
        return L


@register_expert("field_format")
class FieldFormatExpert(Expert):
    # restricts each position to its allowed alphabet from payload.field_alphabets
    name = "field format"

    def likelihood(self, payload, instance, known):
        L = np.zeros((payload.n, payload.A))
        for i, allowed in enumerate(payload.field_alphabets):
            L[i, np.asarray(allowed, int)] = 1.0
        return L


@register_expert("pattern")
class PatternExpert(Expert):
    # learns the per-position marginal from synthetic samples; temperature > 1 softens, < 1 sharpens
    name = "pattern (trained)"

    def __init__(self, train_payload, train_samples=20000, seed=123, temperature=1.0):
        rng = np.random.default_rng(seed)
        n, A = train_payload.n, train_payload.A
        counts = np.full((n, A), 0.5)
        for _ in range(train_samples):
            x = train_payload.generate(rng)
            for i in range(n):
                counts[i, int(x[i])] += 1
        marg = counts / counts.sum(axis=1, keepdims=True)
        if temperature != 1.0:
            marg = marg ** (1.0 / temperature)
            marg = marg / marg.sum(axis=1, keepdims=True)
        self.marginal = marg
        self.trained_on = train_payload.name
        self.samples = train_samples
        self.temperature = temperature

    def likelihood(self, payload, instance, known):
        return self.marginal.copy()


@register_expert("observed")
class ObservedReadExpert(Expert):
    # observed hardware symbols at per-read reliability p
    name = "observed read"

    def __init__(self, observed, reliability, read_positions=None):
        self.observed = np.asarray(observed, int)
        self.p = reliability
        self.read_positions = None if read_positions is None else set(int(i) for i in read_positions)

    def likelihood(self, payload, instance, known):
        n, A = payload.n, payload.A
        p = np.broadcast_to(np.asarray(self.p, float), (n,))
        rd = set(range(n)) if self.read_positions is None else self.read_positions
        L = np.full((n, A), 1.0 / A)
        for i in range(n):
            if i in rd and i < len(self.observed):
                obs = int(self.observed[i]); pi = float(p[i])
                L[i] = (1.0 - pi) / (A - 1)
                L[i, obs] = pi
        return L


@register_expert("hardware")
class HardwareExpert(Expert):
    # each read position is correct with probability reliability, else uniform
    name = "hardware read"

    def __init__(self, reliability, read_fraction=1.0, seed=0):
        self.p = float(reliability)
        self.f = float(read_fraction)
        self.seed = seed
        self._readset = None

    def _reads(self, n):
        if self._readset is None:
            k = max(0, int(round(self.f * n)))
            self._readset = set(np.random.default_rng(self.seed).permutation(n)[:k].tolist())
        return self._readset

    def likelihood(self, payload, instance, known):
        # correct with probability p else uniformly wrong; blake2b hashes keep reads paired
        rd = self._reads(payload.n)
        A = payload.A
        base = np.ascontiguousarray(instance, dtype=np.int64).tobytes() + int(self.seed).to_bytes(4, "little")
        L = np.full((payload.n, A), 1.0 / A)
        for i in range(payload.n):
            if i in rd:
                true = int(instance[i])
                d = hashlib.blake2b(base + int(i).to_bytes(4, "little"), digest_size=16, person=b"hw-read1").digest()
                u1 = int.from_bytes(d[:8], "little") / 2.0 ** 64
                u2 = int.from_bytes(d[8:], "little") / 2.0 ** 64
                if u1 < self.p:
                    obs = true
                else:
                    w = int(u2 * (A - 1))
                    obs = w + (1 if w >= true else 0)
                L[i] = (1.0 - self.p) / (A - 1)
                L[i, obs] = self.p
        return L
