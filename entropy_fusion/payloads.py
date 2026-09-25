# entropy_fusion.payloads: synthetic parametric payloads
import numpy as np
from .core import Payload
from .registry import register_payload

BINS = ["453201", "510510", "401288", "601100", "371449", "353011"]
FL_AREA = ["305", "786", "561", "407", "813", "904", "850", "352", "239", "321", "754", "727", "941", "386", "863", "689"]


@register_payload("credit_card")
class CreditCard(Payload):
    name = "credit card (16 digits: issuer prefix + Luhn)"
    n = 16
    A = 10

    def __init__(self):
        self._bins = [list(map(int, b)) for b in BINS]
        self._F = self._luhn_table()
        self.checksum = (self._F, 10)

    def generate(self, rng):
        b = self._bins[rng.integers(len(self._bins))]
        acct = list(rng.integers(0, 10, 9))
        p15 = b + acct
        return np.array(p15 + [self._luhn_check(p15)])

    @property
    def field_alphabets(self):
        fa = [np.array(sorted({b[i] for b in self._bins})) for i in range(6)]
        fa += [np.arange(10) for _ in range(9)]
        fa += [np.arange(10)]
        return fa

    def _luhn_table(self):
        F = np.zeros((16, 10), int)
        for i in range(15):
            for d in range(10):
                v = 2 * d if i % 2 == 0 else d
                if i % 2 == 0 and v > 9:
                    v -= 9
                F[i, d] = v % 10
        F[15] = np.arange(10)
        return F

    def _luhn_check(self, p15):
        s = sum(int(self._F[i, p15[i]]) for i in range(15))
        return (10 - (s % 10)) % 10


@register_payload("phone_number")
class PhoneNumber(Payload):
    name = "Florida phone number (area code + NXX-XXXX)"
    n = 10
    A = 10

    def __init__(self):
        self._areas = [list(map(int, a)) for a in FL_AREA]

    def generate(self, rng):
        a = self._areas[rng.integers(len(self._areas))]
        exch = [int(rng.integers(2, 10))] + list(rng.integers(0, 10, 2))
        line = list(rng.integers(0, 10, 4))
        return np.array(a + exch + line)

    @property
    def field_alphabets(self):
        fa = [np.array(sorted({a[i] for a in self._areas})) for i in range(3)]
        fa += [np.arange(2, 10)]
        fa += [np.arange(10) for _ in range(6)]
        return fa

    def marginals(self):
        # exact per-position generating marginal g_i[d] = P(instance[i] == d)
        g = np.zeros((self.n, self.A))
        for i in range(3):
            for a in self._areas:
                g[i, a[i]] += 1.0 / len(self._areas)
        g[3, 2:10] = 1.0 / 8.0
        g[4:10] = 1.0 / 10.0
        return g


@register_payload("device_error")
class DeviceErrorString(Payload):
    name = "device error string (hardware fingerprint)"
    A = 2

    def __init__(self, length=24, hot_frac=0.3, hot_mult=6.0, base=0.05, seed=7):
        # per-bit error rates of the device
        self.n = length
        rng = np.random.default_rng(seed)
        hot = np.ones(length); hot[rng.random(length) < hot_frac] = hot_mult
        self.rates = np.clip(base * hot, 1e-3, 0.45)

    def generate(self, rng):
        return (rng.random(self.n) < self.rates).astype(int)

    @property
    def field_alphabets(self):
        return [np.arange(2) for _ in range(self.n)]

    def marginals(self):
        # exact per-bit generating marginal from the device fingerprint rates
        return np.stack([1.0 - self.rates, self.rates], axis=1)


STANDARD = {"credit_card": CreditCard, "phone_number": PhoneNumber, "device_error": DeviceErrorString}
