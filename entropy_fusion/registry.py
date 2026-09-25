# entropy_fusion.registry: payload and expert registries
import numpy as np
from .core import Payload, Expert

PAYLOADS = {}
EXPERTS = {}


def register_payload(name):
    def deco(cls):
        PAYLOADS[name] = cls
        cls.registry_name = name
        return cls
    return deco


def register_expert(name):
    def deco(cls):
        EXPERTS[name] = cls
        cls.registry_name = name
        return cls
    return deco


class SpecPayload(Payload):
    # payload from a spec
    def __init__(self, name, n, A, generate, field_alphabets=None, checksum=None):
        self.name = name
        self.n = n
        self.A = A
        self._gen = generate
        self._fa = field_alphabets
        self.checksum = checksum

    def generate(self, rng):
        return np.asarray(self._gen(rng), int)

    @property
    def field_alphabets(self):
        return self._fa if self._fa is not None else [np.arange(self.A) for _ in range(self.n)]


def spec_payload(name, n, A, generate, field_alphabets=None, checksum=None, register=True):
    inst_factory = lambda: SpecPayload(name, n, A, generate, field_alphabets, checksum)
    if register:
        PAYLOADS[name] = inst_factory
    return inst_factory


def make_payload(name, **kw):
    if name not in PAYLOADS:
        raise KeyError(f"no payload '{name}'. registered: {sorted(PAYLOADS)}")
    return PAYLOADS[name](**kw) if kw else PAYLOADS[name]()


def make_expert(name, *a, **kw):
    if name not in EXPERTS:
        raise KeyError(f"no expert '{name}'. registered: {sorted(EXPERTS)}")
    return EXPERTS[name](*a, **kw)


def catalog():
    return {"payloads": sorted(PAYLOADS), "experts": sorted(EXPERTS)}
