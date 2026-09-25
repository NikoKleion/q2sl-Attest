# entropy_fusion: residual min-entropy of a structured value given several likelihood sources
from .core import Payload, Expert, fuse, residual_bits, evaluate, shapley
from .registry import (register_payload, register_expert, spec_payload, make_payload, make_expert,
                       catalog, PAYLOADS, EXPERTS, SpecPayload)
from . import payloads, experts
from .fast import evaluate_fast, shapley_fast, generating_marginals


def evaluate_auto(payload, experts_, reveals=0, samples=300, seed=0):
    # exact closed form when it applies, else Monte-Carlo evaluate
    r = evaluate_fast(payload, experts_, reveals)
    return r if r is not None else evaluate(payload, experts_, reveals, samples, seed)


__all__ = ["Payload", "Expert", "fuse", "residual_bits", "evaluate", "shapley",
           "evaluate_fast", "shapley_fast", "generating_marginals", "evaluate_auto",
           "register_payload", "register_expert", "spec_payload", "make_payload", "make_expert",
           "catalog", "PAYLOADS", "EXPERTS", "SpecPayload", "payloads", "experts"]
