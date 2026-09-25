# suite.core: module registry and dispatch
MODULES = {}


class Finding:
    def __init__(self, module, dimension, target, headline, metric, confidence, detail="", vulnerable=True):
        self.module = module
        self.dimension = dimension
        self.target = target
        self.headline = headline
        self.metric = metric
        self.confidence = confidence
        self.detail = detail
        self.vulnerable = vulnerable

    def __str__(self):
        c = "high" if self.confidence > 0.8 else "medium" if self.confidence > 0.5 else "low"
        return (f"[{self.dimension}/{self.target}] {self.module}\n"
                f"    {self.headline}\n"
                f"    metric {self.metric}  confidence {c} ({self.confidence:.2f})"
                + (f"\n    {self.detail}" if self.detail else ""))


class Module:
    name = ""
    dimension = ""
    target = ""
    summary = ""

    def run(self, **kw):
        raise NotImplementedError


def register(cls):
    MODULES[cls.name] = cls
    return cls


def catalog():
    return {n: {"dimension": m.dimension, "target": m.target, "summary": m.summary} for n, m in MODULES.items()}


def find(target=None, dimension=None, modules=None, **kw):
    # run every registered module matching the target or dimension
    hits = []
    for n, m in MODULES.items():
        if modules is not None and n not in modules:
            continue
        if target is not None and m.target != target:
            continue
        if dimension is not None and m.dimension != dimension:
            continue
        hits.append(m)
    return [f for f in (m().run(**kw) for m in hits) if f is not None]


def targets():
    return sorted({m.target for m in MODULES.values()})


def dimensions():
    return sorted({m.dimension for m in MODULES.values()})
