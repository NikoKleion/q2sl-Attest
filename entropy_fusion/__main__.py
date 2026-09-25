# entropy_fusion on three synthetic payloads
import sys
import numpy as np

from .core import evaluate, shapley
from . import payloads as P
from .experts import RevealExpert, FieldFormatExpert, PatternExpert


def run(payload, reveals_grid):
    fmt = FieldFormatExpert()
    pat = PatternExpert(payload)
    rev = RevealExpert()
    blind = payload.blind_bits()
    print(f"\n# payload: {payload.name}")
    print(f"#   positions {payload.n}, alphabet {payload.A}, blind (structural) entropy {blind:.1f} bits")
    print(f"#   {'reveals':>7}  {'field only':>11}  {'+ pattern':>10}  {'+ pattern (acc)':>15}")
    for r in reveals_grid:
        b_fmt, _ = evaluate(payload, [fmt, rev], reveals=r)
        b_all, a_all = evaluate(payload, [fmt, pat, rev], reveals=r)
        print(f"#   {r:>7}  {b_fmt:>11.1f}  {b_all:>10.1f}  {a_all:>15.2f}")
    phi = shapley(payload, [fmt, pat], reveals=reveals_grid[len(reveals_grid) // 2])
    print(f"#   Shapley bits removed at {reveals_grid[len(reveals_grid)//2]} reveals: field {phi[0]:+.1f}, pattern {phi[1]:+.1f}")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else None
    print("# entropy_fusion: residual min-entropy (bits) as experts are added")
    print("# synthetic payloads")
    grids = {
        "credit_card": (0, 4, 8, 12),
        "phone_number": (0, 2, 4, 7),
        "device_error": (0, 4, 8, 12),
    }
    items = [(which, P.STANDARD[which]())] if which in P.STANDARD else [(k, P.STANDARD[k]()) for k in P.STANDARD]
    for key, payload in items:
        run(payload, grids[key])


if __name__ == "__main__":
    main()
