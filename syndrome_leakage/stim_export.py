# the extraction round of hardware.build_extraction_circuits as Stim circuit text
import math

from .css import _rref, css_matrices
from .hardware import select_checks


def twirled_damping(gamma):
    """(px, py, pz) of the Pauli twirl of amplitude damping at strength gamma."""
    s = math.sqrt(1.0 - gamma)
    return gamma / 4.0, gamma / 4.0, (1.0 - s) ** 2 / 4.0


def stim_text(code, prep=0, checks="all", gamma=None, readout=None, flags=False):
    """Stim circuit text: encode, idle under twirled damping, measure the checks, read the data."""
    idx = select_checks(code, checks)
    n, m = code.n, len(idx)
    anc = list(range(n, n + m))
    flag = list(range(n + m, n + 2 * m)) if flags else []
    L = [f"# {code.name}, logical {'1' if prep else '0'}, {m} generators"]
    L.append("R " + " ".join(str(q) for q in range(n + m + len(flag))))
    Hx, _Hz = css_matrices(code.stab_strings)
    R, piv = _rref(Hx)
    for row, p in zip(R, piv):
        L.append(f"H {p}")
        pairs = [f"{p} {j}" for j, bit in enumerate(row) if bit and j != p]
        if pairs:
            L.append("CX " + " ".join(pairs))
    if prep:
        L.append("X " + " ".join(str(q) for q, ch in enumerate(code.xl_str) if ch in "XY"))
    L.append("TICK")
    if gamma:
        px, py, pz = twirled_damping(gamma)
        L.append(f"PAULI_CHANNEL_1({px:.10g}, {py:.10g}, {pz:.10g}) " + " ".join(str(q) for q in range(n)))
        L.append("TICK")
    measured = []
    for j, g in enumerate(idx):
        string = code.stab_strings[g]
        support = [(q, ch) for q, ch in enumerate(string) if ch != "I"]
        z_type = set(string) <= {"I", "Z"}
        flagged = flags and len(support) >= 3
        if flagged and z_type:
            L.append(f"H {flag[j]}")
        if not z_type:
            L.append(f"H {anc[j]}")
        for done, (q, ch) in enumerate(support):
            if z_type:
                L.append(f"CX {q} {anc[j]}")
            else:
                L.append(f"C{ch} {anc[j]} {q}")
            if flagged and done in (0, len(support) - 2):
                L.append(f"CX {flag[j]} {anc[j]}" if z_type else f"CX {anc[j]} {flag[j]}")
        if not z_type:
            L.append(f"H {anc[j]}")
        if flagged and z_type:
            L.append(f"H {flag[j]}")
        measured.append((anc[j], flag[j] if flagged else None, z_type, [q for q, _ch in support]))
    arg = f"({readout:.10g})" if readout else ""
    order = [a for a, _f, _z, _s in measured] + [f for _a, f, _z, _s in measured if f is not None]
    L.append(f"M{arg} " + " ".join(str(q) for q in order))
    for k in range(len(order)):
        L.append(f"DETECTOR rec[{k - len(order)}]")
    L.append(f"M{arg} " + " ".join(str(q) for q in range(n)))
    total = len(order) + n
    pos = {q: i for i, q in enumerate(order)}
    for a, _f, z_type, support in measured:
        if z_type:
            recs = [pos[a] - total] + [len(order) + q - total for q in support]
            L.append("DETECTOR " + " ".join(f"rec[{r}]" for r in recs))
    L.append("OBSERVABLE_INCLUDE(0) " + " ".join(
        f"rec[{len(order) + q - total}]" for q, ch in enumerate(code.zl_str) if ch == "Z"))
    return "\n".join(L) + "\n"


def parse(text):
    """(name, arguments, targets) per instruction line, for checking the text without Stim installed."""
    out = []
    for line in text.splitlines():
        line = line.split("#")[0].strip()
        if not line:
            continue
        head, _, rest = line.partition(" ")
        name, args = head, []
        if "(" in line:
            name = line[:line.index("(")]
            args = [float(x) for x in line[line.index("(") + 1:line.index(")")].split(",")]
            rest = line[line.index(")") + 1:].strip()
        out.append((name.upper(), args, rest.split()))
    return out
