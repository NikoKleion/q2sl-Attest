# circuits and analysis for the syndrome leak on a backend
import math

import numpy as np

from .codes import SHOR_HX, SHOR_HZ, SHOR_X_CHECKS, SHOR_Z_PAIRS, shor_code  # noqa: F401
from .core import tvd
from .eavesdrop import chernoff_exponent, ml_attack

DATA = (0, 1, 2)
ANCILLA = (3, 4)


def build_circuits(delays_s, qubits=None):
    # per idle delay, |0_L> and |1_L> of the repetition code, then ZZI and IZZ onto two ancillas
    from qiskit import QuantumCircuit
    data = list(qubits[:3]) if qubits else list(DATA)
    anc = list(qubits[3:5]) if qubits else list(ANCILLA)
    circs, labels = [], []
    for t in delays_s:
        for prep in (0, 1):
            qc = QuantumCircuit(5, 2)
            if prep == 1:
                for i in range(3):
                    qc.x(i)
            qc.barrier()
            if t > 0:
                for i in range(3):
                    qc.delay(t, i, unit="s")
            qc.barrier()
            qc.cx(0, 3); qc.cx(1, 3)          # ZZI
            qc.cx(1, 4); qc.cx(2, 4)          # IZZ
            qc.measure(3, 0); qc.measure(4, 1)
            circs.append(qc)
            labels.append((t, prep))
    return circs, labels


def syndrome_dist_from_counts(counts, shots, n_bits=2):
    # ancilla bits -> a distribution over the 2^n_bits syndromes, bit j for generator j
    from .estimate import dist_from_counts
    return dist_from_counts(counts, 2 ** n_bits)[0]


def marginal_dist(d, keep_bits):
    """Distribution over a subset of the generators, from a full syndrome distribution."""
    out = np.zeros(2 ** len(keep_bits))
    for s, p in enumerate(d):
        out[sum(((s >> b) & 1) << j for j, b in enumerate(keep_bits))] += p
    return out


def syndrome_weight_dist(d):
    """Distribution of the number of generators that fired, from a syndrome distribution."""
    n_bits = int(round(math.log2(len(d))))
    out = np.zeros(n_bits + 1)
    for s, p in enumerate(d):
        out[bin(s).count("1")] += p
    return out


SHOR_Z_BITS = tuple(range(2, 8))


def shor_encoder(qc, data):
    # input on data[0], the other eight in |0>
    qc.cx(data[0], data[3]); qc.cx(data[0], data[6])
    for q in (0, 3, 6):
        qc.h(data[q])
    for c, t in ((0, 1), (0, 2), (3, 4), (3, 5), (6, 7), (6, 8)):
        qc.cx(data[c], data[t])


# circuit layout derived from QECops by Jithesh Mithra, MIT licensed, doi:10.5281/zenodo.19410365; see NOTICE
def build_shor_circuits(delays_s, prep_errors=None, checks="all"):
    # per idle delay, |0_L> and |1_L> of the Shor code, then one round of its generators
    from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
    x_checks = SHOR_X_CHECKS if checks == "all" else ()
    n_anc = len(x_checks) + len(SHOR_Z_PAIRS)
    circs, labels = [], []
    for t in delays_s:
        for prep in (0, 1):
            data, anc = QuantumRegister(9, "d"), QuantumRegister(n_anc, "a")
            creg = ClassicalRegister(n_anc, "syn")
            qc = QuantumCircuit(data, anc, creg)
            qc.h(data[0])
            if prep == 1:
                qc.z(data[0])
            shor_encoder(qc, data)
            qc.barrier()
            if t > 0:
                for i in range(9):
                    qc.delay(t, data[i], unit="s")
            if prep_errors and (t, prep) in prep_errors:
                q, pauli = prep_errors[(t, prep)]
                getattr(qc, pauli.lower())(data[q])
            qc.barrier()
            for j, group in enumerate(x_checks):
                qc.h(anc[j])
                for q in group:
                    qc.cx(anc[j], data[q])
                qc.h(anc[j])
            for k, (qi, qj) in enumerate(SHOR_Z_PAIRS):
                j = len(x_checks) + k
                qc.cx(data[qi], anc[j]); qc.cx(data[qj], anc[j])
            for j in range(n_anc):
                qc.measure(anc[j], creg[j])
            circs.append(qc)
            labels.append((t, prep))
    return circs, labels


def css_encoder(qc, data, code):
    # |0_L> of a CSS code: the pivot qubits of Hx in |+>, then CX onto the rest of each row
    from .css import _rref, css_matrices
    Hx, _Hz = css_matrices(code.stab_strings)
    R, piv = _rref(Hx)
    for row, p in zip(R, piv):
        qc.h(data[p])
        for j, bit in enumerate(row):
            if bit and j != p:
                qc.cx(data[p], data[j])


def measure_generator(qc, data, anc, creg, g, j, flag=None, fcreg=None, fault=None):
    # one generator onto one ancilla, with an optional flag qubit (Chao and Reichardt 2018, Fig. 2(c))
    support = [(q, ch) for q, ch in enumerate(g) if ch != "I"]
    z_type = set(g) <= {"I", "Z"}
    flagged = flag is not None and len(support) >= 3

    def couple():
        if z_type:
            qc.cx(flag, anc[j])
        else:
            qc.cx(anc[j], flag)

    if flagged and z_type:
        qc.h(flag)
    if not z_type:
        qc.h(anc[j])
    for done, (q, ch) in enumerate(support):
        if z_type:
            qc.cx(data[q], anc[j])
        elif ch == "X":
            qc.cx(anc[j], data[q])
        elif ch == "Z":
            qc.cz(anc[j], data[q])
        else:
            qc.cy(anc[j], data[q])
        if flagged and done in (0, len(support) - 2):
            couple()
        if fault is not None and fault[0] == done + 1:
            getattr(qc, fault[1].lower())(anc[j])
    if not z_type:
        qc.h(anc[j])
    qc.measure(anc[j], creg[j])
    if flagged:
        if z_type:
            qc.h(flag)
        qc.measure(flag, fcreg[j])


def select_checks(code, checks="all"):
    """Indices of the generators a round measures: all of them, the Z-type, the X-type, or a list."""
    if checks == "all":
        return list(range(len(code.stab_strings)))
    if checks in ("z", "x"):
        want = {"I", checks.upper()}
        return [j for j, g in enumerate(code.stab_strings) if set(g) <= want]
    return list(checks)


def build_extraction_circuits(code, delays_s, checks="all", prep_errors=None, flags=False, faults=None,
                              final_data=False):
    # per idle delay, |0_L> and |1_L> of a CSS code, then one round of the selected generators
    from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
    idx = select_checks(code, checks)
    circs, labels = [], []
    for t in delays_s:
        for prep in (0, 1):
            data, anc = QuantumRegister(code.n, "d"), QuantumRegister(len(idx), "a")
            creg = ClassicalRegister(len(idx), "syn")
            regs = [data, anc, creg]
            fq = fcreg = dcreg = None
            if flags:
                fq, fcreg = QuantumRegister(len(idx), "f"), ClassicalRegister(len(idx), "flag")
                regs += [fq, fcreg]
            if final_data:
                dcreg = ClassicalRegister(code.n, "data")
                regs.append(dcreg)
            qc = QuantumCircuit(*regs)
            css_encoder(qc, data, code)
            if prep == 1:
                for q, ch in enumerate(code.xl_str):
                    if ch in "XY":
                        qc.x(data[q])
            qc.barrier()
            if t > 0:
                for q in range(code.n):
                    qc.delay(t, data[q], unit="s")
            if prep_errors and (t, prep) in prep_errors:
                q, pauli = prep_errors[(t, prep)]
                getattr(qc, pauli.lower())(data[q])
            qc.barrier()
            for j, g in enumerate(idx):
                measure_generator(qc, data, anc, creg, code.stab_strings[g], j,
                                  flag=fq[j] if flags else None, fcreg=fcreg,
                                  fault=(faults or {}).get(j))
            if final_data:
                qc.measure(data, dcreg)
            circs.append(qc)
            labels.append((t, prep))
    return circs, labels


def split_registers(counts, widths):
    """Counts keyed by space-separated registers -> one tuple of ints per key, first register added first."""
    out = {}
    for key, n in counts.items():
        fields = str(key).split()
        assert len(fields) == len(widths), f"expected {len(widths)} registers in {key!r}"
        vals = tuple(int(f, 2) for f in reversed(fields))
        out[vals] = out.get(vals, 0) + n
    return out


def analyse(dists, labels, T1=None):
    # dists keyed by (delay, prep); report the leak against idle delay
    out = []
    for t in sorted({t for t, _ in labels}):
        d0, d1 = dists[(t, 0)], dists[(t, 1)]
        leak = tvd(d0, d1)
        C, _s = chernoff_exponent(d0, d1)
        row = {"delay_s": t, "tvd": leak, "chernoff": C,
               "attack_1shot": ml_attack(d0, d1, 1, trials=200000, seed=1),
               "theory_1shot": 0.5 * (1 - leak)}
        if T1:
            row["gamma"] = 1.0 - math.exp(-t / T1)
        out.append(row)
    return out


def report(rows, T1=None, code_name="3-qubit repetition code, one round of Z-stabilizer extraction"):
    L = [f"measured syndrome leakage, {code_name}", ""]
    head = f"  {'delay (us)':>11} {'gamma':>8} {'syndrome TVD':>13} {'Chernoff':>10} {'1-shot attack':>14}"
    L.append(head)
    L.append("  " + "-" * (len(head) - 2))
    for r in rows:
        g = f"{r.get('gamma', float('nan')):8.4f}" if T1 else "       -"
        L.append(f"  {r['delay_s'] * 1e6:11.1f} {g} {r['tvd']:13.4f} {r['chernoff']:10.4f} "
                 f"{r['attack_1shot']:14.4f}")
    L.append("")
    L.append("  a leak that is flat in delay would indicate preparation or readout asymmetry, not damping")
    return "\n".join(L)


def available():
    try:
        import qiskit  # noqa: F401
        return True
    except Exception:
        return False


def gate_duration(backend, name, qubits, default=None):
    # seconds, from the BackendV2 Target first, then BackendV1 properties
    try:
        inst = backend.target[name][tuple(qubits)]
        if inst is not None and inst.duration is not None:
            return float(inst.duration)
    except Exception:
        pass
    try:
        return float(backend.properties().gate_length(name, list(qubits)))
    except Exception:
        return default


def t1t2_from_backend(backend, qubit):
    try:
        qp = backend.target.qubit_properties[qubit]
        return float(qp.t1), float(qp.t2)
    except Exception:
        p = backend.properties()
        return float(p.t1(qubit)), float(p.t2(qubit))


def damping_gamma(T1, elapsed):
    # amplitude-damping strength over `elapsed` seconds at relaxation time T1
    if T1 <= 0.0:
        return 1.0
    return float(1.0 - math.exp(-float(elapsed) / float(T1)))


def channel_from_backend(backend, qubit, gate_time=None):
    # relaxation channel from this qubit's T1 and T2 over one gate or idle time
    from .channels import channel_from_t1t2
    T1, T2 = t1t2_from_backend(backend, qubit)
    if gate_time is None:
        gate_time = gate_duration(backend, "sx", (qubit,), default=1e-7)
    return channel_from_t1t2(T1, T2, gate_time)


def syndrome_leak_from_backend(backend, code="hamming_7", qubits=None, syndrome_time=None, rounds=(1, 10, 100)):
    # damping strength from the worst data-qubit T1 over one extraction cycle
    from . import codes, eavesdrop
    qubits = list(range(backend.num_qubits)) if qubits is None else list(qubits)
    c = codes.STANDARD[code]() if isinstance(code, str) else code
    use = qubits[:c.n]
    t1s = []
    for q in use:
        try:
            t1s.append(t1t2_from_backend(backend, q)[0])
        except Exception:
            pass
    if not t1s:
        return None
    T1 = float(min(t1s))
    if syndrome_time is None:
        two_q = None
        for g in ("cx", "ecr", "cz"):
            two_q = gate_duration(backend, g, (use[0], use[1]))
            if two_q is not None:
                break
        syndrome_time = (two_q if two_q is not None else 1e-6) * 4.0
    gamma = damping_gamma(T1, syndrome_time)
    res = eavesdrop.attack_point(c, gamma=gamma, rounds=rounds)
    name = getattr(backend, "name", "backend")
    res["backend"] = name() if callable(name) else name
    res["T1"] = T1
    res["syndrome_time"] = float(syndrome_time)
    res["qubits"] = use
    return res
