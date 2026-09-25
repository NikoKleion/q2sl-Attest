# a transpiled round on a device model built from per-qubit calibration
import json
import math
import os

import numpy as np


def calibration_from_target(target, qubits, pairs=None):
    """Per-qubit T1, T2, sx and x error and duration, readout error; cz error and duration per pair."""
    qp = target.qubit_properties
    cal = {"qubits": {}, "pairs": {}, "dt": target.dt}
    for q in qubits:
        t1 = float(qp[q].t1)
        row = {"t1": t1, "t2": float(qp[q].t2) if qp[q].t2 is not None else t1}
        for g in ("sx", "x"):
            props = target[g].get((q,)) if g in target.operation_names else None
            row[g] = {"error": float(props.error or 0.0), "duration": float(props.duration or 0.0)} if props else None
        m = target["measure"].get((q,))
        row["readout"] = float(m.error or 0.0) if m else 0.0
        cal["qubits"][str(q)] = row
    for name in ("cz", "ecr", "cx"):
        if name not in target.operation_names:
            continue
        for (a, b), props in target[name].items():
            if a in qubits and b in qubits and (pairs is None or (a, b) in pairs or (b, a) in pairs):
                cal["pairs"][f"{a},{b}"] = {"gate": name, "error": float(props.error or 0.0),
                                           "duration": float(props.duration or 0.0)}
    return cal


def save_calibration(cal, path, **meta):
    json.dump({"meta": meta, "calibration": cal}, open(path, "w"), indent=1)


def load_calibration(path):
    return json.load(open(path))["calibration"]


def load_circuits(path, keys_path=None):
    """Transpiled circuits from a QPY file; a dict when keys_path names their labels."""
    from qiskit import qpy
    with open(path, "rb") as fh:
        try:
            circs = list(qpy.load(fh))
        except Exception as exc:
            raise RuntimeError(
                f"{os.path.basename(path)} was written by a different Qiskit version and this one "
                f"cannot read it: {exc}. The version that wrote it is recorded as qiskit_version in "
                f"results/hardware_shor_ibm_fez_job.json. Install that version, or refetch the "
                f"circuits from the job.") from exc
    if keys_path is None:
        return circs
    with open(keys_path, encoding="utf-8") as fh:
        return dict(zip(json.load(fh), circs))


def compress(tqc):
    """A transpiled circuit rebuilt on the physical qubits it uses; returns (circuit, used qubits)."""
    from qiskit import QuantumCircuit
    used = sorted({tqc.find_bit(q).index for inst in tqc.data for q in inst.qubits
                   if inst.operation.name not in ("barrier", "delay")})
    pos = {q: i for i, q in enumerate(used)}
    out = QuantumCircuit(len(used), tqc.num_clbits)
    for inst in tqc.data:
        if inst.operation.name == "barrier":
            continue
        if any(tqc.find_bit(q).index not in pos for q in inst.qubits):
            continue
        out.append(inst.operation, [pos[tqc.find_bit(q).index] for q in inst.qubits],
                   [tqc.find_bit(c).index for c in inst.clbits])
    return out, used


def _relax(t1, t2, dur):
    from qiskit_aer.noise import thermal_relaxation_error
    return thermal_relaxation_error(t1, min(t2, 2 * t1), max(dur, 0.0))


def with_idle_noise(circ, used, cal):
    """Every delay replaced by a relaxation channel of its own duration on that qubit."""
    out = circ.copy_empty_like()
    dt = cal.get("dt") or 1.0
    for inst in circ.data:
        if inst.operation.name == "delay":
            q = circ.find_bit(inst.qubits[0]).index
            row = cal["qubits"][str(used[q])]
            dur = inst.operation.duration * (dt if inst.operation.unit == "dt" else 1.0)
            if inst.operation.unit == "s":
                dur = inst.operation.duration
            out.append(_relax(row["t1"], row["t2"], dur).to_instruction(), [inst.qubits[0]])
        else:
            out.append(inst.operation, inst.qubits, inst.clbits)
    return out


def noise_model(used, cal):
    """Gate and readout noise for the compressed circuit, from the calibration of the qubits it uses."""
    from qiskit_aer.noise import NoiseModel, ReadoutError, depolarizing_error
    nm = NoiseModel(basis_gates=["rz", "sx", "x", "cz", "ecr", "cx", "id", "kraus"])
    pos = {q: i for i, q in enumerate(used)}
    for q in used:
        row = cal["qubits"][str(q)]
        for g in ("sx", "x"):
            if row.get(g):
                err = depolarizing_error(max(row[g]["error"], 0.0), 1)
                nm.add_quantum_error(err.compose(_relax(row["t1"], row["t2"], row[g]["duration"])), g, [pos[q]])
        if row["readout"]:
            e = row["readout"]
            nm.add_readout_error(ReadoutError([[1 - e, e], [e, 1 - e]]), [pos[q]])
    for key, props in cal["pairs"].items():
        a, b = (int(x) for x in key.split(","))
        if a in pos and b in pos:
            ra, rb = cal["qubits"][str(a)], cal["qubits"][str(b)]
            err = depolarizing_error(max(props["error"], 0.0), 2)
            relax = _relax(ra["t1"], ra["t2"], props["duration"]).expand(_relax(rb["t1"], rb["t2"], props["duration"]))
            nm.add_quantum_error(err.compose(relax), props["gate"], [pos[a], pos[b]])
    return nm


def run_round(tqcs, cal, shots, seeds):
    """Sample each transpiled circuit under the device model; returns one counts dict per circuit."""
    from qiskit_aer import AerSimulator
    out = []
    for circ, seed in zip(tqcs, seeds):
        small, used = compress(circ)
        noisy = with_idle_noise(small, used, cal)
        sim = AerSimulator(method="statevector", noise_model=noise_model(used, cal))
        out.append(sim.run(noisy, shots=shots, seed_simulator=seed).result().get_counts())
    return out


def spaced_seeds(n, shots, base=1_000_003):
    """Seeds farther apart than the shot count, so no two runs share Aer's random stream."""
    return [base * (i + 1) + 40 * base * (shots // 4000) for i in range(n)]
