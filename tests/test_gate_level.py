# the device model: compression, idles as channels, the noise model, and the saved job's files
import json
import os

import numpy as np

from _harness import needs, skip
from syndrome_leakage import gate_level as gl

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _cal(qubits, pairs):
    cal = {"dt": 4e-9, "qubits": {}, "pairs": {}}
    for q in qubits:
        cal["qubits"][str(q)] = {"t1": 100e-6, "t2": 80e-6, "readout": 0.01,
                                 "sx": {"error": 2e-4, "duration": 32e-9}, "x": {"error": 2e-4, "duration": 32e-9}}
    for a, b in pairs:
        cal["pairs"][f"{a},{b}"] = {"gate": "cz", "error": 3e-3, "duration": 68e-9}
    return cal


def test_spaced_seeds_are_farther_apart_than_the_shots():
    for shots in (4000, 16000):
        s = gl.spaced_seeds(8, shots)
        assert min(b - a for a, b in zip(s, s[1:])) > shots
        assert len(set(s)) == 8


def test_compress_keeps_only_gated_qubits_and_their_delays():
    needs("qiskit")
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(6, 1)
    qc.h(2)
    qc.cz(2, 4)
    qc.delay(100, 4, unit="dt")
    qc.delay(100, 0, unit="dt")                 # an idle on a qubit nothing else touches
    qc.measure(4, 0)
    small, used = gl.compress(qc)
    assert used == [2, 4] and small.num_qubits == 2
    assert small.count_ops().get("delay", 0) == 1


def test_idles_become_relaxation_of_their_own_duration():
    needs("qiskit")
    needs("qiskit_aer")
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(1, 1)
    qc.x(0)
    qc.delay(25_000, 0, unit="dt")              # 100 us at dt = 4 ns, one T1
    qc.measure(0, 0)
    small, used = gl.compress(qc)
    noisy = gl.with_idle_noise(small, used, _cal([0], []))
    assert noisy.count_ops().get("delay", 0) == 0
    from qiskit_aer import AerSimulator
    counts = AerSimulator().run(noisy, shots=20000, seed_simulator=7_000_003).result().get_counts()
    p1 = counts.get("1", 0) / 20000
    assert abs(p1 - np.exp(-1.0)) < 0.02, p1


def test_noise_model_covers_the_gates_and_readout():
    needs("qiskit_aer")
    nm = gl.noise_model([2, 4], _cal([2, 4], [(2, 4)]))
    names = set(nm.noise_instructions)
    assert {"sx", "x", "cz", "measure"} <= names, names


def test_saved_job_files_are_consistent():
    job = json.load(open(os.path.join(ROOT, "results", "hardware_shor_ibm_fez_job.json")))
    calfile = json.load(open(os.path.join(ROOT, "results", "ibm_fez_calibration.json")))
    cal = calfile["calibration"]
    used = set(map(int, cal["qubits"]))
    for lay in job["layouts"]:
        assert set(lay["data_qubits"]) | set(lay["ancilla_qubits"]) <= used
    assert calfile["meta"]["job_id"] == job["job_id"]
    for key, props in cal["pairs"].items():
        a, b = map(int, key.split(","))
        assert a in used and b in used and props["gate"] == "cz"
        assert 0 <= props["error"] < 0.1 and props["duration"] > 0
    for q, row in cal["qubits"].items():
        assert row["t1"] > 0 and row["t2"] > 0 and 0 <= row["readout"] < 0.5
    # a real calibration is not tidy: the file read two days after the job carries one qubit at 0.23
    assert max(r["readout"] for r in cal["qubits"].values()) > 0.1


def test_executed_circuits_match_the_job():
    needs("qiskit")
    job = json.load(open(os.path.join(ROOT, "results", "hardware_shor_ibm_fez_job.json")))
    try:
        circs = gl.load_circuits(os.path.join(ROOT, "results", "hardware_shor_ibm_fez_circuits.qpy"))
    except RuntimeError as exc:
        skip(str(exc))
    assert len(circs) == len(job["labels"]) == len(job["layouts"])
    for c, cz, lay in zip(circs, job["cz"], job["layouts"]):
        assert c.count_ops().get("cz", 0) == cz
        _small, used = gl.compress(c)
        assert sorted(used) == sorted(lay["data_qubits"] + lay["ancilla_qubits"])
    first = job["layouts"][0]
    assert first["data_qubits"] == job["data_qubits"] and first["ancilla_qubits"] == job["ancilla_qubits"]
    # the zero-delay pair shares a layout; at every other delay the two states do not
    for (t, _p), lay in zip(job["labels"], job["layouts"]):
        same = sorted(lay["data_qubits"]) == sorted(first["data_qubits"])
        assert same == (t == 0), (t, lay)
