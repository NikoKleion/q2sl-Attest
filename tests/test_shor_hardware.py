# Shor-code hardware circuits on the local Aer simulator
import math

import numpy as np

from _harness import needs
from syndrome_leakage import expectations as ex
from syndrome_leakage import hardware as hw
from syndrome_leakage import protection as pr


def test_generator_order_matches_circuit():
    c = hw.shor_code()
    xs = ["".join("X" if q in g else "I" for q in range(9)) for g in hw.SHOR_X_CHECKS]
    zs = ["".join("Z" if q in g else "I" for q in range(9)) for g in hw.SHOR_Z_PAIRS]
    assert c.stab_strings == xs + zs


def test_marginal_and_weight():
    import math
    from syndrome_leakage.channels import amplitude_damping
    from syndrome_leakage.core import tvd
    code = hw.shor_code()
    leak, d0, d1 = ex.population_leak(code, amplitude_damping(0.14))
    z0, z1 = hw.marginal_dist(d0, hw.SHOR_Z_BITS), hw.marginal_dist(d1, hw.SHOR_Z_BITS)
    x0, x1 = hw.marginal_dist(d0, (0, 1)), hw.marginal_dist(d1, (0, 1))
    assert abs(z0.sum() - 1) < 1e-12 and len(z0) == 64
    assert abs(tvd(z0, z1) - leak) < 1e-12
    assert tvd(x0, x1) < 1e-12
    w = hw.syndrome_weight_dist(z0)
    assert len(w) == 7 and abs(w.sum() - 1) < 1e-12
    assert abs(w[0] - z0[0]) < 1e-12


def test_syndrome_dist_from_counts_bit_order():
    d = hw.syndrome_dist_from_counts({"00000100": 3, "10000000": 1}, 4, n_bits=8)
    assert d[4] == 0.75 and d[128] == 0.25 and d.sum() == 1.0


def test_encoder_prepares_q2sl_logical_states():
    needs("qiskit")
    from qiskit import QuantumCircuit, QuantumRegister
    from qiskit.quantum_info import Statevector
    code = hw.shor_code()
    for prep, t in ((0, 0.0), (1, math.pi)):
        data = QuantumRegister(9, "d")
        qc = QuantumCircuit(data)
        qc.h(0)
        if prep:
            qc.z(0)
        hw.shor_encoder(qc, data)
        v = Statevector(qc).reverse_qargs().data
        assert abs(abs(np.vdot(v, ex.logical_state_vector(code, t, 0.0))) - 1.0) < 1e-12


def test_z_only_round_sees_x_errors_only():
    needs("qiskit")
    needs("qiskit_aer")
    from qiskit import transpile
    from qiskit_aer import AerSimulator
    code = hw.shor_code()
    circs, want = [], []
    for q in range(9):
        for p in "XYZ":
            cs, _ = hw.build_shor_circuits([0.0], prep_errors={(0.0, 0): (q, p)}, checks="z")
            bits = pr.syndrome_of(code, "I" * q + p + "I" * (8 - q))
            circs.append(cs[0])
            want.append(sum(b << j for j, b in enumerate(bits[2:])))
    assert circs[0].num_qubits == 15 and circs[0].num_clbits == 6
    sim = AerSimulator(method="statevector")
    res = sim.run(transpile(circs, sim), shots=20, seed_simulator=3).result()
    for i, s in enumerate(want):
        assert hw.syndrome_dist_from_counts(res.get_counts(i), 20, n_bits=6)[s] == 1.0, (i, s)
    assert want[2] == 0 and want[0] > 0


def test_syndromes_on_aer():
    needs("qiskit")
    needs("qiskit_aer")
    from qiskit import transpile
    from qiskit_aer import AerSimulator
    code = hw.shor_code()
    cases = [None] + [(q, p) for q in range(9) for p in "XYZ"]
    circs, want = [], []
    for e in cases:
        cs, _ = hw.build_shor_circuits([0.0], prep_errors=None if e is None else {(0.0, 0): e, (0.0, 1): e})
        s = 0 if e is None else sum(b << j for j, b in enumerate(pr.syndrome_of(code, "I" * e[0] + e[1] + "I" * (8 - e[0]))))
        circs += cs
        want += [s, s]
    sim = AerSimulator(method="statevector")
    res = sim.run(transpile(circs, sim), shots=20, seed_simulator=1).result()
    for i, s in enumerate(want):
        assert hw.syndrome_dist_from_counts(res.get_counts(i), 20, n_bits=8)[s] == 1.0, (i, s)


def test_measured_shor_records():
    # the saved ibm_fez job: no offset at zero delay, a leak at every other, and the free plan recorded
    import json
    import os
    from syndrome_leakage import estimate as es
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    R = json.load(open(os.path.join(root, "results", "hardware_shor_ibm_fez.json")))
    meta, counts = R["meta"], R["counts"]
    assert meta["plan"] == "open" and meta["pricing_type"] == "free"
    assert meta["usage"]["quantum_seconds"] <= meta["cap_seconds"]
    assert len(meta["data_qubits"]) == 9 and len(meta["ancilla_qubits"]) == 6 and len(meta["t1"]) == 9
    assert len(set(meta["cz"])) == 1                       # both logical states run the same gates
    shots = meta["shots"]
    for t in meta["delays_s"]:
        d0 = es.dist_from_counts(counts[f"{t}_0"], 64)
        d1 = es.dist_from_counts(counts[f"{t}_1"], 64)
        assert d0[1] == shots and d1[1] == shots
        r = es.leak_from_dists(d0[0], shots, d1[0], shots, boots=400, splits=4, seed=1, null_reps=8)
        if t == 0:
            assert r["p_value"] > 0.05 and r["distance_corrected"] < 0.01, r
        else:
            assert r["p_value"] < 0.01 and r["distance_corrected"] > 0.05, (t, r)


def test_pinned_shor_records():
    # the one-layout ibm_fez job: every circuit on the same qubits and gates, the records equal at zero delay
    # and different at every other, and the free plan recorded
    import json
    import os
    from syndrome_leakage import estimate as es
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    R = json.load(open(os.path.join(root, "results", "hardware_shor_pinned_ibm_fez.json")))
    meta, counts = R["meta"], R["counts"]
    assert meta["plan"] == "open" and meta["pricing_type"] == "free"
    assert meta["usage"]["quantum_seconds"] <= meta["cap_s"]
    assert len(meta["data_qubits"]) == 9 and len(meta["ancilla_qubits"]) == 6
    shots = meta["shots"]
    for t in meta["delays_s"]:
        d0 = es.dist_from_counts(counts[f"{t}_0"], 64)
        d1 = es.dist_from_counts(counts[f"{t}_1"], 64)
        assert d0[1] == shots and d1[1] == shots
        r = es.leak_from_dists(d0[0], shots, d1[0], shots, boots=400, splits=4, seed=1, null_reps=8)
        if t == 0:
            assert r["p_value"] > 0.05 and r["distance_corrected"] < 0.01, r
        else:
            assert r["p_value"] < 0.01 and r["distance_corrected"] > 0.01, (t, r)


def test_pinned_circuits_share_one_layout():
    needs("qiskit")
    import os
    from syndrome_leakage import gate_level as gl
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    circs = gl.load_circuits(os.path.join(root, "results", "hardware_shor_pinned_ibm_fez_circuits.qpy"))

    def signature(c):
        idx = lambda q: c.find_bit(q).index
        two = [tuple(idx(q) for q in i.qubits) for i in c.data if len(i.qubits) == 2]
        meas = sorted((idx(i.qubits[0]), c.find_bit(i.clbits[0]).index) for i in c.data if i.operation.name == "measure")
        return two, meas

    assert len(circs) == 8 and all(signature(c) == signature(circs[0]) for c in circs)
