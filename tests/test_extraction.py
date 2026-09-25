# the generic encoder and extraction round, on codes that carry no hand-written circuit
import math

import numpy as np

from _harness import needs
from syndrome_leakage import expectations as ex
from syndrome_leakage import hardware as hw
from syndrome_leakage import protection as pr
from syndrome_leakage.codes import STANDARD


def _codes():
    return [("repetition", STANDARD["repetition"]()), ("steane", STANDARD["steane"]()),
            ("shor", hw.shor_code()), ("surface", ex.surface_code_3())]


def test_select_checks():
    c = ex.surface_code_3()
    assert len(hw.select_checks(c, "all")) == 8
    z = hw.select_checks(c, "z")
    x = hw.select_checks(c, "x")
    assert len(z) == 4 and len(x) == 4 and set(z) | set(x) == set(range(8))
    assert all(set(c.stab_strings[j]) <= {"I", "Z"} for j in z)
    assert hw.select_checks(c, [1, 3]) == [1, 3]


def test_encoder_prepares_the_logical_states():
    needs("qiskit")
    from qiskit import QuantumCircuit, QuantumRegister
    from qiskit.quantum_info import Statevector
    for name, code in _codes():
        for prep, t in ((0, 0.0), (1, math.pi)):
            data = QuantumRegister(code.n, "d")
            qc = QuantumCircuit(data)
            hw.css_encoder(qc, data, code)
            if prep:
                for q, ch in enumerate(code.xl_str):
                    if ch in "XY":
                        qc.x(data[q])
            v = Statevector(qc).reverse_qargs().data
            assert abs(abs(np.vdot(v, ex.logical_state_vector(code, t, 0.0))) - 1.0) < 1e-12, (name, prep)


def test_round_gives_the_predicted_syndrome():
    needs("qiskit")
    needs("qiskit_aer")
    from qiskit import transpile
    from qiskit_aer import AerSimulator
    sim = AerSimulator(method="statevector")
    for name, code in _codes():
        idx = hw.select_checks(code, "all")
        circs, want = hw.build_extraction_circuits(code, [0.0])[0], [0, 0]
        for q in range(code.n):
            for p in "XYZ":
                bits = pr.syndrome_of(code, "I" * q + p + "I" * (code.n - q - 1))
                s = sum(bits[g] << j for j, g in enumerate(idx))
                cs, _ = hw.build_extraction_circuits(code, [0.0], prep_errors={(0.0, 0): (q, p)})
                circs = circs + [cs[0]]
                want.append(s)
        res = sim.run(transpile(circs, sim), shots=16, seed_simulator=5).result()
        for i, s in enumerate(want):
            assert hw.syndrome_dist_from_counts(res.get_counts(i), 16, n_bits=len(idx))[s] == 1.0, (name, i)


def test_generic_matches_the_hand_written_repetition_circuit():
    needs("qiskit")
    a = hw.build_circuits([0.0])[0][0]
    b = hw.build_extraction_circuits(STANDARD["repetition"](), [0.0])[0][0]
    two = lambda c: sum(1 for i in c.data if i.operation.num_qubits == 2)
    assert (a.num_qubits, two(a)) == (b.num_qubits, two(b))


def test_flag_fires_on_the_fault_that_spreads_and_on_nothing_else():
    # Chao and Reichardt 2018, Fig. 2(c): a no-fault round never flags; a fault after the second data
    # gate, which spreads to two data qubits, always does; after the third or fourth it does not
    needs("qiskit")
    needs("qiskit_aer")
    from qiskit import transpile
    from qiskit_aer import AerSimulator
    sim = AerSimulator(method="statevector")
    for name, code in (("steane", STANDARD["steane"]()), ("surface", ex.surface_code_3())):
        m = len(code.stab_strings)
        j = next(i for i, g in enumerate(code.stab_strings) if sum(c != "I" for c in g) == 4)
        pauli = "Z" if set(code.stab_strings[j]) <= {"I", "Z"} else "X"
        clean = hw.build_extraction_circuits(code, [0.0], flags=True)[0][0]
        res = sim.run(transpile(clean, sim), shots=32, seed_simulator=9_000_011).result().get_counts()
        assert all(flag == 0 and syn == 0 for (syn, flag) in hw.split_registers(res, (m, m))), name
        for done, want in ((2, True), (3, False), (4, False)):
            qc = hw.build_extraction_circuits(code, [0.0], flags=True, faults={j: (done, pauli)})[0][0]
            res = sim.run(transpile(qc, sim), shots=32,
                          seed_simulator=50_000_017 + 104_729 * done).result().get_counts()
            fired = {bool((flag >> j) & 1) for (_syn, flag) in hw.split_registers(res, (m, m))}
            assert fired == {want}, (name, done, fired)
