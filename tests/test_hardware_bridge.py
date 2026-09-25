# tests for the device bridge in syndrome_leakage.hardware, on a stub backend (no qiskit needed)
import math

from syndrome_leakage import hardware


class _QubitProps:
    def __init__(self, t1, t2):
        self.t1, self.t2 = t1, t2


class _Inst:
    def __init__(self, duration):
        self.duration = duration


class _GateMap:
    def __init__(self, duration):
        self.duration = duration

    def __getitem__(self, qubits):
        return _Inst(self.duration)


class _Target:
    def __init__(self, t1s, cx=5e-7, sx=3.5e-8):
        self.qubit_properties = [_QubitProps(t, 0.6 * t) for t in t1s]
        self._gates = {"cx": _GateMap(cx), "sx": _GateMap(sx)}

    def __getitem__(self, name):
        return self._gates[name]


class _Backend:
    name = "stub"

    def __init__(self, t1s):
        self.target = _Target(t1s)
        self.num_qubits = len(t1s)


def test_damping_gamma_limits():
    assert hardware.damping_gamma(100e-6, 0.0) == 0.0
    assert abs(hardware.damping_gamma(100e-6, 100e-6) - (1.0 - math.exp(-1.0))) < 1e-12
    assert hardware.damping_gamma(0.0, 1e-6) == 1.0


def test_t1t2_and_gate_duration_from_target():
    b = _Backend([120e-6, 90e-6, 150e-6])
    T1, T2 = hardware.t1t2_from_backend(b, 1)
    assert T1 == 90e-6 and abs(T2 - 54e-6) < 1e-18
    assert hardware.gate_duration(b, "cx", (0, 1)) == 5e-7
    assert hardware.gate_duration(b, "ecr", (0, 1), default=1e-6) == 1e-6


def test_channel_from_backend_is_trace_preserving():
    b = _Backend([120e-6, 90e-6])
    K = hardware.channel_from_backend(b, 0)
    s = sum(k.conj().T @ k for k in K)
    assert abs(s[0, 0] - 1) < 1e-12 and abs(s[1, 1] - 1) < 1e-12 and abs(s[0, 1]) < 1e-12


def test_syndrome_leak_from_backend_uses_worst_t1():
    b = _Backend([120e-6, 90e-6, 150e-6, 200e-6])
    r = hardware.syndrome_leak_from_backend(b, code="repetition", rounds=(1, 10))
    assert r["backend"] == "stub" and r["qubits"] == [0, 1, 2]
    assert r["T1"] == 90e-6 and r["syndrome_time"] == 4 * 5e-7
    assert abs(r["gamma"] - hardware.damping_gamma(90e-6, 2e-6)) < 1e-15
    assert r["tvd"] > 0 and [row["rounds"] for row in r["rows"]] == [1, 10]
