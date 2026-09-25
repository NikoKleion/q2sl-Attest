# device run: leak under channel_from_t1t2 at ibm_marrakesh median T1 and T2, with the T1-removed control
import syndrome_leakage as sl
from syndrome_leakage.channels import channel_from_t1t2
from syndrome_leakage.analyze import population_leak

T1, T2, GT = 163.3e-6, 77.0e-6, 1e-6  # ibm_marrakesh medians, 1 us syndrome cycle


def test_repetition_leaks_at_device_params():
    leak, _, _ = population_leak(sl.codes.STANDARD["repetition"](), channel_from_t1t2(T1, T2, GT))
    assert leak > 1e-3


def test_distance3_codes_protected_at_device_params():
    for name in ("steane", "five_qubit"):
        leak, _, _ = population_leak(sl.codes.STANDARD[name](), channel_from_t1t2(T1, T2, GT))
        assert leak < 1e-12


def test_leak_is_t1_driven():
    # remove amplitude damping by sending T1 to infinity; only phase damping remains and it must not leak
    leak, _, _ = population_leak(sl.codes.STANDARD["repetition"](), channel_from_t1t2(1e9, T2, GT))
    assert leak < 1e-12
