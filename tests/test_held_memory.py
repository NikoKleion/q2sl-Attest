# held_memory run: repeated extraction on the repetition code at device parameters, with a Pauli control
import math

import syndrome_leakage as sl
from syndrome_leakage.channels import channel_from_t1t2, dephasing
from syndrome_leakage.eavesdrop import repeated_extraction

K = channel_from_t1t2(163.3e-6, 77.0e-6, 1e-6)  # ibm_marrakesh medians, 1 us cycle


def test_repetition_held_memory_leaks_at_device_params():
    res = repeated_extraction(sl.codes.STANDARD["repetition"](), kraus_override=K, rounds=200, correct=True)
    total = sum(r["chernoff"] for r in res["rounds"])
    assert total > 1.0
    assert 0.5 * math.exp(-total) < 0.05


def test_pauli_channel_held_memory_leaks_nothing():
    res = repeated_extraction(sl.codes.STANDARD["repetition"](), kraus_override=dephasing(0.05),
                              rounds=20, correct=True)
    assert max(r["tvd"] for r in res["rounds"]) < 1e-12
