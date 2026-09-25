# pauli_boundary run: syndrome distributions across logical states under Pauli and non-Pauli channels
import syndrome_leakage as sl
from syndrome_leakage.channels import amplitude_damping, dephasing, depolarizing
from syndrome_leakage.core import tvd
from syndrome_leakage.eavesdrop import state_grid

GRID = state_grid(7, 6)


def _max_pairwise_syndrome_tvd(code, kraus):
    ds = [code.syndrome_dist(code.apply(code.logical_state(t, p), kraus)) for t, p in GRID]
    return max(tvd(ds[i], ds[j]) for i in range(len(ds)) for j in range(i + 1, len(ds)))


def test_pauli_channels_carry_no_state_information():
    for name in ("repetition", "code_4_1_2"):
        code = sl.codes.STANDARD[name]()
        for kraus in (depolarizing(0.1), dephasing(0.1)):
            assert _max_pairwise_syndrome_tvd(code, kraus) < 1e-12


def test_non_pauli_channel_does_carry_state_information():
    code = sl.codes.STANDARD["repetition"]()
    assert _max_pairwise_syndrome_tvd(code, amplitude_damping(0.2)) > 1e-3


def test_protected_codes_hold_at_every_noise_strength():
    # protection is not a small-gamma effect: no growth from gamma 0.1 to 0.99
    for name in ("steane", "five_qubit"):
        code = sl.codes.STANDARD[name]()
        for gamma in (0.1, 0.3, 0.5, 0.7, 0.9, 0.99):
            leak = _max_pairwise_syndrome_tvd(code, amplitude_damping(gamma))
            assert leak < 1e-12, f"{name} at gamma {gamma}: {leak}"
