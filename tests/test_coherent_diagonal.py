# coherent run: coherent diagonal noise leaks through X-type stabilizers at second order in the angle
import syndrome_leakage as sl
from syndrome_leakage.channels import coherent_diagonal, dephasing, depolarizing
from syndrome_leakage.analyze import population_leak


def test_x_stabilizer_code_leaks_coherent_diagonal():
    c = sl.codes.STANDARD["code_4_1_2"]()
    leak, _, _ = population_leak(c, coherent_diagonal(0.2))
    assert leak > 1e-2


def test_z_only_code_does_not_leak():
    c = sl.codes.STANDARD["repetition"]()
    leak, _, _ = population_leak(c, coherent_diagonal(0.3))
    assert leak < 1e-12


def test_pauli_controls_zero():
    c = sl.codes.STANDARD["code_4_1_2"]()
    for K in (dephasing(0.1), depolarizing(0.1)):
        leak, _, _ = population_leak(c, K)
        assert leak < 1e-12


def test_second_order_in_theta():
    c = sl.codes.STANDARD["code_4_1_2"]()
    a, _, _ = population_leak(c, coherent_diagonal(0.05))
    b, _, _ = population_leak(c, coherent_diagonal(0.10))
    assert abs(b / a - 4.0) < 0.3
