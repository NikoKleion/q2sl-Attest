# worst_pair run: the |0_L>, |1_L> pair over a Bloch grid under amplitude damping and coherent diagonal noise
import syndrome_leakage as sl
from syndrome_leakage.eavesdrop import worst_case_pair, state_grid

GRID = state_grid(9, 8)


def test_axis_is_worst_amplitude_damping():
    for name in ("repetition", "code_4_1_2"):
        w = worst_case_pair(sl.codes.STANDARD[name](), gamma=0.2, channel="amplitude_damping", grid=GRID)
        assert w["axis_is_worst"]


def test_axis_is_worst_coherent_diagonal():
    # [[4,1,2]] leaks under coherent diagonal noise; the worst pair is still the population axis
    w = worst_case_pair(sl.codes.STANDARD["code_4_1_2"](), gamma=0.3, channel="coherent_diagonal", grid=GRID)
    assert w["axis_is_worst"]
    assert w["chernoff"] > 0.1
