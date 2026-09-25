# Knill-Laflamme branch matrices and syndrome-resolved blocks
import numpy as np

import syndrome_leakage as sl
from syndrome_leakage.analyze import population_leak
from syndrome_leakage.approx_qec import (branch_matrices, branch_spread, fitted_order, population_split,
                                         syndrome_blocks)
from syndrome_leakage.channels import (amplitude_damping, channel_from_t1t2, coherent_diagonal, dephasing,
                                       depolarizing)

STD = sl.codes.STANDARD


def test_no_jump_branch_matches_leung_et_al_eq_39():
    # Leung, Nielsen, Chuang and Yamamoto (1997) eq. 39: eigenvalues (1-g)^2 and (1 + (1-g)^4) / 2
    code = STD["code_4_1_2"]()
    for g in (0.1, 0.3):
        M = branch_matrices(code, amplitude_damping(g))[(0, 0, 0, 0)]
        ev = sorted(np.linalg.eigvalsh(M))
        ref = sorted([(1 - g) ** 2, 0.5 * (1 + (1 - g) ** 4)])
        assert np.allclose(ev, ref, atol=1e-12), (g, ev, ref)


def test_blocks_sum_to_identity():
    for name in ("repetition", "code_4_1_2", "steane"):
        code = STD[name]()
        for K in (amplitude_damping(0.2), coherent_diagonal(0.3), channel_from_t1t2(163.3e-6, 77e-6, 1e-6)):
            total = sum(syndrome_blocks(code, K).values())
            assert np.allclose(total, np.eye(2), atol=1e-12), name


def test_block_diagonals_give_the_population_leak():
    for name in ("repetition", "code_4_1_2", "hamming_7"):
        code = STD[name]()
        K = amplitude_damping(0.2)
        assert abs(population_split(syndrome_blocks(code, K)) - population_leak(code, K)[0]) < 1e-12, name


def test_pauli_blocks_have_equal_diagonals():
    for name in ("repetition", "code_4_1_2", "hamming_7"):
        code = STD[name]()
        for K in (depolarizing(0.1), dephasing(0.1)):
            for B in syndrome_blocks(code, K).values():
                assert abs(B[0, 0] - B[1, 1]) < 1e-12, name


def test_branch_spread_orders():
    # repetition 1 - (1-g)^3 is order 1; [[4,1,2]] is order 2 in Leung et al.; Steane no-jump is order 3
    for name, order in (("repetition", 1), ("code_4_1_2", 2), ("steane", 3)):
        code = STD[name]()
        slope = fitted_order(lambda g: branch_spread(code, amplitude_damping(g))[0])
        assert abs(slope - order) < 0.1, (name, slope)
