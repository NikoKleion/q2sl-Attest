# stabilizer codes; each constructor returns a Code
import numpy as np
from .core import Code


def _toric_hz_strings(L):
    n = 2 * L * L
    Hh = lambda r, c: (r % L) * L + (c % L)
    Vt = lambda r, c: L * L + (r % L) * L + (c % L)
    out = []
    for r in range(L):
        for c in range(L):
            row = ["I"] * n
            for idx in (Hh(r, c), Hh(r + 1, c), Vt(r, c), Vt(r, c + 1)):
                row[idx] = "Z"
            out.append("".join(row))
    return out


def repetition():
    # 3-qubit bit-flip code
    return Code("3-qubit repetition", 3, ["ZZI", "IZZ"], "ZII", "XXX")


def five_qubit():
    # [[5,1,3]] code
    return Code("5-qubit [[5,1,3]]", 5, ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"], "ZZZZZ", "XXXXX")


def steane():
    # [[7,1,3]] Steane code
    return Code("Steane [[7,1,3]]", 7,
                ["IIIXXXX", "IXXIIXX", "XIXIXIX", "IIIZZZZ", "IZZIIZZ", "ZIZIZIZ"], "ZZZZZZZ", "XXXXXXX")


def code_4_1_2():
    # [[4,1,2]] code
    return Code("[[4,1,2]]", 4, ["XXXX", "ZZII", "IIZZ"], "ZIZI", "XXII")


# Shor [[9,1,3]]: blocks (0,1,2), (3,4,5), (6,7,8).
SHOR_HX = np.array([[1, 1, 1, 1, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 1, 1, 1, 1]], np.uint8)
SHOR_X_CHECKS = ((0, 1, 2, 3, 4, 5), (3, 4, 5, 6, 7, 8))
SHOR_Z_PAIRS = ((0, 1), (1, 2), (3, 4), (4, 5), (6, 7), (7, 8))
SHOR_HZ = np.zeros((6, 9), np.uint8)
for _i, _pair in enumerate(SHOR_Z_PAIRS):
    SHOR_HZ[_i, list(_pair)] = 1


def shor_code():
    # [[9,1,3]] Shor code as generator strings, in the extraction circuit's classical-bit order
    from .css import css_strings
    return css_strings(SHOR_HX, SHOR_HZ, "Shor [[9,1,3]]")


def hamming_7():
    # CSS code with the Hamming [7,4,3] X side, distance 2
    return Code("Hamming [[7,1,2]]", 7,
                ["XIXIXIX", "IXXIIXX", "IIIXXXX", "ZZZZZZZ", "ZIIZZII", "IZIZIZI"], "IIZIZZI", "XXXIXXX")


STANDARD = {
    "repetition": repetition,
    "five_qubit": five_qubit,
    "steane": steane,
    "code_4_1_2": code_4_1_2,
    "hamming_7": hamming_7,
}

# codes held as generator strings only; the command line reports them the way it reports a code file
STRINGS = {
    "shor": shor_code,
}
