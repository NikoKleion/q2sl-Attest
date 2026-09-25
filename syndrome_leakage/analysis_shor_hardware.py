# analysis of the Shor-code run on ibm_fez of 2026-09-20
import json
import math
import os

import numpy as np

from syndrome_leakage import estimate as es
from syndrome_leakage import expectations as ex
from syndrome_leakage import hardware as hw
from syndrome_leakage.channels import amplitude_damping, per_qubit_t1t2
from syndrome_leakage.core import tvd


def main():
    HERE = os.path.dirname(os.path.abspath(__file__))
    R = json.load(open(os.path.join(os.path.dirname(HERE), "results", "hardware_shor_ibm_fez.json")))
    meta, counts = R["meta"], R["counts"]
    SHOTS = meta["shots"]
    delays = meta["delays_s"]
    code = hw.shor_code()
    Z = hw.SHOR_Z_BITS

    print(f"job {meta['job_id']} on {meta['backend']}, {SHOTS} shots per circuit, "
          f"{meta['usage']['quantum_seconds']} quantum seconds, plan {meta['plan']} and {meta['pricing_type']}")
    print(f"data qubits {meta['data_qubits']}")
    print(f"  T1 us {[round(x * 1e6) for x in meta['t1']]}")
    print(f"  T2 us {[round(x * 1e6) for x in meta['t2']]}")
    print(f"ancilla qubits {meta['ancilla_qubits']}, readout error "
          f"{[round(x, 4) for x in meta['ancilla_readout_error']]}")
    print(f"transpiled depth {sorted(set(meta['depth']))}, cz {sorted(set(meta['cz']))} on every circuit")
    print()

    d = {}
    for t in delays:
        for prep in (0, 1):
            d[(t, prep)] = es.dist_from_counts(counts[f"{t}_{prep}"], 64)[0]

    print("1. The records through estimate.py")
    print()
    rows = []
    for t in delays:
        r = es.leak_from_dists(d[(t, 0)], SHOTS, d[(t, 1)], SHOTS, boots=2000, splits=16, seed=1)
        rows.append((f"{t * 1e6:.0f} us", r))
    print(es.report(rows, "   six Z generators, 64 outcomes"))
    print()
    rows_w = []
    for t in delays:
        r = es.leak_from_dists(d[(t, 0)], SHOTS, d[(t, 1)], SHOTS, boots=2000, splits=16, seed=2,
                               transform=hw.syndrome_weight_dist)
        rows_w.append((f"{t * 1e6:.0f} us", r))
    print(es.report(rows_w, "   number of generators that fired, 7 outcomes"))

    print()
    print("2. Measured against models of the same round")
    print("   uniform: one damping channel at the median T1 of the nine data qubits")
    print("   per qubit: each data qubit's own T1 and T2, in the layout the job ran on")
    print("   readout: the per-qubit model, each syndrome bit flipped at its ancilla's readout error")
    print("   exposure: the readout model over the delay plus the time the scheduled circuit itself takes")
    print()
    t1s, t2s = meta["t1"], [min(b, 2 * a) for a, b in zip(meta["t1"], meta["t2"])]
    med = float(np.median(t1s))
    flips = [0.0, 0.0] + list(meta["ancilla_readout_error"])
    CIRCUIT = 3.8e-6


    def z_leak(kraus, readout=None):
        _l, a, b = ex.population_leak(code, kraus, readout=readout)
        return tvd(hw.marginal_dist(a, Z), hw.marginal_dist(b, Z))


    print(f"   {'delay us':>8} {'measured':>9} {'uniform':>9} {'per qubit':>10} {'readout':>9} {'exposure':>9} "
          f"{'measured / exposure':>20}")
    for t in delays:
        m = tvd(d[(t, 0)], d[(t, 1)])
        expo = z_leak(per_qubit_t1t2(t1s, t2s, t + CIRCUIT), flips)
        if t == 0:
            print(f"   {t * 1e6:8.1f} {m:9.4f} {0.0:9.4f} {0.0:10.4f} {0.0:9.4f} {expo:9.4f} {'under the floor':>20}")
            continue
        uni = z_leak(amplitude_damping(1.0 - math.exp(-t / med)))
        per = z_leak(per_qubit_t1t2(t1s, t2s, t))
        ro = z_leak(per_qubit_t1t2(t1s, t2s, t), flips)
        print(f"   {t * 1e6:8.1f} {m:9.4f} {uni:9.4f} {per:10.4f} {ro:9.4f} {expo:9.4f} {m / expo:20.3f}")

    print()
    print("3. Trivial-syndrome probability, measured")
    print()
    print(f"   {'delay us':>8} {'0_L':>8} {'1_L':>8}")
    for t in delays:
        print(f"   {t * 1e6:8.1f} {d[(t, 0)][0]:8.4f} {d[(t, 1)][0]:8.4f}")


if __name__ == "__main__":
    main()
