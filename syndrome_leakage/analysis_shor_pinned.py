# analysis of the pinned Shor-code run on ibm_fez
import json
import os

from syndrome_leakage import estimate as es
from syndrome_leakage import gate_level as gl
from syndrome_leakage import hardware as hw
from syndrome_leakage.core import tvd

MODEL_SHOTS = 16000


def _signature(tqc):
    idx = lambda q: tqc.find_bit(q).index
    two = [(i.operation.name, tuple(idx(q) for q in i.qubits)) for i in tqc.data if len(i.qubits) == 2]
    meas = sorted((idx(i.qubits[0]), tqc.find_bit(i.clbits[0]).index) for i in tqc.data
                  if i.operation.name == "measure")
    used = sorted({idx(q) for i in tqc.data for q in i.qubits if i.operation.name not in ("barrier", "delay")})
    return two, meas, used


def main():
    root = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
    R = json.load(open(os.path.join(root, "hardware_shor_pinned_ibm_fez.json")))
    meta, counts = R["meta"], R["counts"]
    calfile = json.load(open(os.path.join(root, "hardware_shor_pinned_ibm_fez_calibration.json")))
    cal, cmeta = calfile["calibration"], calfile["meta"]
    circs = gl.load_circuits(os.path.join(root, "hardware_shor_pinned_ibm_fez_circuits.qpy"))
    shots, delays = meta["shots"], meta["delays_s"]
    labels = [tuple(x) for x in meta["labels"]]
    sigs = [_signature(c) for c in circs]

    print(f"job {meta['job_id']} on {meta['backend']}, {shots} shots per circuit, "
          f"{meta['usage']['quantum_seconds']} quantum seconds, plan {meta['plan']} and {meta['pricing_type']}")
    print(f"submitted {meta['submitted_utc']}, calibration read {cmeta['fetched_utc']} ({cmeta['note']})")
    q = cal["qubits"]
    print(f"data qubits {meta['data_qubits']}")
    print(f"  T1 us {[round(q[str(i)]['t1'] * 1e6) for i in meta['data_qubits']]}")
    print(f"  T2 us {[round(q[str(i)]['t2'] * 1e6) for i in meta['data_qubits']]}")
    print(f"ancilla qubits {meta['ancilla_qubits']}, readout error "
          f"{[round(q[str(i)]['readout'], 4) for i in meta['ancilla_qubits']]}")
    print(f"{len({tuple(s[2]) for s in sigs})} set of qubits across the {len(circs)} circuits; "
          f"two-qubit gates and measurement map identical on every circuit: {all(s == sigs[0] for s in sigs)}; "
          f"{len(sigs[0][0])} cz each")
    print()

    d = {lab: es.dist_from_counts(counts[f"{lab[0]}_{lab[1]}"], 64)[0] for lab in labels}
    print("1. The records through estimate.py")
    print()
    rows = [(f"{t * 1e6:.0f} us", es.leak_from_dists(d[(t, 0)], shots, d[(t, 1)], shots, boots=2000, splits=16,
                                                    seed=1)) for t in delays]
    print(es.report(rows, "   six Z generators, 64 outcomes"))
    print()
    rows = [(f"{t * 1e6:.0f} us", es.leak_from_dists(d[(t, 0)], shots, d[(t, 1)], shots, boots=2000, splits=16,
                                                    seed=2, transform=hw.syndrome_weight_dist)) for t in delays]
    print(es.report(rows, "   number of generators that fired, 7 outcomes"))
    print()

    print(f"2. Against the gate-level model built from the calibration read before submission, {MODEL_SHOTS} "
          f"model shots per circuit")
    print()
    sim = gl.run_round(circs, cal, MODEL_SHOTS, gl.spaced_seeds(len(circs), MODEL_SHOTS))
    m = {lab: es.dist_from_counts(c, 64)[0] for lab, c in zip(labels, sim)}
    print(f"   {'delay us':>8} {'measured':>9} {'model':>9} {'measured / model':>17} "
          f"{'trivial 0_L':>12} {'model':>7} {'trivial 1_L':>12} {'model':>7}")
    for t in delays:
        meas, mod = tvd(d[(t, 0)], d[(t, 1)]), tvd(m[(t, 0)], m[(t, 1)])
        ratio = f"{meas / mod:17.3f}" if t > 0 else f"{'':>17}"
        print(f"   {t * 1e6:8.0f} {meas:9.4f} {mod:9.4f} {ratio} {d[(t, 0)][0]:12.4f} {m[(t, 0)][0]:7.4f} "
              f"{d[(t, 1)][0]:12.4f} {m[(t, 1)][0]:7.4f}")
    print()

    print("3. Beside the 2026-09-20 run of the same design, whose circuits the transpiler placed on three layouts")
    print()
    S = json.load(open(os.path.join(root, "hardware_shor_ibm_fez.json")))
    old = {tuple(x): es.dist_from_counts(S["counts"][f"{x[0]}_{x[1]}"], 64)[0] for x in S["meta"]["labels"]}
    print(f"   {'delay us':>8} {'pinned':>9} {'three layouts':>14}")
    for t in delays:
        print(f"   {t * 1e6:8.0f} {tvd(d[(t, 0)], d[(t, 1)]):9.4f} {tvd(old[(t, 0)], old[(t, 1)]):14.4f}")


if __name__ == "__main__":
    main()
