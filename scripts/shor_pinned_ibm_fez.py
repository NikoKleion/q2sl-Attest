# the Shor Z-check round on ibm_fez with every circuit pinned to one layout
import json
import math
import os
import sys
import time


REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
OUT = os.path.join(REPO, "build", "pinned_dry_run")
BACKEND = "ibm_fez"
DELAYS = [0.0, 20e-6, 50e-6, 100e-6]
SHOTS = 16000
SEARCH_SEEDS = range(24)
CAP_SECONDS = 120          # hard cap on the job's execution time
PREVIOUS_SECONDS = 41      # the 2026-09-20 job: the same eight circuits at 16000 shots
# the run of 2026-09-25 was job daquif6ekp0c73arbd70: 41 quantum seconds, one layout, transpiler seed 16


def guards(svc, need_seconds):
    insts = svc.instances()
    if len(insts) != 1:
        raise SystemExit(f"expected exactly one instance, found {len(insts)}; stopping")
    inst = insts[0]
    if inst.get("plan") != "open" or inst.get("pricing_type") != "free":
        raise SystemExit(f"instance is not the free open plan ({inst.get('plan')}, {inst.get('pricing_type')}); stopping")
    u = svc.usage()
    if u.get("usage_limit_reached") or u.get("usage_remaining_seconds", 0) < need_seconds:
        raise SystemExit(f"allowance too low: {u.get('usage_remaining_seconds')} s left, {need_seconds} s needed; stopping")
    return inst, u


def layout_of(tqc):
    return list(tqc.layout.initial_index_layout(filter_ancillas=True))


def signature(tqc):
    """The two-qubit gates in order with their physical qubits, and the measurement map."""
    idx = lambda q: tqc.find_bit(q).index
    two = [(i.operation.name, tuple(idx(q) for q in i.qubits)) for i in tqc.data if len(i.qubits) == 2]
    meas = sorted((idx(i.qubits[0]), tqc.find_bit(i.clbits[0]).index) for i in tqc.data if i.operation.name == "measure")
    used = sorted({idx(q) for i in tqc.data for q in i.qubits if i.operation.name not in ("barrier", "delay")})
    return two, meas, used


def score(tqc, target):
    two, meas, used = signature(tqc)
    props = lambda name, q: target[name].get(tuple(q)) or target[name].get(tuple(reversed(q)))
    err = sum((props(name, q).error or 0.0) for name, q in two
              if name in target.operation_names and props(name, q) is not None)
    err += sum((target["measure"][(p,)].error or 0.0) for p, _c in meas)
    return (len(two), err)


def dry():
    from qiskit import qpy
    from qiskit.transpiler import generate_preset_pass_manager
    from qiskit_ibm_runtime import QiskitRuntimeService
    from syndrome_leakage import estimate as es
    from syndrome_leakage import gate_level as gl
    from syndrome_leakage import hardware as hw

    svc = QiskitRuntimeService()
    inst, usage = guards(svc, 2 * PREVIOUS_SECONDS)
    backend = svc.backend(BACKEND)
    status = backend.status()
    target = backend.target
    print(f"instance {inst.get('name')}: plan {inst.get('plan')}, pricing {inst.get('pricing_type')}; "
          f"{usage['usage_remaining_seconds']} s of {usage['usage_limit_seconds']} left this period")
    print(f"{BACKEND}: operational {status.operational}, pending jobs {status.pending_jobs}")

    circs, labels = hw.build_shor_circuits(DELAYS, checks="z")
    ref = circs[0]
    best = None
    for seed in SEARCH_SEEDS:
        pm = generate_preset_pass_manager(optimization_level=3, backend=backend, seed_transpiler=seed,
                                          scheduling_method="alap")
        t = pm.run(ref)
        s = score(t, target)
        if best is None or s < best[0]:
            best = (s, seed, layout_of(t))
    (n2, err), seed, layout = best
    print(f"layout search over {len(SEARCH_SEEDS)} seeds: best seed {seed}, {n2} two-qubit gates, summed "
          f"two-qubit and readout error {err:.4f}")
    print(f"virtual qubits d0..d8, a0..a5 -> physical {layout}")

    pm = generate_preset_pass_manager(optimization_level=3, backend=backend, seed_transpiler=seed,
                                      initial_layout=layout, scheduling_method="alap")
    tqcs = [pm.run(c) for c in circs]
    sigs = [signature(t) for t in tqcs]
    for (t, prep), sig, tq in zip(labels, sigs, tqcs):
        assert layout_of(tq) == layout, f"circuit {(t, prep)} left the pinned layout"
        assert sig == sigs[0], f"circuit {(t, prep)} differs from the first in qubits, two-qubit gates or measurement"
    two, meas, used = sigs[0]
    print(f"all {len(tqcs)} circuits: one layout, the same {len(two)} two-qubit gates on the same pairs, the same "
          f"measurement map, qubits {used}")
    dt = target.dt
    durs = [t.estimate_duration(target, unit='s') for t in tqcs]
    rep = backend.configuration().default_rep_delay
    rep = rep * 1e-6 if rep and rep > 1e-3 else (rep or 250e-6)
    est = sum(SHOTS * (d + rep) for d in durs)
    print(f"scheduled durations {[round(d * 1e6, 1) for d in durs]} us; rep delay {rep * 1e6:.0f} us")
    print(f"estimate: {est:.0f} s of shots plus overhead; the 2026-09-20 job ran the same eight circuits at "
          f"{SHOTS} shots in {PREVIOUS_SECONDS} quantum seconds; cap {CAP_SECONDS} s")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "circuits.qpy"), "wb") as fh:
        qpy.dump(tqcs, fh)
    cal = gl.calibration_from_target(target, used)
    gl.save_calibration(cal, os.path.join(OUT, "calibration.json"), backend=BACKEND,
                        fetched_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                        note="read from the live backend target at the dry run")
    plan = {"backend": BACKEND, "delays_s": DELAYS, "shots": SHOTS, "labels": labels, "layout": layout,
            "data_qubits": layout[:9], "ancilla_qubits": layout[9:], "seed_transpiler": seed,
            "two_qubit_gates": len(two), "used_qubits": used, "durations_s": durs, "estimate_s": est,
            "cap_s": CAP_SECONDS, "made_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    json.dump(plan, open(os.path.join(OUT, "plan.json"), "w"), indent=1)

    print()
    print("rehearsal on the device model built from this calibration, 16000 shots per circuit")
    counts = gl.run_round(tqcs, cal, SHOTS, gl.spaced_seeds(len(tqcs), SHOTS))
    d = {lab: es.dist_from_counts(c, 64)[0] for lab, c in zip(labels, counts)}
    print(f"  {'delay us':>8} {'distance':>9} {'floor':>7} {'p':>6} {'achieved':>9}")
    for t in DELAYS:
        r = es.leak_from_dists(d[(t, 0)], SHOTS, d[(t, 1)], SHOTS, boots=1000, splits=8, seed=1, null_reps=10)
        print(f"  {t * 1e6:8.0f} {r['tvd']:9.4f} {r['null']:7.4f} {r['p_value']:6.3f} {r['distance_corrected']:9.4f}")
    print()
    print(f"saved to {OUT}; nothing was submitted")


def submit():
    from qiskit import qpy
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    import qiskit
    import qiskit_ibm_runtime

    plan = json.load(open(os.path.join(OUT, "plan.json")))
    with open(os.path.join(OUT, "circuits.qpy"), "rb") as fh:
        tqcs = list(qpy.load(fh))
    need = max(2 * PREVIOUS_SECONDS, math.ceil(2 * plan["estimate_s"]))
    svc = QiskitRuntimeService()
    inst, usage = guards(svc, need)
    backend = svc.backend(plan["backend"])
    if not backend.status().operational:
        raise SystemExit(f"{plan['backend']} is not operational; stopping")
    target = backend.target
    for q in plan["used_qubits"]:
        if target.qubit_properties[q] is None or (q,) not in target["measure"]:
            raise SystemExit(f"qubit {q} is no longer usable on {plan['backend']}; stopping")
    sampler = SamplerV2(mode=backend)
    sampler.options.default_shots = plan["shots"]
    sampler.options.max_execution_time = CAP_SECONDS
    sampler.options.dynamical_decoupling.enable = False
    sampler.options.twirling.enable_gates = False
    sampler.options.twirling.enable_measure = False
    import shutil
    from syndrome_leakage import gate_level as gl
    base = os.path.join(REPO, "results", "hardware_shor_pinned_ibm_fez")
    gl.save_calibration(gl.calibration_from_target(target, plan["used_qubits"]), base + "_calibration.json",
                        backend=plan["backend"], fetched_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                        note="read from the live backend target immediately before submission")
    shutil.copyfile(os.path.join(OUT, "circuits.qpy"), base + "_circuits.qpy")
    inst, usage = guards(svc, need)          # once more, immediately before the job
    job = sampler.run(tqcs)
    job_id = job.job_id()
    stub = os.path.join(REPO, "results", "hardware_shor_pinned_ibm_fez_job.json")
    json.dump({"job_id": job_id, "backend": plan["backend"], "submitted_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "plan": inst.get("plan"), "pricing_type": inst.get("pricing_type"),
               "shots": plan["shots"], "delays_s": plan["delays_s"], "labels": plan["labels"],
               "layout": plan["layout"], "data_qubits": plan["data_qubits"], "ancilla_qubits": plan["ancilla_qubits"],
               "seed_transpiler": plan["seed_transpiler"], "cap_s": CAP_SECONDS,
               "qiskit_version": qiskit.__version__, "qiskit_ibm_runtime_version": qiskit_ibm_runtime.__version__},
              open(stub, "w"), indent=1)
    print(f"submitted job {job_id}; waiting")
    result = job.result()
    counts = {}
    for (t, prep), pub in zip(plan["labels"], result):
        counts[f"{t}_{prep}"] = pub.data.syn.get_counts()
    metrics = job.metrics()
    out = {"meta": {**json.load(open(stub)), "usage": metrics.get("usage", {}),
                    "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())},
           "counts": counts}
    json.dump(out, open(os.path.join(REPO, "results", "hardware_shor_pinned_ibm_fez.json"), "w"), indent=1)
    print(f"done: {metrics.get('usage', {})}")


if __name__ == "__main__":
    {"dry": dry, "submit": submit}[sys.argv[1]]()
