# analysis of the measured syndrome leak
import json
import math
import os

import numpy as np

import syndrome_leakage as sl
from syndrome_leakage.analyze import population_leak
from syndrome_leakage.channels import amplitude_damping
from syndrome_leakage.core import tvd


def main():
    HERE = os.path.dirname(os.path.abspath(__file__))
    R = json.load(open(os.path.join(os.path.dirname(HERE), "results", "hardware_ibm_fez.json")))
    meta = R["meta"]
    SHOTS = meta["shots"]
    T1_DEV = meta["T1"]
    rng = np.random.default_rng(0)
    code = sl.codes.STANDARD["repetition"]()

    delays = sorted({float(k.split("_")[0]) for k in R["dists"]})
    d = {(float(k.split("_")[0]), int(k.split("_")[1])): np.array(v) for k, v in R["dists"].items()}

    print("1. Sampling error on the measured distance")
    print("   bootstrap: resample both distributions at the shot count and recompute")
    print()
    boot = {}
    for t in delays:
        vals = []
        for _ in range(2000):
            a = rng.multinomial(SHOTS, d[(t, 0)]) / SHOTS
            b = rng.multinomial(SHOTS, d[(t, 1)]) / SHOTS
            vals.append(tvd(a, b))
        boot[t] = (float(np.mean(vals)), float(np.std(vals)))
        print(f"   delay {t*1e6:6.1f} us: measured {tvd(d[(t,0)], d[(t,1)]):.4f}  bootstrap SE {np.std(vals):.4f}")

    print()
    print("2. Null floor")
    print("   distance between two independent samples from one distribution")
    print()
    null = []
    for _ in range(4000):
        a = rng.multinomial(SHOTS, d[(0.0, 0)]) / SHOTS
        b = rng.multinomial(SHOTS, d[(0.0, 0)]) / SHOTS
        null.append(tvd(a, b))
    null = np.array(null)
    z0 = tvd(d[(0.0, 0)], d[(0.0, 1)])
    print(f"   null mean {null.mean():.4f}, 95th pct {np.percentile(null,95):.4f}, max {null.max():.4f}")
    print(f"   measured zero-delay distance {z0:.4f}  ->  {z0/null.mean():.1f}x the null mean")
    print(f"   p(null >= measured) = {(null >= z0).mean():.4f}")

    print()
    print("3. T1 fitted from the leak, against the device value")
    print("   model: measured = a * simulated(gamma(t, T1)), fitted on the nonzero delays")
    print()
    nz = [t for t in delays if t > 0]
    meas = np.array([tvd(d[(t, 0)], d[(t, 1)]) for t in nz])
    best = None
    for T1 in np.linspace(40e-6, 400e-6, 721):
        sim = np.array([population_leak(code, amplitude_damping(1 - math.exp(-t / T1)))[0] for t in nz])
        a = float(np.dot(sim, meas) / np.dot(sim, sim))
        resid = float(np.sum((meas - a * sim) ** 2))
        if best is None or resid < best[0]:
            best = (resid, T1, a, sim)
    resid, T1_fit, a_fit, sim_at_fit = best
    print(f"   fitted T1        {T1_fit*1e6:8.1f} us")
    print(f"   device median T1 {T1_DEV*1e6:8.1f} us")
    print(f"   ratio            {T1_fit/T1_DEV:8.3f}")
    print(f"   attenuation a    {a_fit:8.3f}   (1.0 = no readout blurring)")
    print(f"   rms residual     {math.sqrt(resid/len(nz)):8.4f}")

    print()
    print("4. MEASURED vs SIMULATED AT THE DEVICE T1")
    print()
    print(f"   {'delay us':>9} {'measured':>10} {'+/-SE':>7} {'simulated':>11} {'diff':>8} {'ratio':>7}")
    for t in nz:
        g = 1 - math.exp(-t / T1_DEV)
        s = population_leak(code, amplitude_damping(g))[0]
        m = tvd(d[(t, 0)], d[(t, 1)])
        print(f"   {t*1e6:9.1f} {m:10.4f} {boot[t][1]:7.4f} {s:11.4f} {m-s:+8.4f} {m/s:7.3f}")

    print()
    print("5. THE SAME RECORDS THROUGH estimate.py")
    print("   permutation test for any difference, and the distance a fitted likelihood ratio achieves")
    print("   on held-out shots")
    print()
    from syndrome_leakage import estimate as es  # noqa: E402

    rows = [(f"{t*1e6:.0f} us", es.leak_from_dists(d[(t, 0)], SHOTS, d[(t, 1)], SHOTS, boots=2000,
                                                   splits=16, seed=1)) for t in delays]
    print(es.report(rows, "   ibm_fez records"))


if __name__ == "__main__":
    main()
