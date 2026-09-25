# figures for the documentation, drawn from the saved files in results/
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .estimate import dist_from_counts, null_floor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(ROOT, "results")
FIGDIR = os.path.join(RESDIR, "figures")

INDIGO, TEAL, AMBER, SLATE = "#1F3A5F", "#6B7280", "#111827", "#CBD5E1"
DELAYS = (0.0, 2e-05, 5e-05, 0.0001)
OUTCOMES = 64


def _counts(store, delay, prep):
    for key in (f"({delay}, {prep})", f"{delay}_{prep}"):
        if key in store:
            return store[key]
    raise KeyError((delay, prep))


def _tvd(a, b):
    return 0.5 * float(np.abs(np.asarray(a) - np.asarray(b)).sum())


def _dist(store, delay, prep):
    return dist_from_counts(_counts(store, delay, prep), OUTCOMES)[0]


def _distance(store, delay, prep0=0, prep1=1):
    return _tvd(_dist(store, delay, prep0), _dist(store, delay, prep1))


def _load():
    def counts(name):
        return json.load(open(os.path.join(RESDIR, name)))["counts"]
    return counts("hardware_shor_ibm_fez.json"), counts("hardware_shor_model_today.json"), \
        counts("hardware_shor_model_job_day.json")


def hardware_figure(path):
    """Measured distance against the device model, and the model's state and layout split."""
    measured, today, job_day = _load()
    us = [d * 1e6 for d in DELAYS]
    m = [_distance(measured, d) for d in DELAYS]
    t = [_distance(today, d) for d in DELAYS]
    j = [_distance(job_day, d) for d in DELAYS]
    floor = [null_floor(_dist(measured, d, 0), 16000) for d in DELAYS]

    nonzero = [d for d in DELAYS if d]
    state, layout = [], []
    for d in nonzero:
        twin = dist_from_counts(today[f"{d}_0_on_1"], OUTCOMES)[0]
        state.append(_tvd(twin, _dist(today, d, 1)))
        layout.append(_tvd(_dist(today, d, 0), twin))

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.1))

    ax.fill_between(us, 0, floor, color=SLATE, alpha=0.30, lw=0, label="sampling floor, 16000 shots")
    ax.plot(us, m, "o-", color=AMBER, lw=2.4, ms=7, label="measured on ibm_fez", zorder=3)
    ax.plot(us, t, "s--", color=INDIGO, lw=2, ms=6, label="device model, calibration of the day")
    ax.plot(us, j, "^--", color=TEAL, lw=2, ms=6, label="device model, calibration of the job")
    ax.set_xlabel("idle delay (us)")
    ax.set_ylabel("syndrome distance")
    ax.set_title("Shor [[9,1,3]], six Z generators", fontsize=11)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8.5, frameon=False)
    ax.grid(alpha=0.25, lw=0.6)

    idx = np.arange(len(nonzero))
    bx.bar(idx - 0.19, state, 0.36, color=INDIGO, label="state, two logical states on one layout")
    bx.bar(idx + 0.19, layout, 0.36, color=SLATE, label="layout, one state on two layouts")
    bx.set_xticks(idx, [f"{d * 1e6:.0f}" for d in nonzero])
    bx.set_xlabel("idle delay (us)")
    bx.set_ylabel("syndrome distance")
    bx.set_title("state and layout contributions, under the model", fontsize=11)
    bx.legend(fontsize=8.5, frameon=False)
    bx.grid(alpha=0.25, lw=0.6, axis="y")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return m, t, j


def order_figure(path, gammas=np.logspace(-4, -1, 16)):
    """Population leak against damping strength: the slope is the leak order."""
    from . import expectations as ex
    from .channels import amplitude_damping
    from .experiments import STD
    from .hardware import shor_code

    codes = [STD["repetition"](), STD["code_4_1_2"](), STD["hamming_7"](),
             shor_code(), ex.surface_code_3()]
    colors = ["#9CA3AF", "#6B7280", "#374151", "#111827", "#1F3A5F"]

    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    for code, color in zip(codes, colors):
        leak = [ex.population_leak(code, amplitude_damping(g))[0] for g in gammas]
        slope = np.polyfit(np.log(gammas), np.log(leak), 1)[0]
        ax.loglog(gammas, leak, "o-", color=color, ms=4, lw=1.9,
                  label=f"{code.name}, order {round(slope)}")
    ax.set_xlabel("amplitude damping strength gamma")
    ax.set_ylabel("population leak")
    ax.set_title("leak against damping strength", fontsize=11)
    ax.legend(fontsize=8.5, frameon=False)
    ax.grid(alpha=0.25, lw=0.6, which="both")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _table(name, header, width):
    """The numbers on each line of a saved results table after `header`, until the first blank line after
    rows begin; lines with fewer than `width` numbers are skipped."""
    import re
    number = re.compile(r"-?\d+(?:\.\d+)?(?:e[-+]?\d+)?")
    rows, on = [], False
    for line in open(os.path.join(RESDIR, name), encoding="utf-8"):
        if header in line:
            on = True
            continue
        if not on:
            continue
        if not line.strip():
            if rows:
                break
            continue
        vals = [float(x) for x in number.findall(line)]
        if len(vals) >= width:
            rows.append(vals)
    return rows


def pinned_figure(path):
    """The one-layout Shor run: measured distance, model and floor, and the achieved distance."""
    rec = _table("hardware_shor_pinned_analysis.txt", "six Z generators, 64 outcomes", 10)
    model = _table("hardware_shor_pinned_analysis.txt", "delay us  measured     model", 3)
    us = [r[0] for r in rec]
    measured, floor = [r[3] for r in rec], [r[4] for r in rec]
    achieved, lo, hi = [r[7] for r in rec], [r[8] for r in rec], [r[9] for r in rec]
    mod = [r[2] for r in model]

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.1))
    ax.fill_between(us, 0, floor, color=SLATE, alpha=0.6, lw=0, label="sampling floor, 16000 shots")
    ax.plot(us, measured, "o-", color=INDIGO, lw=2.2, ms=7, label="measured on ibm_fez", zorder=3)
    ax.plot(us, mod, "s--", color=TEAL, lw=1.8, ms=6, label="gate-level device model")
    ax.set_xlabel("idle delay (us)")
    ax.set_ylabel("syndrome distance")
    ax.set_title("Shor [[9,1,3]], six Z generators, one layout", fontsize=11)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8.5, frameon=False)
    ax.grid(alpha=0.25, lw=0.6)

    err = [np.subtract(achieved, lo), np.subtract(hi, achieved)]
    bx.errorbar(us, achieved, yerr=err, fmt="o", color=INDIGO, ms=7, capsize=4, lw=1.6)
    bx.set_xlabel("idle delay (us)")
    bx.set_ylabel("distance achieved on held-out shots")
    bx.set_title("achieved distance with 95% interval", fontsize=11)
    bx.set_ylim(bottom=0)
    bx.grid(alpha=0.25, lw=0.6)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return measured, mod, achieved


def precision_figure(path):
    """L2 distance over gamma cubed on the Hamming code, from the exact engine and the contraction."""
    rows = _table("tensor.txt", "gamma   exact engine    contraction", 2)
    g_exact = [r[0] for r in rows if len(r) == 3]
    exact = [r[1] for r in rows if len(r) == 3]
    g_all = [r[0] for r in rows]
    contracted = [r[-1] for r in rows]
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.semilogx(g_exact, exact, "s--", color=TEAL, lw=1.8, ms=6, label="exact engine")
    ax.semilogx(g_all, contracted, "o-", color=INDIGO, lw=2.2, ms=6, label="tensor network contraction")
    ax.set_yscale("log")
    ax.set_xlabel("amplitude damping strength gamma")
    ax.set_ylabel("L2 distance / gamma^3")
    ax.set_title("Hamming [[7,1,2]]: the ratio at small damping", fontsize=11)
    ax.legend(fontsize=8.5, frameon=False)
    ax.grid(alpha=0.25, lw=0.6, which="both")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return contracted


def records_figure(path):
    """Distance between the two states' Z-check records against idle time, with its floor, 4 million shots."""
    rows = []
    for line in open(os.path.join(RESDIR, "zchecks.txt"), encoding="utf-8"):
        if line.strip().startswith("rotated surface") and "4000000" in line:
            parts = line.split()
            d = 3 if "[[9,1,3]]" in line else int(line.split("d=")[1].split()[0])
            rows.append((d, float(parts[-7]), float(parts[-4]), float(parts[-3])))
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    for d, color in ((3, INDIGO), (5, TEAL)):
        pts = [r for r in rows if r[0] == d]
        ax.plot([r[1] for r in pts], [r[2] for r in pts], "o-", color=color, lw=2, ms=6,
                label=f"rotated surface d={d}")
        ax.plot([r[1] for r in pts], [r[3] for r in pts], ":", color=color, lw=1.2, label=f"d={d}, floor")
    ax.set_yscale("log")
    ax.set_xlabel("idle time (us), T1 163.3 us")
    ax.set_ylabel("distance, number of Z checks that fire")
    ax.set_title("the leak read from 4 million records per state", fontsize=11)
    ax.legend(fontsize=8.5, frameon=False)
    ax.grid(alpha=0.25, lw=0.6, which="both")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return rows


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    m, t, j = hardware_figure(os.path.join(FIGDIR, "fig1_hardware_shor.png"))
    print("fig1_hardware_shor.png")
    print("  measured      ", " ".join(f"{v:.4f}" for v in m))
    print("  model today   ", " ".join(f"{v:.4f}" for v in t))
    print("  model job day ", " ".join(f"{v:.4f}" for v in j))
    order_figure(os.path.join(FIGDIR, "fig2_leak_order.png"))
    print("fig2_leak_order.png")
    m, mod, a = pinned_figure(os.path.join(FIGDIR, "hardware_pinned.png"))
    print("hardware_pinned.png")
    print("  measured ", " ".join(f"{v:.4f}" for v in m))
    print("  model    ", " ".join(f"{v:.4f}" for v in mod))
    print("  achieved ", " ".join(f"{v:.4f}" for v in a))
    c = precision_figure(os.path.join(FIGDIR, "precision.png"))
    print("precision.png", " ".join(f"{v:.6f}" for v in c))
    rows = records_figure(os.path.join(FIGDIR, "records.png"))
    print("records.png", rows)


if __name__ == "__main__":
    main()
