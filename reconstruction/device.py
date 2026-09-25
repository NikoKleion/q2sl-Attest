# simulated toric code device with per-edge noise rates
import numpy as np


def toric_Hz(L):
    # Z plaquette checks: 2L^2 edge qubits, L^2 plaquettes, rank L^2-1
    n = 2 * L * L
    Hh = lambda r, c: (r % L) * L + (c % L)
    Vt = lambda r, c: L * L + (r % L) * L + (c % L)
    rows = []
    for r in range(L):
        for c in range(L):
            row = np.zeros(n, dtype=int)
            for idx in (Hh(r, c), Hh(r + 1, c), Vt(r, c), Vt(r, c + 1)):
                row[idx] ^= 1
            rows.append(row)
    return np.array(rows, dtype=int)


def toric_Hx(L):
    # X vertex checks; vertex (r, c) touches its 4 edges
    n = 2 * L * L
    Hh = lambda r, c: (r % L) * L + (c % L)
    Vt = lambda r, c: L * L + (r % L) * L + (c % L)
    rows = []
    for r in range(L):
        for c in range(L):
            row = np.zeros(n, dtype=int)
            for idx in (Hh(r, c), Hh(r, c - 1), Vt(r, c), Vt(r - 1, c)):
                row[idx] ^= 1
            rows.append(row)
    return np.array(rows, dtype=int)


class Device:
    # same seed gives the same device
    def __init__(self, L=4, p_x=0.03, eta=10.0, hot_frac=0.25, hot_mult=8.0, seed=20260816):
        self.L = L
        self.n = 2 * L * L
        self.Hz = toric_Hz(L)
        self.Hx = toric_Hx(L)
        self.eta = eta
        rng = np.random.default_rng(seed)
        hot = np.ones(self.n)
        hot[rng.random(self.n) < hot_frac] = hot_mult
        self.hot = hot
        self.rx = np.clip(p_x * hot, 0, 0.45)
        self.rz = np.clip(p_x * eta * hot, 0, 0.45)
        self.seed = seed

    def fingerprint(self):
        return {"L": self.L, "eta": self.eta, "rx": self.rx.copy(), "rz": self.rz.copy(),
                "hot_edges": np.where(self.hot > 1)[0]}

    def simulate(self, n_samples, sector="x", rng=None):
        # (syndrome, error) pairs drawn from the device rates
        rng = rng or np.random.default_rng(0)
        H = self.Hz if sector == "x" else self.Hx
        r = self.rx if sector == "x" else self.rz
        E = (rng.random((n_samples, self.n)) < r).astype(np.int8)
        S = (E @ H.T) % 2
        return S.astype(np.int8), E


def device_from_backend(backend, L=4, seed=0):
    # per-edge X and Z rates from a Qiskit backend calibration, needs qiskit
    try:
        props = backend.properties()
        nq = backend.num_qubits
    except Exception as e:
        raise RuntimeError(f"device_from_backend needs a Qiskit backend with calibration: {e}")
    gate_err, read_err = [], []
    for q in range(nq):
        try:
            gate_err.append(float(props.gate_error("sx", q)))
        except Exception:
            gate_err.append(0.01)
        try:
            read_err.append(float(props.readout_error(q)))
        except Exception:
            read_err.append(0.02)
    d = Device(L=L, seed=seed)
    idx = np.arange(d.n) % nq
    d.rx = np.clip(np.array(gate_err)[idx], 1e-4, 0.45)
    d.rz = np.clip(np.array(read_err)[idx], 1e-4, 0.45)
    d.hot = np.ones(d.n)
    d.calibrated_from = getattr(backend, "name", None) or "qiskit-backend"
    return d


def reproduce():
    # checks: Hz and Hx commute, rates show in the data
    d = Device()
    css_ok = int((d.Hz @ d.Hx.T % 2).sum()) == 0
    rng = np.random.default_rng(1)
    Sx, Ex = d.simulate(20000, "x", rng)
    Sz, Ez = d.simulate(20000, "z", rng)
    hot = d.hot > 1
    hot_rate = Ex[:, hot].mean(); cold_rate = Ex[:, ~hot].mean()
    bias_seen = Ez.mean() / max(Ex.mean(), 1e-9)
    ok = css_ok and (hot_rate > 4 * cold_rate) and (bias_seen > 4)
    label = f"device (CSS={css_ok}, hot/cold={hot_rate/max(cold_rate,1e-9):.1f}x, Z/X bias={bias_seen:.1f}x)"
    return (label, round(float(bias_seen), 1), "CSS + hot + bias visible", bool(ok))


if __name__ == "__main__":
    print(reproduce())
    d = Device()
    fp = d.fingerprint()
    print(f"# device L={fp['L']} eta={fp['eta']} hot_edges={len(fp['hot_edges'])}/{d.n}")
    print(f"# rx range {d.rx.min():.3f}-{d.rx.max():.3f}, rz range {d.rz.min():.3f}-{d.rz.max():.3f}")
