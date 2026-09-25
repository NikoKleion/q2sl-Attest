# BP+OSD decoders on the simulated device with uniform, learned, true and neural priors
import os
import numpy as np

from .device import Device

try:
    import torch
    import torch.nn as nn
    _HAVE_TORCH = True
except Exception:
    _HAVE_TORCH = False
try:
    from ldpc import bposd_decoder
    _HAVE_LDPC = True
except Exception:
    try:
        from ldpc.bposd_decoder import bposd_decoder
        _HAVE_LDPC = True
    except Exception:
        _HAVE_LDPC = False

HERE = os.path.dirname(os.path.abspath(__file__))
DEV = Device()
H = DEV.Hz; M, N = H.shape
if _HAVE_TORCH:
    torch.manual_seed(0)


def make_dec(channel_probs):
    if not _HAVE_LDPC:
        raise RuntimeError("reconstruction decoder needs ldpc (pip install ldpc)")
    cp = np.clip(np.asarray(channel_probs, float), 1e-4, 0.45)
    return bposd_decoder(H, channel_probs=cp.tolist(), max_iter=48,
                         bp_method="ms", osd_method="osd_cs", osd_order=7)


def sim(nsamp, seed):
    return DEV.simulate(nsamp, "x", np.random.default_rng(seed))


def exact_rate_static(channel_probs, S, E):
    # exact-recovery rate of one BP+OSD decoder with a fixed per-qubit prior
    dec = make_dec(channel_probs)
    ok = 0
    for i in range(len(S)):
        ok += int(np.array_equal(dec.decode(S[i].astype(int)), E[i]))
    return ok / len(S)


if _HAVE_TORCH:
    class DecNet(nn.Module):
        def __init__(self):
            super().__init__()
            self.net = nn.Sequential(nn.Linear(M, 128), nn.ReLU(),
                                     nn.Linear(128, 128), nn.ReLU(), nn.Linear(128, N))

        def forward(self, x):
            return self.net(x)

    def train_nn(steps=1500, B=512):
        net = DecNet(); opt = torch.optim.Adam(net.parameters(), 1e-3); lossf = nn.BCEWithLogitsLoss()
        for t in range(steps):
            S, E = sim(B, 100000 + t)
            x = torch.tensor(S.tolist(), dtype=torch.float32)
            y = torch.tensor(E.tolist(), dtype=torch.float32)
            loss = lossf(net(x), y)
            opt.zero_grad(); loss.backward(); opt.step()
        return net

    def nn_probs(net, S):
        x = torch.tensor(S.tolist(), dtype=torch.float32)
        with torch.no_grad():
            p = torch.sigmoid(net(x))
        return np.array(p.tolist())

    def exact_rate_neural(net, S, E):
        # per-sample priors from the MLP, then BP+OSD
        P = nn_probs(net, S)
        ok = 0
        for i in range(len(S)):
            ok += int(np.array_equal(make_dec(P[i]).decode(S[i].astype(int)), E[i]))
        return ok / len(S)


def main(n_eval=600):
    if not _HAVE_LDPC:
        print("# reconstruction decoder needs ldpc (pip install ldpc)"); return
    print(f"# device: toric L={DEV.L}, eta={DEV.eta}, {len(DEV.fingerprint()['hot_edges'])}/{N} hot edges")
    Str, Etr = sim(40000, 7)
    r_hat = Etr.mean(axis=0)
    r_true = DEV.rx
    mean_rate = float(r_true.mean())
    k = len(DEV.fingerprint()["hot_edges"])
    hot_true = set(np.argsort(r_true)[-k:].tolist())
    hot_learned = set(np.argsort(r_hat)[-k:].tolist())
    print(f"# learned hot edges match true: {len(hot_true & hot_learned)}/{k}")
    Se, Ee = sim(n_eval, 999)
    g = exact_rate_static(np.full(N, mean_rate), Se, Ee)
    l = exact_rate_static(r_hat, Se, Ee)
    o = exact_rate_static(r_true, Se, Ee)
    print("# exact error-reconstruction rate on this device's noise:")
    print(f"#   generic  (device-agnostic, uniform) : {g:.3f}")
    print(f"#   learned  (fingerprint from sim)      : {l:.3f}")
    if _HAVE_TORCH:
        net = train_nn()
        torch.save(net.state_dict(), os.path.join(HERE, "models", "device_model.pt"))
        nrl = exact_rate_neural(net, Se, Ee)
        print(f"#   neural   (trained MLP + BP cleanup)  : {nrl:.3f}")
    print(f"#   oracle   (true rates, upper bound)   : {o:.3f}")


if __name__ == "__main__":
    main()
