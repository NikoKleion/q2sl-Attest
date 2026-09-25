# masked-digit predictor for 16 digit numbers, needs torch; luhn_check and gen_pan
import numpy as np

try:
    import torch
    import torch.nn as nn
    _HAVE_TORCH = True
except Exception:
    _HAVE_TORCH = False

rng = np.random.default_rng(0)
if _HAVE_TORCH:
    torch.manual_seed(0)

BINS = [list(map(int, s)) for s in
        ["453201", "510510", "401288", "601100", "371449", "353011", "456789", "412345"]]


def luhn_check(p15):
    s = 0
    for i, d in enumerate(reversed(list(p15))):
        d = int(d)
        if i % 2 == 0:
            d *= 2
            if d > 9:
                d -= 9
        s += d
    return (10 - (s % 10)) % 10


def gen_pan():
    b = BINS[rng.integers(len(BINS))]
    acct = list(rng.integers(0, 10, 9))
    p15 = b + [int(a) for a in acct]
    return np.array(p15 + [luhn_check(p15)], dtype=int)


def dataset(N):
    return np.stack([gen_pan() for _ in range(N)])


def _need_torch(*a, **k):
    raise RuntimeError("reconstruction pattern net needs torch (pip install torch)")


if _HAVE_TORCH:
    class Net(nn.Module):
        def __init__(self):
            super().__init__()
            self.net = nn.Sequential(nn.Linear(16 * 10 + 16, 256), nn.ReLU(),
                                     nn.Linear(256, 256), nn.ReLU(), nn.Linear(256, 16 * 10))

        def forward(self, x):
            return self.net(x).view(-1, 16, 10)

    def encode(batch, known):
        # (B,16) int, (B,16) bool -> (x, digits, known_tensor); one-hot revealed digits plus the known mask
        digits = torch.tensor(batch.tolist(), dtype=torch.long)
        kn = torch.tensor(known.tolist(), dtype=torch.bool)
        oneh = torch.zeros(digits.shape[0], 16, 10)
        oneh.scatter_(2, digits.unsqueeze(-1), 1.0)
        oneh = oneh * kn.unsqueeze(-1).float()
        x = torch.cat([oneh.reshape(digits.shape[0], -1), kn.float()], dim=1)
        return x, digits, kn

    def train(steps=3500, B=256):
        net = Net(); opt = torch.optim.Adam(net.parameters(), lr=1e-3); lossf = nn.CrossEntropyLoss()
        for t in range(steps):
            batch = dataset(B)
            known = rng.random((B, 16)) < rng.uniform(0.1, 0.9)
            x, digits, kn = encode(batch, known)
            unk = ~kn
            if unk.sum() == 0:
                continue
            loss = lossf(net(x)[unk], digits[unk])
            opt.zero_grad(); loss.backward(); opt.step()
        return net

    def predict(net, digits_row, known_row):
        # (16,) int, (16,) bool -> (16,10) probability table
        x, _, _ = encode(digits_row[None, :], known_row[None, :])
        with torch.no_grad():
            probs = torch.softmax(net(x)[0], dim=1)
        return np.array(probs.tolist())
else:
    Net = encode = train = predict = _need_torch


def _train_and_save():
    import os
    net = train()
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models", "pattern_model.pt")
    torch.save(net.state_dict(), path)
    print(f"# saved {path}")


if __name__ == "__main__":
    _train_and_save()
