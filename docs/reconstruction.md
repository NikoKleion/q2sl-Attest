# reconstruction

One posterior over a synthetic 16 digit number from several sources. The full scenario needs torch and
ldpc; `device`, `magic` and `seed` run with numpy.

- `device.py` a toric code device (Kitaev 2003) with per-edge X and Z error rates. `device_from_backend`
  sets the rates from a Qiskit backend calibration.
- `decoder.py` BP+OSD decoding (Panteleev and Kalachev 2021; Roffe, White, Burton and Campbell 2020, through
  the ldpc package) with a uniform, a learned, the true, and a neural per-qubit prior.
- `pattern_net.py` a torch masked-digit predictor.
- `seed.py` residual bits of a 9 digit account drawn from a b bit seed, by brute force and as
  min(29.9, b) - 3.32 per revealed digit.
- `magic.py` stabilizer Renyi entropy M2 (Leone, Oliviero and Hamma 2022).
- `scenario.py` the sources fused on one posterior, with a weight from M2 on the device reads.
- `models/` trained weights for the pattern net and the neural decoder.

Checked against independent answers. `toric_Hx` and `toric_Hz` give two logical qubits and distance L for
L=3 and L=4, the distance by integer programming. Minimum-weight matching (PyMatching) on those checks under
bit-flip noise does better at L=8 than at L=4 at 7 percent and worse at 14 percent, bracketing the published
threshold near 10.3 percent (Wang, Harrington and Preskill 2003). M2 is zero on stabilizer states up to a
three-qubit GHZ state, log2(4/3) on T|+>, additive over a product, and unchanged by a Clifford circuit. The
scenario and seed figures are written to the working directory.

```bash
python -m reconstruction          # needs torch and ldpc
python -m reconstruction.decoder
python -m reconstruction.seed
```
