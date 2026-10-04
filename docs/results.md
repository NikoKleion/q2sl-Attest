# Results

Validation runs for the tools in `syndrome_leakage/`. Each section runs one tool and checks it against a
control: a closed form, a published equation, a second engine, or a hardware measurement.

```bash
python -m syndrome_leakage.experiments                 # every run
python -m syndrome_leakage.experiments device coherent # named runs
python -m syndrome_leakage.experiments --save          # write the output to results/
python -m syndrome_leakage.figures                     # the figures, from the saved output
```

Run names: `selftest`, `attack`, `leak_order`, `coherent`, `structure`, `device`, `held_memory`,
`worst_pair`, `pauli_boundary`, `audit`, `kl_blocks`, `surface`, `protection`, `shor_circuits`,
`shor_noise`, `estimator`, `disorder`, `readout`, `rounds`, `shor_device_model`, `tensor`, `sampled`,
`zchecks`, `decoded`, `regions`, `circuit_parts`, `circuit_parts_surface`, `drop_sign`. `shor_circuits` needs qiskit and qiskit-aer; `shor_noise` also needs qiskit-ibm-runtime for the
fake backend and takes about seven minutes. `shor_device_model` reads the circuits and calibration saved
with the `hardware_shor` measurement. `sampled` needs stim and takes about 25 minutes, most of it on the two
bivariate bicycle codes. `zchecks` uses scipy for the integer programs, ldpc for the decoded parity and stim
for the comparison with `sampled`, and takes about a quarter of an hour. `regions` uses scipy for the integer
programs and takes about eight minutes. `circuit_parts` reads the circuits and calibration saved with the
`hardware_shor_pinned` measurement and needs qiskit-aer; it samples for about half an hour where
`results/circuit_parts_counts.json` is absent and computes from that file where it is present.
`circuit_parts_surface` does the same for the surface code with qiskit-ibm-runtime's fake backend.

All simulation is exact and deterministic except the attack simulation, which samples syndrome records with
a fixed seed. Where a simulated run uses device parameters they are the medians of IBM `ibm_marrakesh`,
T1 163.3 us and T2 77.0 us, with a 1 us syndrome-extraction cycle.

Three sections, `hardware_shor`, `hardware_shor_pinned` and `hardware`, are measured on IBM `ibm_fez`. They
are not part of `experiments.py`, since they need a provider account and queue time; their data is saved
under `results/`.

---

## Setup

### Synthetic

Exact density-matrix simulation. One single-qubit Kraus channel is applied independently to every physical
qubit of the code. Syndrome probabilities are `Tr(P_s rho)` with the stabilizer projectors constructed in
full, so the stabilizer measurement is ideal: no ancillas, no gate error, no measurement error. Logical
states come from `code.logical_state(theta, phi)`. Every quantity is deterministic; the only sampling in
the package is the attack simulation, which draws syndrome records and carries a fixed seed.

This is the setting for `selftest`, `attack`, `leak_order`, `coherent`, `structure`, `device`,
`held_memory`, `worst_pair`, `pauli_boundary`, `audit`, `kl_blocks`, `surface`, `protection` and
`estimator`, and for the first part of `regions`. The `device` run differs from the others only in taking its channel from measured T1 and T2,
where the others use a free parameter. The `disorder` run is the one place where the channel differs from
qubit to qubit, taking each one's T1 from a backend. The `tensor` run keeps the same channel model and
computes the L2 distance between the two syndrome distributions by contraction instead of simulating a
state. The `sampled` run draws syndrome records from a stabilizer simulator under T1 and T2 relaxation with
T2 <= T1, with fixed seeds, and compares them through a statistic. The `zchecks` run treats the Z checks under
relaxation as a classical process, which needs no quantum state and no condition on T2.

### Aer

`shor_circuits` runs the Shor-code circuits of `hardware.py` on the local Qiskit Aer simulator, with ideal
gates and readout. `shor_noise` runs the same circuits under the calibration of the `FakeFez` backend of
qiskit-ibm-runtime: the transpiler picks the layout and routing on the 156-qubit device, the circuit is
rebuilt on the qubits it uses, and each of those qubits carries its own T1, T2, gate errors and readout
error. Both simulators run locally. `circuit_parts` runs the circuits of the pinned ibm_fez job under the
model of its own calibration, with the encoder, the window, the extraction or the readout made noiseless.
`circuit_parts_surface` does it for the distance 3 surface code on the `FakeFez` calibration.

### Hardware

IBM `ibm_fez`, 4000 shots per circuit. A 3-qubit repetition code on five qubits: three data qubits and two
ancillas. Each circuit prepares the logical state, idles for a fixed delay, extracts both Z-stabilizers,
and measures the two ancillas.

| step | operation |
|---|---|
| prepare | identity for \|0_L>, X on all three data qubits for \|1_L> |
| idle | `delay(t)` on the data qubits, t in {0, 20, 50, 100} us |
| extract ZZI | CX(d0, a0), CX(d1, a0) |
| extract IZZ | CX(d1, a1), CX(d2, a1) |
| measure | both ancillas, giving one of four syndromes |

Transpiled at optimization level 1 with `scheduling_method="alap"` and `seed_transpiler=7`; scheduling is
required because the circuits carry explicit delays. Layout and routing come from the transpiler.

The hardware setting adds gate error on the four CX gates, readout error on the ancillas, and a
transpiler-chosen layout, all of which are absent from the synthetic setting. The delay sweep varies the
amplitude-damping contribution with idle time and leaves the additions above fixed.

---

## selftest

The analytic leak order comes from the stabilizer weight structure; the exact density-matrix
simulation is the reference. This checks one against the other on every standard code.

Control: the phase leak must stay below 1e-9, and codes that protect the logical population must report no
leak.

| code | verdict | analytic order | measured slope |
|---|---|---|---|
| repetition | leaks | 1 | 0.96 |
| [[4,1,2]] | leaks | 2 | 1.93 |
| Hamming [[7,1,2]] | leaks | 3 | 2.90 |
| Steane [[7,1,3]] | protected | - | - |
| five-qubit [[5,1,3]] | protected | - | - |

The measured slope of the leak against noise strength matches the analytic order on every leaking code.

---

## pauli_boundary

The maximum pairwise syndrome TVD over 32 logical states spanning the Bloch sphere. Zero means the syndrome
distribution is the same whichever state is encoded. For a Pauli channel, Wagner, Kampermann, Bruss and
Kliesch (Quantum 6, 809, 2022, eq. 40) write each stabilizer expectation as a function of the Pauli error
distribution alone.

| code | depolarizing | dephasing | amplitude damping | coherent diagonal |
|---|---|---|---|---|
| repetition | 9.7e-17 | 0 | 4.80e-01 | 0 |
| [[4,1,2]] | 1.3e-16 | 8.3e-17 | 1.02e-01 | 3.19e-01 |
| Hamming [[7,1,2]] | 3.9e-16 | 2.6e-16 | 1.74e-02 | 7.63e-03 |

Pauli channels give machine-precision zero across the sphere; the two non-Pauli channels give the values
in the table.

Steane and the five-qubit code also give machine-precision zero under amplitude damping. Over the same
32 states their maximum pairwise syndrome TVD stays
at or below 2.0e-16 at gamma 0.1, 0.3, 0.5, 0.7, 0.9 and 0.99. The test
`test_protected_codes_hold_at_every_noise_strength` pins this.

---

## leak_order

The analytic leak order beside the amplitude-damping population distance, the minimum weight w for which
some weight-w number operator has a different expectation on the two codewords. The second quantity comes
from the codewords alone and is independent of the stabilizer-group algebra the first uses.

For the [[4,1,2]] code, Leung, Nielsen, Chuang and Yamamoto (1997) find the detection probabilities vary
across the codespace at order gamma^2. The analytic order for that code is 2.

| code | analytic order | AD population distance | equal |
|---|---|---|---|
| repetition | 1 | 1 | yes |
| [[4,1,2]] | 2 | 2 | yes |
| Hamming [[7,1,2]] | 3 | 3 | yes |
| Steane [[7,1,3]] | none | 3 | no |
| five-qubit [[5,1,3]] | none | 5 | no |

The two quantities agree on the three leaking codes and diverge on Steane and the five-qubit code, whose
codewords separate at finite weight while the syndrome distribution does not change. The analytic order
and the codeword distance are different quantities.

The coefficient of the leading term, leak divided by gamma^order at gamma = 1e-4, beside the number of
minimal sets the order comes from:

| code | order | coefficient | minimal sets |
|---|---|---|---|
| 3-qubit repetition | 1 | 2.9997 | 3 |
| [[4,1,2]] | 2 | 3.9992 | 4 |
| Hamming [[7,1,2]] | 3 | 3.9988 | 4 |
| Shor [[9,1,3]] | 3 | 26.9917 | 27 |
| rotated surface [[9,1,3]] | 3 | 6.9981 | 8 |

The two columns agree on four of these five codes and disagree on the surface code, which has eight
minimal sets and a coefficient of seven. On the [[4,2,2]] code the coefficient is 4.0 against 2 sets, and
on the L=2 toric code it is 2.0 against 4, so the agreement above is a coincidence of those codes rather
than a rule. The order is available from the stabilizer strings at any size; the coefficient is not, and
the column is here as a target for a derivation that has not been done.

---

## attack

The error rate of a likelihood-ratio test on the syndrome record, beside the Chernoff rate and the
Bhattacharyya bound. The likelihood-ratio test is the most powerful test for independent shots (Neyman
and Pearson, 1933).
Hamming [[7,1,2]] under amplitude damping at gamma 0.2, syndrome TVD 1.741e-02, Chernoff exponent
0.001002 nats at s* = 0.485, 200000 trials per point.

Control: the one-shot closed form. For a single shot at equal priors the minimum error is exactly
(1 - TVD)/2 = 0.4913, and the simulation returns 0.4908 +/- 0.0011.

| rounds | attack error (200k trials, 1 SE) | Chernoff rate | Bhattacharyya bound |
|---|---|---|---|
| 1 | 0.4908 +/- 0.0011 | 0.4995 | 0.4995 |
| 10 | 0.4518 +/- 0.0011 | 0.4950 | 0.4950 |
| 50 | 0.3782 +/- 0.0011 | 0.4756 | 0.4756 |
| 200 | 0.2625 +/- 0.0010 | 0.4092 | 0.4093 |

At 200 rounds the bound gives 0.409 and the attack reaches 0.2625. Rounds of syndrome data to reach 1
percent error:

| code | gamma 0.1 | gamma 0.2 | gamma 0.3 |
|---|---|---|---|
| repetition | 15 | 8 | 5 |
| [[4,1,2]] | 282 | 88 | 50 |
| Hamming [[7,1,2]] | 26510 | 4595 | 1991 |

Steane and the five-qubit code reach no round count.

Precision: the standard error of the simulated attack is sqrt(e(1-e)/trials). At 200000 trials it is about
0.001, so three decimals are reported.

---

## worst_pair

The leakage measures compare |0_L> against |1_L>. This searches a Bloch grid of 58 logical states and
maximises the Chernoff exponent over all pairs.

Control: protected codes return numerical noise only.

| channel | code | worst pair | population axis | axis is worst |
|---|---|---|---|---|
| amplitude damping | repetition | 0.654 | 0.654 | yes |
| amplitude damping | [[4,1,2]] | 0.053 | 0.053 | yes |
| amplitude damping | Hamming | 0.0010 | 0.0010 | yes |
| coherent diagonal | [[4,1,2]] | 0.384 | 0.384 | yes |
| coherent diagonal | Hamming | 0.0003 | 0.0003 | yes |

On this grid the population axis gives the largest value on every leaking code under both channels.

---

## coherent

The same measures under a coherent Z-rotation applied to every qubit. Population TVD at theta = 0.3:

| code | TVD |
|---|---|
| repetition | 0 |
| [[4,1,2]] | 3.19e-01 |
| Hamming [[7,1,2]] | 7.63e-03 |
| Steane, five-qubit | ~1e-16 |

Repetition has only Z-type stabilizers; the diagonal error commutes with them and the syndrome does not
move.

Controls: dephasing and depolarizing both give exactly 0. The leak is second order in the angle,
TVD/theta^2 holding near 3.99 as theta goes to zero (0.00997 at 0.05, 0.0395 at 0.1, 0.152 at 0.2).

For [[4,1,2]] one of the two logical states never flips XXXX, so the TVD is sin^2(2 theta) and the
Chernoff exponent is -2 ln cos(2 theta): 0.1645 nats per shot at theta = 0.2, where `rounds_for_error`
gives 29 rounds to 1 percent error. `tests/test_closed_forms.py` checks both closed forms, and equation 93
of Hu, Liang and Calderbank for Steane.

A coherent Z-rotation is a diagonal unitary, so this case falls under Hu, Liang and Calderbank
(arXiv:2109.13481), who give these probabilities in closed form and the condition for them not to depend
on the encoded state. The amplitude-damping runs are outside that setting. `hlc.py` implements their
equations 45, 48 and 91, and `tests/test_hlc.py` checks the closed form against this exact simulation on the
four CSS codes at theta 0.1, 0.3 and 0.7 for four logical states, to 1e-12.

---

## structure

Leak order against code distance, for CSS codes. Distances from `css.code_distance`.

| code | d | analytic order | measured slope |
|---|---|---|---|
| [[4,1,2]] | 2 | 2 | - |
| Hamming [[7,1,2]] | 2 | 3 | - |
| Steane [[7,1,3]] | 3 | none | - |
| Hamming-CSS self-dual | 3 | none | - |
| Shor [[9,1,3]] | 3 | 3 | 2.89 |
| rotated surface [[9,1,3]] | 3 | 3 | 2.91 |

At distance 3, Steane and the self-dual Hamming-CSS code show no leak, and Shor and the rotated surface
code leak at order 3. At distance 2, [[4,1,2]] leaks at order 2 and the built-in Hamming code at order 3.
The distance does not set the leak order.

Steane and the self-dual Hamming-CSS code have Z generators equal to their X generators with X replaced by
Z. The built-in Hamming code has X generators of weights 4, 4, 4 and Z generators of weights 7, 3, 3; X on
qubits 0 and 4 commutes with every generator and maps |0_L> to |1_L>, so its distance is 2. Shor pairs two
weight-6 X generators with six weight-2 Z generators.

The last two rows come from the stabilizer strings alone, with no state simulated, so they reach codes the
engines cannot hold:

| code | n | generators | distance | order |
|---|---|---|---|---|
| rotated surface d=3 | 9 | 8 | 3 | 3 |
| rotated surface d=5 | 25 | 24 | 5 | 5 |

The order equals the distance on both, so a larger surface code leaks at a higher order in the noise
strength. `rotated_surface_code(d)` builds the family for any odd d, and the search over the Z-side span
reaches d=5 in under a second; d=7 has a span of 2^25 elements and is refused rather than attempted.
The `tensor` run carries the order to d=11 through the L2 distance.

The nine-qubit codes are past the projector limit, so their measured slopes come from the expectation
engine of the `surface` run.

---

## tensor

The L2 distance between the syndrome distributions of |0_L> and |1_L>, computed by tensor network
contraction in `tensor.l2_leak` rather than by simulating a state. Half of it bounds the leak from below.
Requires cotengra.

Control: the exact engine, `analyze.l2_leak`, on every code it reaches. L2 distance at amplitude
damping 0.05:

| code | exact | contracted |
|---|---|---|
| 3-qubit repetition | 1.645448e-01 | 1.645448e-01 |
| [[4,1,2]] | 6.381639e-03 | 6.381639e-03 |
| Hamming [[7,1,2]] | 1.211094e-04 | 1.211094e-04 |
| Shor [[9,1,3]] | 5.568814e-04 | 5.568814e-04 |
| rotated surface [[9,1,3]] | 2.024389e-04 | 2.024389e-04 |
| toric L=2, logical 0 | 2.556069e-03 | 2.556069e-03 |
| toric L=2, logical 1 | 2.556069e-03 | 2.556069e-03 |

Under amplitude damping at 0.2 and 0.05 every code agrees with the exact engine to a relative 1e-12.
Depolarizing noise at 0.1 gives exactly 0.0 on Shor.

A channel that moves the X part of a Pauli operator takes the general network, with separate X variables
on the three copies. A coherent X rotation of 0.3 is such a channel and gives zero in both engines on
every code. Amplitude damping at 0.2 followed by that rotation gives a nonzero distance:

| code | exact | contracted |
|---|---|---|
| 3-qubit repetition | 3.775477e-01 | 3.775477e-01 |
| [[4,1,2]] | 3.439424e-02 | 3.439424e-02 |
| Hamming [[7,1,2]] | 2.272611e-03 | 2.272611e-03 |
| Shor [[9,1,3]] | 6.727061e-03 | 6.727061e-03 |
| rotated surface [[9,1,3]] | 1.894037e-03 | 1.894037e-03 |
| toric L=2, logical 0 | 1.191234e-02 | 1.191234e-02 |
| toric L=2, logical 1 | 1.191234e-02 | 1.191234e-02 |

Every code agrees to a relative 1e-12.

Precision. On the Hamming code the L2 distance over gamma cubed must settle to a constant as the damping
falls:

| damping | exact engine | contraction |
|---|---|---|
| 1e-3 | 1.114906 | 1.114906 |
| 1e-4 | 1.117687 | 1.117721 |
| 1e-5 | 1.252118 | 1.118003 |
| 1e-6 | 15.222137 | 1.118031 |
| 1e-8 | | 1.118034 |
| 1e-10 | | 1.118034 |

The exact engine subtracts two distributions that agree to near machine precision. Its ratio has moved
in the fourth digit at 1e-4 and is wrong in the first at 1e-6. The contraction keeps only the coset holding the active logical Z, so the difference is
one object, and the ratio settles to 1.118034.

Reach. L2 distance at damping 0.05, and the exponent of the L2 distance against the damping, fitted over
1e-4, 1e-5 and 1e-6:

| code | n | k | L2 at 0.05 | exponent |
|---|---|---|---|---|
| rotated surface d=5 | 25 | 1 | 2.6684e-07 | 5.000 |
| rotated surface d=7 | 49 | 1 | 2.5047e-10 | 7.000 |
| rotated surface d=9 | 81 | 1 | 1.8129e-13 | 8.999 |
| rotated surface d=11 | 121 | 1 | 1.0384e-16 | 10.999 |
| hypergraph product rep4 | 25 | 1 | 1.3616e-06 | 4.000 |
| hypergraph product rep5 | 41 | 1 | 2.5993e-08 | 5.000 |
| toric L=3, logical 0 | 18 | 2 | 5.3895e-05 | 3.000 |
| rotated surface d=13 | 169 | 1 | out of reach | |
| toric L=4, logical 0 | 32 | 2 | out of reach | |
| bivariate bicycle [[18,4,4]], logical 0 | 18 | 4 | out of reach | |

The exact engine stops at about 14 qubits, so every row here is past it. When the difference of the two
distributions starts at order w in the damping, both its L1 and its L2 norm start at order w, so the fitted
exponent is the leak order. On every surface code it equals the distance, which the `structure` run
establishes by enumeration only up to d=5. The hypergraph product of two length-r repetition codes is the
planar surface code of distance r, and its exponent is r; toric L=3 has distance 3 and exponent 3.

The last three rows need an intermediate tensor of more than 2^26 entries, the default limit, and are not
contracted. A code with more than one logical qubit carries every other logical Z as a variable across its
support in two of the three copies, and a periodic boundary widens the network again. The contraction plan
comes from a seeded search of fixed length, so a given cotengra version returns the same plan, and the same
verdict, on every run.

---

## sampled

Syndrome records drawn from a stabilizer simulator under T1 and T2 relaxation, in `sampled.py`, and the
leak read from them through a statistic. With T2 <= T1 the relaxation channel is exactly reset to |0>
with probability 1 - exp(-t/T1) and Z with probability (1 - exp(-t (1/T2 - 1/T1))) / 2, so Stim runs it
on codes of any size. Every random outcome comes from numpy with a fixed seed. Requires stim.

Control: the channel. The reset and Z mixture equals `channel_from_t1t2` at four settings of T1, T2 and
idle time, with the largest Choi difference below 1e-14. T2 > T1 is refused.

Control: the exact engine. T1 100 us, T2 60 us and 30 us of idle give reset probability 0.2592 and Z
probability 0.0906. With 200000 shots per state, distances between distributions of the full syndrome:

| code | exact leak | sampled leak | 0_L from exact | 1_L from exact | floor | p, 0_L | p, 1_L |
|---|---|---|---|---|---|---|---|
| 3-qubit repetition | 0.5760 | 0.5758 | 0.0000 | 0.0019 | 0.0015 | 1.000 | 0.236 |
| [[4,1,2]] | 0.1475 | 0.1468 | 0.0007 | 0.0012 | 0.0020 | 0.997 | 0.859 |
| 5-qubit [[5,1,3]] | 0.0000 | 0.0047 | 0.0035 | 0.0036 | 0.0033 | 0.419 | 0.305 |
| Shor [[9,1,3]] | 0.1911 | 0.1917 | 0.0098 | 0.0107 | 0.0098 | 0.102 | 0.103 |
| rotated surface [[9,1,3]] | 0.0496 | 0.0519 | 0.0122 | 0.0118 | 0.0119 | 0.323 | 0.554 |
| toric L=2, logical 1 | 0.0737 | 0.0714 | 0.0059 | 0.0075 | 0.0066 | 0.865 | 0.094 |

The floor is the mean distance between a sample of 200000 shots and the distribution it came from, and
each record lies at about that distance from its exact distribution. p is the chance that such a sample
lies at least as far; the smallest of the twelve is 0.094. The repetition code's |0_L> is |000>, which
relaxation leaves in place, so its record is exact. The five-qubit code is not CSS and is protected: its
two records differ by 0.0047, which is the reading two samples of one distribution give at this size.

Past the exact engines. T1 163.3 us and T2 77.0 us at three idle times, 20000 shots per state, read
through the number of Z generators that read -1. `achieved` is the distance an observer reading that
statistic achieves on held-out shots, with its null subtracted:

| code | n | k | idle us | reset | tvd | null | p | achieved | exact leak | L2 / 2 |
|---|---|---|---|---|---|---|---|---|---|---|
| rotated surface d=3 | 9 | 1 | 20 | 0.1153 | 0.0028 | 0.0059 | 0.883 | 0.0000 | 0.0074 | 7.55e-04 |
| rotated surface d=3 | 9 | 1 | 60 | 0.3075 | 0.0412 | 0.0070 | 0.000 | 0.0374 | 0.0676 | 5.66e-03 |
| rotated surface d=3 | 9 | 1 | 120 | 0.5204 | 0.0877 | 0.0073 | 0.000 | 0.0852 | 0.1088 | 8.75e-03 |
| rotated surface d=5 | 25 | 1 | 20 | 0.1153 | 0.0080 | 0.0094 | 0.639 | 0.0000 | | 1.10e-06 |
| rotated surface d=5 | 25 | 1 | 60 | 0.3075 | 0.0092 | 0.0107 | 0.659 | 0.0000 | | 1.01e-05 |
| rotated surface d=5 | 25 | 1 | 120 | 0.5204 | 0.0118 | 0.0106 | 0.328 | 0.0000 | | 1.83e-05 |
| rotated surface d=7 | 49 | 1 | 20 | 0.1153 | 0.0147 | 0.0121 | 0.181 | 0.0017 | | 4.24e-10 |
| rotated surface d=7 | 49 | 1 | 60 | 0.3075 | 0.0111 | 0.0135 | 0.771 | 0.0000 | | 1.27e-09 |
| rotated surface d=7 | 49 | 1 | 120 | 0.5204 | 0.0111 | 0.0131 | 0.758 | 0.0000 | | 2.29e-09 |
| BB [[72,12,6]], logical 0 | 72 | 12 | 20 | 0.1153 | 0.0131 | 0.0153 | 0.773 | 0.0000 | | |
| BB [[72,12,6]], logical 0 | 72 | 12 | 60 | 0.3075 | 0.0146 | 0.0155 | 0.601 | 0.0000 | | |
| BB [[72,12,6]], logical 0 | 72 | 12 | 120 | 0.5204 | 0.0192 | 0.0148 | 0.091 | 0.0078 | | |
| BB [[144,12,12]], logical 0 | 144 | 12 | 20 | 0.1153 | 0.0184 | 0.0196 | 0.655 | 0.0000 | | |
| BB [[144,12,12]], logical 0 | 144 | 12 | 60 | 0.3075 | 0.0171 | 0.0185 | 0.670 | 0.0000 | | |
| BB [[144,12,12]], logical 0 | 144 | 12 | 120 | 0.5204 | 0.0201 | 0.0178 | 0.206 | 0.0008 | | |

On the d=3 surface code the statistic shows the leak at 60 and 120 us, where the observer achieves 0.0374
and 0.0852 against exact leaks of 0.0676 and 0.1088. At 20 us, where the exact leak is 0.0074, the
statistic's distance of 0.0028 is under its null. At 60 and 120 us the contracted L2 / 2 is 5.66e-03 and
8.75e-03.

On surface d=5 and d=7 and on both bivariate bicycle codes the statistic shows no difference between the
two logical states at any idle time; the smallest p-value is 0.091, on [[72,12,6]] at 120 us. This
describes an observer reading this statistic from 20000 records per state, and does not bound the leak
from above. The contracted L2 / 2 on d=5 and d=7 is at most 1.83e-05, far below what 20000 records
resolve. 120 us of idle is 0.73 T1; over a 1 us syndrome round the reset probability is 0.0061.

---

## zchecks

The Z checks of a CSS code under relaxation, in `zchecks.py`, treated as a classical process: the resets
leave the X error x AND R, for x a Z-basis outcome of the logical state and R the reset set, and the Z
syndrome is Hz (x AND R). The leak order comes from an integer program.

Control: T2. At T1 100 us and 30 us of idle, the Z-check syndromes under relaxation with T2 of 60, 100,
150 and 200 us equal those under reset alone on eight codes, to 1e-12. The Z checks see only populations,
and relaxation moves populations the same way at every T2.

Control: the exact engine, the enumeration and the Stim sampler. At reset probability 0.3 the process
equals the exact engine's Z-check marginal to 1e-12 on each code, and 200000 sampled shots follow the
enumeration of the process:

| code | exact to 1e-12 | p, 0_L | p, 1_L |
|---|---|---|---|
| 3-qubit repetition | True | 1.000 | 0.219 |
| [[4,1,2]] | True | 0.152 | 0.091 |
| Hamming [[7,1,2]] | True | 0.808 | 0.731 |
| Steane [[7,1,3]] | True | 0.098 | 0.683 |
| Shor [[9,1,3]] | True | 0.828 | 0.259 |
| rotated surface [[9,1,3]] | True | 0.877 | 0.921 |
| toric L=2, logical 0 | True | 0.613 | 0.580 |
| toric L=2, logical 1 | True | 0.613 | 0.042 |

One of the sixteen p-values falls below 0.05, as one in twenty does by chance. Against the Stim sampler of
the `sampled` run on [[72,12,6]], at reset 0.3 with 20000 shots each, the distributions of the number of Z
checks that fire give permutation p-values of 0.075 and 0.869 for the two states. The process takes 0.12 s
for those shots and the Stim sampler 35 s.

Control: the enumeration of the Z-side span. The integer-program leak order against `analytic_leak`:

| code | integer program | enumeration |
|---|---|---|
| 3-qubit repetition | 1 | 1 |
| [[4,1,2]] | 2 | 2 |
| Hamming [[7,1,2]] | 3 | 3 |
| Steane [[7,1,3]] | protected | protected |
| Shor [[9,1,3]] | 3 | 3 |
| rotated surface [[9,1,3]] | 3 | 3 |
| toric L=2, logical 0 | 2 | 2 |
| toric L=2, logical 1 | 2 | 2 |
| rotated surface d=5 | 5 | 5 |
| toric L=3, logical 0 | 3 | 3 |
| toric L=4, logical 0 | 4 | 4 |
| hypergraph product rep4 | 4 | 4 |
| BB [[18,4,4]], logical 0 | 4 | 4 |

Equal on every code, including Steane, which both find protected, and the Hamming [[7,1,2]] code, whose
order of 3 exceeds its distance of 2.

Control: known and published distances. The Z distance by integer programming, each proven optimal by the
solver:

| code | integer program | known |
|---|---|---|
| rotated surface d=7 | 7 | 7 |
| toric L=4 | 4 | 4 |
| hypergraph product rep5 | 5 | 5 |
| BB [[72,12,6]] | 6 | 6 |

The [[72,12,6]] distance is the published one [Bravyi et al. 2024]. qLDPC's `BBCode` builds the same check
matrices as the construction here for both bicycle codes, and its randomized distance bounds give 6 and 12.

The leak order past every enumeration, proven optimal:

| code | leak order |
|---|---|
| BB [[72,12,6]], logical 0 | 6 |
| BB [[144,12,12]], logical 0 | 12 |

Which statistic carries the leak. Exact distances between the two states' distributions of the full Z
syndrome and of three statistics of it, by enumerating the process; `parity` is the parity of the error
BP+OSD decodes against the active Z logical, and `both` is the parity together with the number of Z checks
that fire:

| code | reset | full | z_weight | parity | both |
|---|---|---|---|---|---|
| Shor [[9,1,3]] | 0.1 | 0.01968 | 0.01166 | 0.00073 | 0.01385 |
| Shor [[9,1,3]] | 0.3 | 0.25005 | 0.14818 | 0.00926 | 0.17596 |
| Shor [[9,1,3]] | 0.5 | 0.42188 | 0.25000 | 0.01562 | 0.29688 |
| rotated surface [[9,1,3]] | 0.1 | 0.00510 | 0.00365 | 0.00073 | 0.00365 |
| rotated surface [[9,1,3]] | 0.3 | 0.06483 | 0.04630 | 0.00926 | 0.04630 |
| rotated surface [[9,1,3]] | 0.5 | 0.10938 | 0.07812 | 0.01562 | 0.07812 |
| toric L=2, logical 0 | 0.1 | 0.01620 | 0.01381 | 0.00120 | 0.01381 |
| toric L=2, logical 0 | 0.3 | 0.08820 | 0.06671 | 0.01074 | 0.06671 |
| toric L=2, logical 0 | 0.5 | 0.12500 | 0.09375 | 0.01562 | 0.09375 |

The number of Z checks that fire carries 59 to 85 percent of the Z-syndrome leak on these codes, and the
decoded parity 4 to 14 percent. The parity registers the leak only when the resets cover a whole Z logical
and the decoder then lands on the right side of it.

The leak read from records through the number of Z checks that fire, at T1 163.3 us; T2 does not enter:

| code | order | idle us | reset | shots | tvd | null | p | achieved |
|---|---|---|---|---|---|---|---|---|
| rotated surface [[9,1,3]] | 3 | 20 | 0.1153 | 4000000 | 0.00521 | 0.00042 | 0.000 | 0.00518 |
| rotated surface [[9,1,3]] | 3 | 60 | 0.3075 | 4000000 | 0.04816 | 0.00051 | 0.000 | 0.04794 |
| rotated surface [[9,1,3]] | 3 | 120 | 0.5204 | 4000000 | 0.07798 | 0.00051 | 0.000 | 0.07773 |
| rotated surface d=5 | 5 | 20 | 0.1153 | 4000000 | 0.00039 | 0.00067 | 0.917 | 0.00000 |
| rotated surface d=5 | 5 | 60 | 0.3075 | 4000000 | 0.00285 | 0.00076 | 0.000 | 0.00249 |
| rotated surface d=5 | 5 | 120 | 0.5204 | 4000000 | 0.00532 | 0.00075 | 0.000 | 0.00501 |
| rotated surface d=7 | 7 | 120 | 0.5204 | 10000000 | 0.00066 | 0.00059 | 0.278 | 0.00007 |
| BB [[72,12,6]], logical 0 | 6 | 120 | 0.5204 | 40000000 | 0.00034 | 0.00033 | 0.407 | 0.00001 |
| BB [[144,12,12]], logical 0 | 12 | 120 | 0.5204 | 20000000 | 0.00074 | 0.00057 | 0.040 | 0.00017 |

The d=5 surface code shows the leak at 60 and 120 us. There the observer achieves 0.00249 and 0.00501 on
held-out shots, where the contracted L2 / 2 under the same channel is 2.16e-05 and 2.40e-05. Surface d=7
at 10 million shots and [[72,12,6]] at 40 million show no difference.

The [[144,12,12]] row reads p = 0.040. Three other runs of the same code at 120 us, each with its own seed,
read p = 0.098 at 4 million shots, 0.453 at 20 million and 0.812 at 40 million, and [[72,12,6]], whose
order is half as high, shows nothing at 40 million. Four rows of the table resolve no leak, and among four
such tests a p-value of 0.040 or less arises by chance about 15 percent of the time. It is not a detection.

These describe an observer counting the Z checks that fire in that many records. They do not bound the
leak from above.

---

## device

Leak and rounds-to-error with the channel built from measured device parameters.

Control: setting T1 to infinity, leaving only the T2 phase-damping part, drops the repetition leak to
3.00e-15.

| code | population TVD per round | rounds to 1 percent |
|---|---|---|
| repetition | 1.820e-02 | 251 |
| [[4,1,2]] | 1.473e-04 | 62539 |
| Hamming [[7,1,2]] | 8.936e-07 | 9.6e9 |
| Steane, five-qubit | ~1e-16 | none |

These counts are in the independent-shot regime, where each round is a fresh preparation of the same
logical state.

---

## held_memory

A logical memory holds one state across rounds and extracts syndromes from it each round.
This evolves the ensemble exactly across rounds, with recovery applied each round, at device parameters.

Control: a dephasing channel gives per-round TVD exactly 0 at every round.

| code | leak order | partial sum at 300 rounds | attack error bound |
|---|---|---|---|
| repetition | 1 | 5.42 nats | 0.0022 |
| [[4,1,2]] | 2 | 0.058 nats | 0.472 |
| Hamming [[7,1,2]] | 3 | 3.2e-07 nats | 0.500 |

For repetition the per-round leak goes from 1.82e-2 to 1.54e-2 across 1500 rounds. The re-preparation
case reaches 1 percent at 251 rounds.

The cumulative exponent does not converge on this horizon: the increment over the first 100 rounds is
1.83 nats and over the last 100 rounds of 1500 it is 1.56 nats. Every figure here is a partial sum at
the stated round count. The `rounds` section carries the same quantity to 3000 rounds on five codes and
shows which of them stop growing.

---

## audit

Two checks on the machinery.

The attack against its closed form, one shot, 200000 trials:

| case | theory (1 - TVD)/2 | simulated |
|---|---|---|
| repetition, amplitude damping | 0.26000 | 0.26056 |
| [[4,1,2]], amplitude damping | 0.44880 | 0.45066 |
| [[4,1,2]], coherent diagonal | 0.34059 | 0.34012 |
| Hamming, amplitude damping | 0.49130 | 0.49080 |

The Chernoff exponent is a minimisation over a grid. Varying the grid from 51 to 2001 points moves the
value from 0.00100229 to 0.00100234. The default is 201 points.

---

## kl_blocks

The Knill-Laflamme matrix of each code under amplitude damping, computed two ways. The branch matrices are
P_C A_e^dagger A_e P_C for every Kraus product A_e, as in Leung, Nielsen, Chuang and Yamamoto (1997); their
largest eigenvalue gap over branches is the branch spread. The syndrome blocks resolve the same matrix by
syndrome, B_s[i, j] = Tr(Pi_s E(|j_L><i_L|)), and half the sum of |B_s[0,0] - B_s[1,1]| is the population
TVD. Orders are log-log slopes over gamma 0.005, 0.01 and 0.02.

Control: the no-jump branch of the [[4,1,2]] code at gamma 0.1 against equation 39 of Leung et al.

| source | eigenvalues of P_C A^dagger A P_C |
|---|---|
| this package | 0.810000000000, 0.828050000000 |
| Leung et al. eq. 39 | 0.810000000000, 0.828050000000 |

| code | analytic order | syndrome slope | branch spread slope |
|---|---|---|---|
| repetition | 1 | 0.99 | 0.99 |
| [[4,1,2]] | 2 | 1.98 | 1.99 |
| Hamming [[7,1,2]] | 3 | 2.97 | 2.98 |
| Steane [[7,1,3]] | none | none | 2.98 |
| five-qubit [[5,1,3]] | none | none | 5.00 |

On the three leaking codes the three columns agree. On Steane and the five-qubit code the branch spread is
nonzero, at orders 3 and 5, and the syndrome blocks have equal diagonals. The branch spread orders equal the
amplitude-damping population distances in `leak_order`.

---

## surface

The rotated [[9,1,3]] surface code on `expectations.py`, which computes the syndrome distribution from
the expectations of the stabilizer group elements, P(s) = 2^-m sum_a (-1)^(s.a) <S^a>, on a state vector
and the adjoint channel (Kobori and Todo 2025, eq. 15; Blume-Kohout and Young 2025, eqs. 7 and 34). It
holds no syndrome projectors, so the surface code runs where `Code` cannot be built.

Control 1: on the five standard codes, three channels and three logical states, the two engines agree to
5.6e-16.

Control 2: at amplitude damping 0.2, the density matrix is evolved and the eight projector factors are
applied one at a time to it, nothing stored. Difference from the expectation engine 1.7e-16.


Analytic leak order under amplitude damping: 3. Z_L is ZZZ on the first row and X_L is X on the first
column.

| channel | population TVD |
|---|---|
| amplitude damping 0.1 | 5.317650e-03 |
| amplitude damping 0.2 | 3.128320e-02 |
| coherent diagonal 0.3 | 0.0e+00 |
| depolarizing 0.1 | 0.0e+00 |

Every Z check has even weight and Z_L has odd weight, the condition of Hu, Liang and Calderbank (example
1) under which no transversal Z rotation moves the syndrome distribution. Under depolarizing the
stabilizer expectations themselves are the same for every logical state, to 2.2e-16 (Wagner et al.
2022, eq. 40).

---

## protection

Logical error after one round of ideal syndrome measurement and recovery, from `protection.py`, beside the
population leak. The recovery for syndrome s is the first single-qubit Pauli with that syndrome, or the
identity when none has it. For a pure logical state psi, R_s psi lies in the syndrome-s space, so the
fidelity after recovery is sum_s <R_s psi| E(rho) |R_s psi>, with no syndrome projectors. The distance d is
the smallest weight of a Pauli that commutes with every generator and anticommutes with Z_L or X_L.

Control 1: on the five standard codes, three channels and three logical states, the fidelity form and the
projector form sum_s R_s Pi_s E(rho) Pi_s R_s agree to 1.0e-15.

Control 2: on the four distance-3 codes, every single-qubit Pauli on three logical states is corrected to
4.4e-16.

Amplitude damping (AD) and depolarizing at strength 0.01. Slopes are of log value against log strength,
from 0.01 to 0.02. The logical error is the largest over |0_L>, |1_L> and |+_L>.

| code | d | leak AD | slope | logical AD | slope | logical depolarizing | slope |
|---|---|---|---|---|---|---|---|
| 3-qubit repetition | 1 | 2.970e-02 | 0.99 | 7.481e-03 | 1.00 | 1.973e-02 | 0.98 |
| [[4,1,2]] | 2 | 3.920e-04 | 1.97 | 1.000e-02 | 1.00 | 1.333e-02 | 1.00 |
| Hamming [[7,1,2]] | 2 | 3.891e-06 | 2.96 | 7.795e-03 | 1.05 | 1.122e-02 | 1.14 |
| 5-qubit [[5,1,3]] | 3 | 0.000e+00 | - | 1.856e-04 | 1.99 | 6.520e-04 | 1.97 |
| Steane [[7,1,3]] | 3 | 0.000e+00 | - | 5.172e-04 | 1.98 | 1.804e-03 | 1.95 |
| Shor [[9,1,3]] | 3 | 2.620e-05 | 2.96 | 1.095e-03 | 1.96 | 2.677e-03 | 1.94 |
| rotated surface [[9,1,3]] | 3 | 6.819e-06 | 2.96 | 8.795e-04 | 1.97 | 2.931e-03 | 1.93 |

The recovery above is the built-in table, the first single-qubit Pauli with each syndrome. Any decoder can
take its place. BP+OSD of `ldpc`, at depolarizing 0.01:

| code | d | table | BP+OSD | ratio |
|---|---|---|---|---|
| 5-qubit [[5,1,3]] | 3 | 6.520e-04 | - | not CSS |
| Steane [[7,1,3]] | 3 | 1.804e-03 | 9.047e-04 | 0.501 |
| Shor [[9,1,3]] | 3 | 2.677e-03 | 1.153e-03 | 0.431 |
| rotated surface [[9,1,3]] | 3 | 2.931e-03 | 7.795e-04 | 0.266 |

At distance 3 the table corrects every single-qubit error, so the factor of two to four is the weight-2
errors that the table has no entry for. The adapter splits a CSS code, decoding X errors from the Z checks
and Z errors from the X checks, so a non-CSS code has no adapter.

Past distance 3 the table stops being a decoder at all. Random errors on the d=5 rotated surface code, 200
of each weight, scored on syndromes alone so the 25 qubits never have to be simulated:

| weight | decoder | corrected | logical | unexplained |
|---|---|---|---|---|
| 1 | table | 200 | 0 | 0 |
| 1 | BP+OSD | 200 | 0 | 0 |
| 2 | table | 4 | 0 | 196 |
| 2 | BP+OSD | 200 | 0 | 0 |

Unexplained means the correction leaves a syndrome behind, which is what the table does when it holds no
entry for what was measured. Every logical-error number above this line uses the table, so at distance 3 it
is an upper bound on what a real decoder would leave.

The distance-3 codes have logical-error slope 2 under both channels. Of these, Steane and the five-qubit
code show no leak, and Shor and the rotated surface code leak with slope 3. The distance-2 and distance-1
codes have logical-error slope 1; for the built-in Hamming code the slope-1 term is its weight-2 logical X,
which moves |0_L> and |1_L> and leaves |+_L>.

---

## decoded

The leak that remains when an observer sees a decoder's output in place of the syndrome record, from
`decoded.py`. A view is a function of the syndrome: the correction a decoder returns, the qubits it acts on
(support), their number (weight), whether it acts at all (acted), and the two frame bits, which say whether
the correction anticommutes with Z_L (Z frame) or with X_L (X frame). The frame bit is the output Higgott
and Gidney describe for a matching decoder, a prediction of which logical observable measurements were
flipped. The distance between the two logical states' distributions of a view is at most the syndrome leak,
since the statistical distance does not increase under a function (Hastad et al., Proposition 2.2; Shen and
Zhong state it for a decoder), and equals it when the difference of the two syndrome distributions has one
sign on every set of syndromes the view merges. States |0_L> and |1_L>, one round, ideal syndrome measurement.

Three decoders. The least weight decoder returns, for each syndrome, the first Pauli of least weight that
has it. The coset decoder is the one of `rounds.py`. The single-qubit table is the default of
`protection.py` and returns the identity for a syndrome that no single-qubit Pauli has.

Control 1: under depolarizing noise of 0.1, on seven codes and three decoders, the largest value of any view
is 0.0e+00.

Control 2: under amplitude damping no view exceeds the syndrome leak; the largest excess is 2.8e-17.

Control 3: the coset and least weight decoders return a correction with the measured syndrome on every
code, and a different correction for every syndrome. The single-qubit table does not.

Share of the syndrome leak in each view of the least weight correction, amplitude damping 0.05:

| code | leak | correction | support | weight | acted | Z frame | X frame |
|---|---|---|---|---|---|---|---|
| 3-qubit repetition | 1.425e-01 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.3333 | 0.0000 |
| [[4,1,2]] | 9.025e-03 | 1.0000 | 0.7500 | 0.7500 | 0.2500 | 1.0000 | 0.0000 |
| Hamming [[7,1,2]] | 4.343e-04 | 1.0000 | 0.4075 | 0.2532 | 0.0633 | 0.0617 | 0.0925 |
| Shor [[9,1,3]] | 2.894e-03 | 1.0000 | 0.6528 | 0.4375 | 0.0625 | 0.0370 | 0.0625 |
| rotated surface [[9,1,3]] | 7.654e-04 | 1.0000 | 0.4660 | 0.2371 | 0.1051 | 0.1541 | 0.0705 |
| 5-qubit [[5,1,3]] | 0.000e+00 | - | - | - | - | - | - |
| Steane [[7,1,3]] | 0.000e+00 | - | - | - | - | - | - |

At amplitude damping 0.2 the shares of the Shor code are the same to four decimals, and those of the surface
code are 0.4730, 0.1980, 0.0990, 0.1849 and 0.0683.

The correction under the three decoders, amplitude damping 0.05:

| code | decoder | corrections | syndromes | correction | Z frame |
|---|---|---|---|---|---|
| Hamming [[7,1,2]] | single-qubit table | 19 | 64 | 0.4399 | 0.1558 |
| Hamming [[7,1,2]] | coset | 64 | 64 | 1.0000 | 0.0000 |
| Hamming [[7,1,2]] | least weight | 64 | 64 | 1.0000 | 0.0617 |
| Shor [[9,1,3]] | single-qubit table | 22 | 256 | 0.3750 | 0.1250 |
| Shor [[9,1,3]] | coset | 256 | 256 | 1.0000 | 0.2963 |
| Shor [[9,1,3]] | least weight | 256 | 256 | 1.0000 | 0.0370 |
| rotated surface [[9,1,3]] | single-qubit table | 24 | 256 | 0.5522 | 0.1724 |
| rotated surface [[9,1,3]] | coset | 256 | 256 | 1.0000 | 0.1474 |
| rotated surface [[9,1,3]] | least weight | 256 | 256 | 1.0000 | 0.1541 |

The frame bit depends on which representative of Z_L it is taken against. Multiplying Z_L by the stabilizer
with generator mask a adds a.s to the bit, which is equation 1 of Higgott and Gidney. Share of the syndrome
leak in the Z frame bit over the representatives, least weight decoder, amplitude damping 0.05; "Z type"
counts the representatives written with Z alone:

| code | stored | representatives | least | largest | Z type | least | largest | largest at |
|---|---|---|---|---|---|---|---|---|
| 3-qubit repetition | 0.3333 | 4 | 0.3333 | 1.0000 | 4 | 0.3333 | 1.0000 | ZZZ |
| [[4,1,2]] | 1.0000 | 8 | 0.0000 | 1.0000 | 4 | 0.0000 | 1.0000 | ZIZI |
| Hamming [[7,1,2]] | 0.0617 | 64 | 0.0000 | 0.5584 | 8 | 0.0000 | 0.4318 | ZZIZIIZ |
| Shor [[9,1,3]] | 0.0370 | 256 | 0.0000 | 1.0000 | 64 | 0.0370 | 1.0000 | ZZZZZZZZZ |
| rotated surface [[9,1,3]] | 0.1541 | 256 | 0.0000 | 0.9537 | 16 | 0.1383 | 0.9537 | ZZZZZZZZZ |

The frame bit of Z on every qubit is the parity of the number of X corrections. On the Shor code it carries
the whole leak, and the leak has a closed form: with p_b = 1 - prod(1 - gamma_q) - prod(gamma_q) over the
three qubits of block b, the chance that the Z checks of an excited block fire, the leak is p_1 p_2 p_3.
With one strength, p_b is eq. 2 of Shen and Zhong for the three-qubit repetition code at one round.
The engine, the frame bit and the closed form agree to 1.5e-16 at amplitude damping 0.2, 0.05, 0.01 and
0.005 and at strengths from 0.02 to 0.10 across the qubits.

Order in the damping strength, slope of log value against log strength from 0.005 to 0.01:

| code | syndrome | correction | support | weight | acted | Z frame | X frame |
|---|---|---|---|---|---|---|---|
| 3-qubit repetition | 0.99 | 0.99 | 0.99 | 0.99 | 0.99 | 0.99 | - |
| [[4,1,2]] | 1.99 | 1.99 | 1.99 | 1.99 | 1.99 | 1.99 | - |
| Hamming [[7,1,2]] | 2.98 | 2.98 | 2.98 | 2.98 | 2.98 | 2.98 | 2.98 |
| Shor [[9,1,3]] | 2.98 | 2.98 | 2.98 | 2.98 | 2.98 | 2.98 | 2.98 |
| rotated surface [[9,1,3]] | 2.98 | 2.98 | 2.98 | 2.97 | 2.98 | 2.99 | 2.98 |

No view raises the order.

The records of the `hardware_shor_pinned` run on `ibm_fez` through the same views. The record is the six Z
generators, and the least weight decoder returns 64 corrections for its 64 values, X on at most one qubit
of each block. Distance is the plug-in reading; achieved is the likelihood ratio fitted on half of each
record and scored on the other half.

| view | outcomes | distance, 20 us | 50 us | 100 us | achieved, 20 us | 50 us | 100 us |
|---|---|---|---|---|---|---|---|
| the record | 64 | 0.0438 | 0.1673 | 0.2337 | 0.0263 | 0.1655 | 0.2329 |
| correction | 64 | 0.0438 | 0.1673 | 0.2337 | 0.0271 | 0.1663 | 0.2336 |
| number of X corrections | 4 | 0.0367 | 0.1671 | 0.2337 | 0.0369 | 0.1671 | 0.2351 |
| their parity | 2 | 0.0367 | 0.1671 | 0.2337 | 0.0357 | 0.1668 | 0.2331 |
| Z frame, stored representative | 2 | 0.0094 | 0.0084 | 0.0033 | 0.0083 | 0.0061 | 0.0000 |
| acted | 2 | 0.0184 | 0.0407 | 0.0569 | 0.0193 | 0.0394 | 0.0566 |

At zero delay every view reads a permutation p of 0.4930 or more. At 50 and 100 us the one bit holds the
distance of the 64 outcomes to 0.0002. At 20 us the bit achieves 0.0357 where the 64-outcome record achieves
0.0263: the plug-in floor of the record is 0.0238 against 0.0044 for the bit, and the fit on 64 outcomes
loses more to its own noise than the bit loses by merging. Under amplitude damping over the delay with each
data qubit's calibrated T1 and nothing else, the number of X corrections and its parity carry 1.0000 of the
Z-record leak, the stored representative 0.0398, 0.0348 and 0.0271 at the three delays, and acted 0.2500;
the measured acted share of the achieved distance is 0.238 and 0.243 at 50 and 100 us.

The Z record at distance 5, under the reset process of `zchecks.py`. The characteristic function of the
syndrome, E[(-1)^(a.s)], is a mean over the stabilizer codewords of a product over their qubits, and its
Walsh-Hadamard transform is the distribution. The cost is the number of codewords times the number of Z
syndromes, 4096 by 4096 at distance 5. Control: on the Shor code and the distance 3 surface code, at five
reset settings, the distribution equals the enumeration of the process to 4.1e-15. At distance 5 the slope of
log leak against log reset from 0.01 to 0.02 is 4.88, for a leak order of 5.

Shares of the Z-record leak at the resets of the `zchecks` run. Fired is the number of Z checks that fire,
weight the least weight of an X error with the syndrome, parity the parity of that weight, which is the frame
bit of Z on every qubit; the last three columns are the Z frame bit of the stored representative and its
range over the representatives made of Z generators.

| code | reset | leak | fired | weight | parity | stored | least | largest |
|---|---|---|---|---|---|---|---|---|
| Shor [[9,1,3]] | 0.1153 | 2.866e-02 | 0.5926 | 1.0000 | 1.0000 | 0.0370 | 0.0370 | 1.0000 |
| Shor [[9,1,3]] | 0.3075 | 2.607e-01 | 0.5926 | 1.0000 | 1.0000 | 0.0370 | 0.0370 | 1.0000 |
| Shor [[9,1,3]] | 0.5204 | 4.198e-01 | 0.5926 | 1.0000 | 1.0000 | 0.0370 | 0.0370 | 1.0000 |
| rotated surface [[9,1,3]] | 0.1153 | 7.430e-03 | 0.7143 | 1.0000 | 1.0000 | 0.1429 | 0.1429 | 1.0000 |
| rotated surface [[9,1,3]] | 0.3075 | 6.759e-02 | 0.7143 | 1.0000 | 1.0000 | 0.1429 | 0.1429 | 1.0000 |
| rotated surface [[9,1,3]] | 0.5204 | 1.088e-01 | 0.7143 | 1.0000 | 1.0000 | 0.1429 | 0.1429 | 1.0000 |
| rotated surface d=5 | 0.1153 | 4.044e-04 | 0.2436 | 0.2323 | 0.1987 | 0.0395 | 0.0001 | 0.3289 |
| rotated surface d=5 | 0.3075 | 1.326e-02 | 0.1778 | 0.2720 | 0.0759 | 0.0585 | 0.0001 | 0.2832 |
| rotated surface d=5 | 0.5204 | 2.920e-02 | 0.1621 | 0.2967 | 0.0472 | 0.0626 | 0.0001 | 0.2727 |

On both distance 3 codes the parity of the least weight carries the whole Z-record leak. At distance 5 it
carries 0.0472 to 0.1987, the count of fired checks 0.1621 to 0.2436, and no frame bit more than 0.3289. The
`zchecks` run reads the distance 5 code through the count of fired checks and achieves 0.00249 and 0.00501
at resets 0.3075 and 0.5204; the count's exact value there is 0.00236 and 0.00473, and the Z-record leak is
1.326e-02 and 2.920e-02.

Limits:

- One round with ideal syndrome measurement in the exact engine, |0_L> against |1_L>, amplitude damping,
  codes of at most nine qubits and one logical qubit. The distance 5 rows are the Z record under the reset
  process, which the `zchecks` section checks against the exact engine on the smaller codes.
- The frame bit depends on the decoder's choice among corrections of equal weight and on the representative.
  The tables give one decoder's choice and the range over representatives, not a worst case over decoders.
- A decoder in a repeated experiment works on detection events across rounds. That case is not covered.
- The hardware rows reuse the records of one job: one code, one backend, the Z generators only.

---

## regions

The leak of a part of the record, from `regions.py`: the generators an observer sees, the fewest generators
whose joint record leaks, and the leak with one qubit made noiseless. The record of a set of generators is a
function of the whole record, so its distance is at most the whole record's; Shen and Zhong state this for a
coarse-grained log. A part therefore gives a lower bound, and a part that reads zero does not show that a code
is protected. States |0_L> and |1_L>, one round, ideal syndrome measurement, amplitude damping.

Control 1: the record of a set of generators, computed from the products of those generators alone, equals the
whole record summed over the other bits on every set of generators of four codes, to 1e-12 (test suite).

Control 2: under relaxation a set of Z generators leaks exactly when some product of them has a support that
holds a Z logical. The leak comes from the classical process of `zchecks`, the cover from the integer program:

| code | sets of Z generators | leak above 1e-12 | product covers a logical | disagree |
|---|---|---|---|---|
| Shor [[9,1,3]] | 63 | 27 | 27 | 0 |
| rotated surface [[9,1,3]] | 15 | 8 | 8 | 0 |

Sets of generators by size, exact engine, amplitude damping 0.1. The Shor code, whole record 0.019683, its Z
generators 0.019683, its X generators 0.000000:

| generators in the set | sets | that leak | largest leak |
|---|---|---|---|
| 1 | 8 | 0 | 0.000000 |
| 2 | 28 | 0 | 0.000000 |
| 3 | 56 | 8 | 0.005832 |
| 4 | 70 | 28 | 0.008748 |
| 5 | 56 | 38 | 0.013122 |
| 6 | 28 | 25 | 0.019683 |
| 7 | 8 | 8 | 0.019683 |
| 8 | 1 | 1 | 0.019683 |

The eight sets of three are one Z generator from each block, each with 0.005832. The rotated surface code,
whole record 0.005318, its Z generators 0.005103, its X generators 0.000000:

| generators in the set | sets | that leak | largest leak |
|---|---|---|---|
| 1 | 8 | 0 | 0.000000 |
| 2 | 28 | 3 | 0.002916 |
| 3 | 56 | 16 | 0.002916 |
| 4 | 70 | 35 | 0.005103 |
| 5 | 56 | 40 | 0.005103 |
| 6 | 28 | 25 | 0.005184 |
| 7 | 8 | 8 | 0.005233 |
| 8 | 1 | 1 | 0.005318 |

The three pairs are Z generators: (4, 6) and (5, 7) with 0.002916, and (4, 5) with 0.001458. The X generators
carry nothing alone and add to the Z generators' 0.005103 when seen with them.

One qubit made noiseless, the others at 0.1:

| code | whole record | least leak without one qubit | largest | least drop | largest drop | sum of the drops |
|---|---|---|---|---|---|---|
| Shor [[9,1,3]] | 0.019683 | 0.013851 | 0.013851 | 0.005832 | 0.005832 | 0.052488 |
| rotated surface [[9,1,3]] | 0.005318 | 0.003219 | 0.004071 | 0.001247 | 0.002098 | 0.012297 |

The largest drop of the surface code is the centre qubit's. Noise on fewer qubits than a Z logical holds gives
no leak, so a drop is not a share: the drops sum to more than the leak.

Past the exact engine, the fewest Z generators that leak, by integer programming, and the leak of that set from
the classical process restricted to the qubits it touches, damping 0.05. d/l* is the distance over the largest
generator weight. Shen and Zhong bound the region a fault must occupy to carry an input-dependent contribution
by d/l* cells, with l* the number of data qubits one detecting cell covers; that region is where the fault
acts, and the count here is of the generators an observer has to read:

| code | n | Z generators | d/l* | fewest | proven | qubits touched | leak of the set | whole Z record |
|---|---|---|---|---|---|---|---|---|
| rotated surface d=3 | 9 | 4 | 0.75 | 2 | True | 6 | 4.2869e-04 | 7.5020e-04 |
| rotated surface d=5 | 25 | 12 | 1.25 | 3 | True | 10 | 1.9345e-06 | 1.0657e-05 |
| rotated surface d=7 | 49 | 24 | 1.75 | 4 | True | 15 | 4.3646e-09 | - |
| rotated surface d=9 | 81 | 40 | 2.25 | 5 | True | 18 | 3.9391e-11 | - |
| rotated surface d=11 | 121 | 60 | 2.75 | 6 | True | 22 | 1.7773e-13 | - |

The fewest are (d + 1)/2 at each distance. The leak of the set is a lower bound on the leak of the whole record.
The `tensor` section gives another, the L2 distance at the same damping, 2.6684e-07 to 1.0384e-16 over distances
5 to 11; the set's leak is 7, 17, 217 and 1712 times that.

The distance 5 code, damping 0.05: a region grown one Z generator at a time from the best set of three, each
time by the generator that adds most. The whole Z record is 1.0657e-05:

| generators | leak | of the whole |
|---|---|---|
| 3 | 1.9345e-06 | 0.1815 |
| 4 | 2.3678e-06 | 0.2222 |
| 5 | 3.0534e-06 | 0.2865 |
| 6 | 5.6941e-06 | 0.5343 |
| 7 | 6.0878e-06 | 0.5713 |
| 8 | 6.6546e-06 | 0.6244 |
| 9 | 7.4278e-06 | 0.6970 |
| 10 | 7.6412e-06 | 0.7170 |
| 11 | 8.2701e-06 | 0.7760 |
| 12 | 1.0657e-05 | 1.0000 |

The same code with one qubit made noiseless: the drop as a fraction of the whole Z record, qubit 5 r + c at row
r and column c:

| row | column 0 | column 1 | column 2 | column 3 | column 4 |
|---|---|---|---|---|---|
| 0 | 0.0877 | 0.0877 | 0.0905 | 0.0905 | 0.1573 |
| 1 | 0.2193 | 0.1929 | 0.1858 | 0.2191 | 0.1573 |
| 2 | 0.2193 | 0.2281 | 0.2563 | 0.2281 | 0.2193 |
| 3 | 0.1573 | 0.2191 | 0.1858 | 0.1929 | 0.2193 |
| 4 | 0.1573 | 0.0905 | 0.0905 | 0.0877 | 0.0877 |

The least is 0.0877, the largest 0.2563 at the centre, and the 25 fractions sum to 4.1272.

X and Z generators together past the exact engine, from `regions.css_region_leak`: the distance 5 code,
damping 0.05. A product of generators is X on a set B and Z on a set A; the X and Y factors shrink by the
coherence factor and each Z outside B becomes (1 - p) Z + p I, so the record needs the Z stabilizers and Z
logicals that lie on the touched qubits and off B. The last row is the whole record of the code, 24 generators
on 25 qubits:

| Z generators | X generators | outcomes | leak | over the Z record |
|---|---|---|---|---|
| 12 | 0 | 4096 | 1.0657e-05 | 1.0000 |
| 12 | 4 | 65536 | 1.1472e-05 | 1.0764 |
| 12 | 8 | 1048576 | 1.2306e-05 | 1.1547 |
| 12 | 12 | 16777216 | 1.3193e-05 | 1.2379 |

Control: on the distance 3 code the same engine gives 7.654328e-04 for the whole record and the exact engine
7.654328e-04. The X generators add 0.2379 to the Z record at distance 5; at distance 3 the Z record is 0.005103
of 0.005318 at damping 0.1.

X and Y shrink by exp(-t / T2), so a set that holds X generators depends on T2. Decay probability 0.05, the
twelve Z and eight X generators:

| T2 over T1 | coherence | leak |
|---|---|---|
| 2.0 | 0.9747 | 1.2306e-05 |
| 1.0 | 0.9500 | 1.1332e-05 |
| 0.5 | 0.9025 | 1.0919e-05 |

The leak over damping^d as the damping falls, d the distance, which tends to the coefficient of the leading
order:

| code | record | 1e-2 | 1e-3 | 1e-4 |
|---|---|---|---|---|
| rotated surface d=3 | Z | 6.7921 | 6.9790 | 6.9979 |
| rotated surface d=3 | whole | 6.8189 | 6.9818 | 6.9982 |
| rotated surface d=5 | Z | 47.8156 | 51.5662 | 51.9565 |
| rotated surface d=5 | whole | 50.0272 | 51.8029 | 51.9803 |

Shen and Zhong give the number of Z logicals of least weight, 8 at distance 3 and 52 at distance 5, state that
this number is not the coefficient, and leave the coefficient at distance 5 undetermined. The values tend to 7
at distance 3 and to 52 at distance 5, for the Z record and for the whole record: the order is attained at
distance 5, and the X generators do not change the leading coefficient. Between the two distances the leak at
damping 0.01 falls by a factor of 1363.

A drop is signed. The least and the largest drop as a fraction of the leak, and the number of qubits whose
drop is negative, by damping strength:

| code | damping | leak | least | largest | negative |
|---|---|---|---|---|---|
| Shor [[9,1,3]] | 0.05 | 2.8936e-03 | 0.3158 | 0.3158 | 0 |
| Shor [[9,1,3]] | 0.20 | 1.1059e-01 | 0.2500 | 0.2500 | 0 |
| Shor [[9,1,3]] | 0.35 | 3.1791e-01 | 0.1538 | 0.1538 | 0 |
| Shor [[9,1,3]] | 0.45 | 4.0934e-01 | 0.0606 | 0.0606 | 0 |
| Shor [[9,1,3]] | 0.50 | 4.2188e-01 | 0.0000 | 0.0000 | 0 |
| Shor [[9,1,3]] | 0.60 | 3.7325e-01 | -0.1667 | -0.1667 | 9 |
| rotated surface [[9,1,3]] | 0.05 | 7.6543e-04 | 0.2612 | 0.4123 | 0 |
| rotated surface [[9,1,3]] | 0.20 | 3.1283e-02 | 0.1727 | 0.3540 | 0 |
| rotated surface [[9,1,3]] | 0.35 | 9.7638e-02 | 0.0499 | 0.2759 | 0 |
| rotated surface [[9,1,3]] | 0.45 | 1.3466e-01 | -0.0686 | 0.2076 | 8 |
| rotated surface [[9,1,3]] | 0.50 | 1.4453e-01 | -0.1461 | 0.1668 | 8 |
| rotated surface [[9,1,3]] | 0.60 | 1.4550e-01 | -0.3218 | 0.0935 | 8 |
| rotated surface d=5, Z record | 0.05 | 1.0657e-05 | 0.0877 | 0.2563 | 0 |
| rotated surface d=5, Z record | 0.20 | 3.3526e-03 | 0.0437 | 0.1701 | 0 |
| rotated surface d=5, Z record | 0.45 | 2.7972e-02 | 0.0180 | 0.0499 | 0 |

On the Shor code the leak is F^3 with F = 1 - (1 - g)^3 - g^3, and F^2 (1 - (1 - g)^2) with one qubit
noiseless; the two are equal at damping 0.5. An idle of 100 us at a T1 of 100 to 250 us is a damping of
0.33 to 0.63, so the range where a drop is negative is inside the hardware runs.

Coverage:

- Amplitude damping, one round, ideal syndrome measurement, the two logical basis states.
- The leak of a part is a lower bound on the leak of the whole record. At distance 5 the best three generators
  carry 0.1815 of the Z record and the last generator added takes it from 0.7760 to 1.0000.
- Sets that mix X and Z generators past the exact engine are read under relaxation only, a decay probability
  and a coherence factor per qubit. The whole record at distance 5 is 2^24 outcomes; distance 7 has 48
  generators and is out of reach as a whole.
- The limits 7 and 52 are read off three dampings each; they are not derived.
- The growth by one generator is one path. It is not the best set of each size.
- The drops rank qubits at weak damping. Past a damping of about 0.4 a noiseless qubit raises the leak for
  some qubits of the surface code and, past 0.5, for every qubit of the Shor code.
- Surface codes up to distance 11. The integer program proves each count; at distance 11 it takes 108 s.

---

## shor_circuits

The Shor-code circuits of `hardware.build_shor_circuits` on the local Aer simulator: nine data qubits and
eight ancillas. The encoder takes |+> or |-> on data qubit 0 to |0_L> or |1_L> of `hardware.shor_code()`,
whose Z_L is ZIIZIIZII; the input |0> gives |+_L>. The data idle for the delay, and one round measures the
two X generators and the six Z pairs, generator j into classical bit j. The syndrome-extraction layout is
that of the Shor-code circuits in QECops (Mithra 2026).

Control 1: encoder output against the |0_L> and |1_L> state vectors, max 1 - |overlap| 3.3e-16.

Control 2: 56 of 56 circuits, the error-free circuit and each of the 27 single-qubit Paulis for both
logical states, give the predicted syndrome on all 200 shots.

Amplitude damping moves population, which the Z generators measure. The two X generators carry no
population leak, so a round of the six Z generators alone keeps the leak and drops the circuit from 17
qubits and 32 two-qubit gates to 15 and 20. `build_shor_circuits(checks="z")` builds that round.

Population leak at the `ibm_fez` T1 of the `hardware` section, 133.19 us, against the null floor, the mean
distance between two samples of the |0_L> distribution at that delay. The three statistics are readings of
one syndrome record: all eight generators, the six Z generators, and the number of generators that fired.

| statistic | delay (us) | leak | floor, 4000 shots | floor, 16000 shots |
|---|---|---|---|---|
| full syndrome (256) | 20 | 0.0466 | 0.0580 | 0.0293 |
| full syndrome (256) | 50 | 0.2684 | 0.0860 | 0.0431 |
| full syndrome (256) | 100 | 0.4179 | 0.0963 | 0.0482 |
| Z generators (64) | 20 | 0.0466 | 0.0336 | 0.0173 |
| Z generators (64) | 50 | 0.2684 | 0.0451 | 0.0228 |
| Z generators (64) | 100 | 0.4179 | 0.0479 | 0.0240 |
| number that fired (7) | 20 | 0.0276 | 0.0140 | 0.0071 |
| number that fired (7) | 50 | 0.1591 | 0.0175 | 0.0086 |
| number that fired (7) | 100 | 0.2477 | 0.0168 | 0.0086 |

The floor grows with the number of outcomes that carry weight. On the full syndrome at 4000 shots the 20 us
leak of 0.0466 sits under the floor of 0.0580. Dropping the X generators puts it above. For comparison, the
four outcomes of the repetition code in the `hardware` section have a floor of 0.0052 at the same shot
count.

The ratio of leak to floor ranks the count of generators that fired first, at 2.0, 9.1 and 14.7. That
ranking does not hold up as a test: the `estimator` run measures the power of each statistic directly and
the six Z generators win.

---

## readout

Error on the syndrome bits themselves. Every other run measures the bits without error. A flip on bit j
with probability q multiplies the expectation of every group element that contains generator j by 1 - 2q,
so the engine applies it exactly, in one multiplication. Control against flips applied outcome by outcome:
maximum difference 5.6e-17.

Leak after readout error, over the leak without it:

| code | bits | q=0.005 | q=0.01 | q=0.02 | q=0.05 |
|---|---|---|---|---|---|
| 3-qubit repetition | 2 | 0.9867 | 0.9735 | 0.9472 | 0.8700 |
| Shor, six Z generators | 6 | 0.9606 | 0.9225 | 0.8498 | 0.6585 |

The factor is the same at 20, 50 and 100 us. Flipping bits scales the coefficients of both syndrome
distributions by the same amount whatever the damping is, so the attenuation does not move with the delay.
That is the form the `hardware` fit assumed when it fitted one constant. More bits cost more: six generators
lose 15 percent at q = 0.02 where two lose 5.

The fit needed 0.945, which on two bits is q = 0.0208. The readout error of `FakeFez` has median 0.0076,
quartiles 0.0049 and 0.0135, so readout alone covers 36 percent of the flip rate the fit implies.

The `shor_noise` run measured 1.09, 0.89 and 0.77 of the ideal leak at the three delays. That moves with the
delay, so bit flips alone do not account for it. Extraction also damps the data qubits, which raises the
leak where it is small and compresses it where it is large:

| model | 20 us | 50 us | 100 us |
|---|---|---|---|
| bit flips only, q=0.02 | 0.850 | 0.850 | 0.850 |
| extra damping 0.01, q=0.01 | 1.071 | 0.955 | 0.919 |
| extra damping 0.02, q=0.02 | 1.134 | 0.910 | 0.844 |
| extra damping 0.03, q=0.02 | 1.291 | 0.939 | 0.840 |
| `shor_noise`, full noise model | 1.09 | 0.89 | 0.77 |

Two terms, one flat in delay and one that is not, cover the shape and most of the size. The remainder is
what the gate-level model carries and this one does not: depolarizing on the data qubits, the ancillas
idling through the round, and the routing the transpiler added.

---

## disorder

One relaxation channel per qubit, taken from the T1 of the first nine qubits of `FakeFez`: 49, 256, 274,
220, 191, 190, 214, 140 and 75 us, against a median of 144.9 us over all 156. Every other run in this
document puts one channel on every qubit, and the `hardware` section had to use the device median because
the layout of the measured job was not recorded.

The matched column is the uniform channel at the mean of the nine gammas, so the disordered column differs
from it only by the spread.

| code | delay (us) | uniform at median | matched mean | disordered | over matched |
|---|---|---|---|---|---|
| Shor [[9,1,3]] | 20 | 0.03827 | 0.04395 | 0.04357 | 0.991 |
| Shor [[9,1,3]] | 50 | 0.23843 | 0.24172 | 0.26119 | 1.081 |
| Shor [[9,1,3]] | 100 | 0.42187 | 0.41951 | 0.47655 | 1.136 |
| rotated surface [[9,1,3]] | 20 | 0.01047 | 0.01207 | 0.00894 | 0.741 |
| rotated surface [[9,1,3]] | 50 | 0.07076 | 0.07183 | 0.06198 | 0.863 |
| rotated surface [[9,1,3]] | 100 | 0.14435 | 0.14113 | 0.13539 | 0.959 |

The spread moves the leak in opposite directions on the two codes: up to 1.136 of the matched value on
Shor at 100 us, down to 0.741 on the surface code at 20 us.

Placement is the larger effect. Dealing the same nine T1 values to the nine qubits in 200 random orders,
Shor at 50 us:

| quantity | value |
|---|---|
| smallest | 0.18136 |
| median | 0.24247 |
| largest | 0.27521 |
| spread | 0.09385 |
| spread over the median | 38.7 percent |

The uniform channel at the mean gamma reads 0.24172, above 48 percent of the orders. A single number for
the leak at a given median T1 is therefore a statement about an average layout, and the assignment of
physical qubits to code qubits moves it by a third of its own size.

The code measured in the `hardware` section has three data qubits. Drawing three T1 values at random from
the 156 of the backend, 300 times, against the median-T1 model it was compared with:

| delay (us) | median model | min | median | max | spread over median |
|---|---|---|---|---|---|
| 20 | 0.3370 | 0.2484 | 0.3537 | 0.6530 | 1.14 |
| 50 | 0.6201 | 0.4802 | 0.6401 | 0.8620 | 0.60 |
| 100 | 0.7500 | 0.5599 | 0.7507 | 0.8296 | 0.36 |

The median-T1 model sits near the median of the draws at every delay, and the draws span more than the
model's own value at 20 us. For this code the order does not matter: the same three values in all six
orders give one leak, 0.7350 at 50 us. Which three qubits carry the code is what moves the number.

The logical error after recovery changes little: 1.022 of the matched value on Shor and 1.003 on the
surface code, at 50 us. The leak is the quantity sensitive to the spread.

---

## rounds

A held state through many rounds of extraction and recovery. A decoder that explains every syndrome
returns each branch to the code space, so a round is one channel on the logical qubit: a 4 by 4 transfer
matrix on (trace, x, y, z) and an affine map from those coordinates to the syndrome distribution. Building
that map is the whole cost, and the round count is free after it, which is what puts the nine-qubit codes
in reach. `eavesdrop.repeated_extraction` holds every projector and stops around seven qubits.

Control: against that density-matrix path with the same recovery, the per-round distance agrees to 1.1e-16
on the repetition code and 1.5e-16 on the five-qubit code, over six rounds. A dephasing channel on the
repetition code gives exactly 0 at every round.

Held at the `ibm_marrakesh` medians with a 1 us cycle, recovered by BP+OSD, which explains every syndrome
on all five codes. The sum is of the per-round Chernoff exponents, in nats:

| code | order | round 1 | round 300 | sum 10 | sum 100 | sum 300 | sum 1000 | sum 3000 |
|---|---|---|---|---|---|---|---|---|
| 3-qubit repetition | 1 | 1.820e-02 | 1.761e-02 | 1.836e-01 | 1.827e+00 | 5.420e+00 | 1.738e+01 | 4.679e+01 |
| [[4,1,2]] | 2 | 1.473e-04 | 3.782e-06 | 3.188e-04 | 9.525e-04 | 1.019e-03 | 1.020e-03 | 1.020e-03 |
| Hamming [[7,1,2]] | 3 | 8.936e-07 | 5.894e-10 | 3.884e-09 | 9.952e-09 | 1.003e-08 | 1.003e-08 | 1.003e-08 |
| Shor [[9,1,3]] | 3 | 6.032e-06 | 5.458e-06 | 9.566e-06 | 7.339e-05 | 1.808e-04 | 4.138e-04 | 6.576e-04 |
| rotated surface [[9,1,3]] | 3 | 1.564e-06 | 1.401e-06 | 1.095e-08 | 1.060e-07 | 2.959e-07 | 7.780e-07 | 1.329e-06 |

Two behaviours, and they follow the distance rather than the leak order. The two distance-2 codes lose the
logical state to the noise, the two held states converge, and the sum stops: 1.020e-03 and 1.003e-08. The
codes that correct well keep the state distinguishable, so the per-round leak barely decays and the sum
grows without bound: the repetition code from 1.820e-02 to 1.761e-02 per round, the surface code from
1.564e-06 to 1.401e-06.

A held state is therefore not protected by saturation. It is protected by the rate, which spans four
orders of magnitude across these codes. Read as attack bounds, 0.5 exp(-sum) at 3000 rounds: the repetition code
falls below 1e-20, while Shor sits at 0.4997 and the surface code at 0.499999.

The sum also depends on the decoder, because a correction that damages the state ends the leak with it:

| decoder | round 1 | round 300 | sum 300 | sum 3000 |
|---|---|---|---|---|
| single-qubit table | 1.820e-02 | 1.761e-02 | 5.420e+00 | 4.679e+01 |
| BP+OSD | 1.820e-02 | 1.761e-02 | 5.420e+00 | 4.679e+01 |
| any coset representative | 1.820e-02 | 2.917e-03 | 2.520e+00 | 2.995e+00 |

On the repetition code the table and BP+OSD agree, since both correct single bit flips. A decoder that
returns any Pauli carrying the measured syndrome corrects the syndrome and not the error, the state decays,
and the accumulated leak saturates at 2.995 instead of growing. Every accumulated number belongs to a
stated decoder.

---

## estimator

The finite-sample tools of `estimate.py` against known leaks. A record is a set of shots, so every reading
carries sampling noise: the plug-in distance between two sampled distributions sits above the floor of its
shot count even when the two states are identical, which is why the run reports a permutation test for any
difference and a distance that a fitted rule achieves on held-out shots.

Power is the fraction of draws whose permutation p falls below 0.05 at the stated leak, over 60 draws. The
false-positive column runs the same test on two draws of one state and must sit near 0.05. The corrected
column is the achieved distance with its own null taken off.

| delay (us) | shots | statistic | leak | power | false positive | corrected |
|---|---|---|---|---|---|---|
| 20 | 4000 | full syndrome (256) | 0.0466 | 0.62 | 0.05 | 0.0138 |
| 20 | 4000 | Z generators (64) | 0.0466 | 0.78 | 0.05 | 0.0195 |
| 20 | 4000 | number that fired (7) | 0.0276 | 0.68 | 0.08 | 0.0161 |
| 20 | 16000 | full syndrome (256) | 0.0466 | 1.00 | 0.03 | 0.0278 |
| 20 | 16000 | Z generators (64) | 0.0466 | 1.00 | 0.05 | 0.0369 |
| 20 | 16000 | number that fired (7) | 0.0276 | 1.00 | 0.08 | 0.0227 |
| 50 | 4000 | full syndrome (256) | 0.2684 | 1.00 | 0.00 | 0.2230 |
| 50 | 4000 | Z generators (64) | 0.2684 | 1.00 | 0.08 | 0.2603 |
| 50 | 4000 | number that fired (7) | 0.1591 | 1.00 | 0.03 | 0.1518 |
| 50 | 16000 | full syndrome (256) | 0.2684 | 1.00 | 0.02 | 0.2635 |
| 50 | 16000 | Z generators (64) | 0.2684 | 1.00 | 0.03 | 0.2664 |
| 50 | 16000 | number that fired (7) | 0.1591 | 1.00 | 0.05 | 0.1565 |

Every statistic holds its false-positive rate near the nominal 0.05. At the 20 us leak and 4000 shots the
six Z generators reach power 0.78, against 0.68 for the count that fired and 0.62 for the full syndrome. At
16000 shots all three reach 1.00, and the Z generators give the largest corrected distance, 0.0369 against
a leak of 0.0466.

The corrected distance sits below the leak when the shots are few, since the rule is fitted on half of
them. At 50 us and 16000 shots it reads 0.2664 against 0.2684.

---

## shor_noise

The Shor Z-check circuits under the calibration of `FakeFez`, the rehearsal for a device run. The
transpiler places and routes them on the 156-qubit device at optimization level 2, giving depth 64 and 52
cz gates on 15 qubits. The circuit is then rebuilt on those 15 qubits, each carrying its own T1, T2, gate
errors and readout error.

The null column is two runs of |0_L> at the same delay, which is the floor a device run measures when
there is no leak. Aer draws its shots from `seed_simulator` and close seeds share random numbers, so every
run here has a seed spaced by more than the shot count.

| shots | delay (us) | leak | null | ideal | measured / ideal | fired leak | fired null |
|---|---|---|---|---|---|---|---|
| 4000 | 0 | 0.0198 | 0.0248 | 0 | - | 0.0038 | 0.0145 |
| 4000 | 20 | 0.0588 | 0.0422 | 0.0466 | 1.260 | 0.0290 | 0.0152 |
| 4000 | 50 | 0.2528 | 0.0513 | 0.2684 | 0.942 | 0.1270 | 0.0273 |
| 4000 | 100 | 0.3205 | 0.0447 | 0.4179 | 0.767 | 0.1922 | 0.0155 |
| 16000 | 0 | 0.0082 | 0.0076 | 0 | - | 0.0026 | 0.0018 |
| 16000 | 20 | 0.0509 | 0.0241 | 0.0466 | 1.091 | 0.0296 | 0.0112 |
| 16000 | 50 | 0.2396 | 0.0258 | 0.2684 | 0.893 | 0.1394 | 0.0053 |
| 16000 | 100 | 0.3222 | 0.0268 | 0.4179 | 0.771 | 0.1839 | 0.0076 |

At zero delay the leak sits at the null at both shot counts, so this pair of circuits carries no
preparation offset. |1_L> differs from |0_L> by a Z on one data qubit, which the transpiler turns into a
frame change, so both preparations run the same gates. The repetition code of the `hardware` section
prepares |1_L> with an X on every data qubit and measures an offset of 0.0140 at zero delay.

Gate and readout error hold the measured leak below the ideal channel at the long delays, 0.77 of it at
100 us, and they put 17 percent of the zero-delay shots on a nonzero syndrome.

At 4000 shots the 20 us leak of 0.0588 stands against a null of 0.0422. At 16000 it is 0.0509 against
0.0241, and every delay separates from its null.

---

## hardware_shor

Measured on IBM `ibm_fez`, job `dao2md5r85ps73fen5lg`. The Shor [[9,1,3]] code with one round of its six Z
generators onto six ancillas, prepared in |0_L> and |1_L>, at four idle delays, 16000 shots per circuit,
41 quantum seconds on the open plan. Raw records in `results/hardware_shor_ibm_fez.json`, analysis in
`results/hardware_shor_analysis.txt` from `analysis_shor_hardware.py`.

The layout is recorded with the job, which the `hardware` run did not do: data qubits 111, 91, 93, 107,
118, 110, 90, 78, 68 with T1 from 63 to 182 us, ancillas 98, 92, 108, 109, 89, 69 with readout error from
0.0045 to 0.0146. Those are the qubits of the zero-delay pair. The transpiler placed the eight circuits on
three layouts, listed per circuit under `layouts` in the job file: the zero-delay pair on one, the three
nonzero-delay |0_L> circuits on a second and the three nonzero-delay |1_L> circuits on a third, so the
two states of a nonzero-delay pair ran on different qubits. Every circuit transpiled to depth 63 or 64
and 52 cz, the same for both logical states, since |1_L> differs from |0_L> by a Z on one data qubit,
which becomes a frame change.

The records through `estimate.py`, over the 64 outcomes of the six generators:

| delay (us) | distance | floor | p | achieved | corrected | 95% interval |
|---|---|---|---|---|---|---|
| 0 | 0.0108 | 0.0133 | 0.8385 | 0.0000 | 0.0000 | (0.0000, 0.0133) |
| 20 | 0.0988 | 0.0256 | 0.0000 | 0.0884 | 0.0860 | (0.0730, 0.1039) |
| 50 | 0.2233 | 0.0292 | 0.0000 | 0.2219 | 0.2202 | (0.2068, 0.2370) |
| 100 | 0.2049 | 0.0286 | 0.0000 | 0.2006 | 0.1986 | (0.1854, 0.2158) |

At zero delay the distance sits under its floor with p = 0.8385, so this pair of circuits carries no
preparation offset. The repetition code of the `hardware` section measured 0.0140 there with p = 0.0065.
At every nonzero delay the permutation test gives p = 0.0000 and the achieved distance is 0.0884 or more:
the records of the two circuits differ. At those delays the two circuits also ran on different qubits, so
the distance holds a layout term beside the state term; the `shor_device_model` section separates the two
under a model.

The same records read as the number of generators that fired, 7 outcomes:

| delay (us) | distance | floor | p | achieved | corrected |
|---|---|---|---|---|---|
| 0 | 0.0026 | 0.0057 | 0.8860 | 0.0000 | 0.0000 |
| 20 | 0.0554 | 0.0087 | 0.0000 | 0.0539 | 0.0523 |
| 50 | 0.1254 | 0.0093 | 0.0000 | 0.1218 | 0.1203 |
| 100 | 0.1084 | 0.0091 | 0.0000 | 0.1069 | 0.1046 |

Against models of the same round. Uniform is one damping channel at the median T1 of the nine data qubits;
per qubit gives each data qubit its own T1 and T2 in the layout of the zero-delay pair; readout adds a flip on each
syndrome bit at its ancilla's readout error; exposure runs the readout model over the delay plus the
3.8 us the scheduled circuit itself takes.

| delay (us) | measured | uniform | per qubit | readout | exposure | measured / exposure |
|---|---|---|---|---|---|---|
| 0 | 0.0108 | 0 | 0 | 0 | 0.0006 | under the floor |
| 20 | 0.0988 | 0.0411 | 0.0501 | 0.0466 | 0.0686 | 1.440 |
| 50 | 0.2233 | 0.2492 | 0.2745 | 0.2554 | 0.2792 | 0.800 |
| 100 | 0.2049 | 0.4214 | 0.4071 | 0.3788 | 0.3739 | 0.548 |

None of the models reproduces the magnitudes. The measured leak is above every model at 20 us and below
every model at 100 us, where it is also smaller than at 50 us. Counting the circuit's own 3.8 us, which
matters because the leak goes as the cube of the damping strength, moves the 20 us ratio from 2.122 to
1.440 and predicts 0.0006 at zero delay, under the floor of 0.0133, as measured. It does not move the long
delays.

The direction of both deviations is the one the `shor_noise` run gave under a gate-level noise model,
1.09, 0.89 and 0.77, and the device is further out in both. The `readout` section names the mechanism that
the analytic models lack: error on the data qubits during the 52 cz of the round. It acts as extra damping,
which raises the leak where it is small, and it scrambles both distributions toward uniform, which lowers
the distance where it is large.

Trivial-syndrome probability, measured:

| delay (us) | \|0_L> | \|1_L> |
|---|---|---|
| 0 | 0.7143 | 0.7158 |
| 20 | 0.3047 | 0.3039 |
| 50 | 0.2311 | 0.1729 |
| 100 | 0.2449 | 0.1761 |

At zero delay 29 percent of shots land on a nonzero syndrome from gate and readout error alone, against
17 percent under the `FakeFez` calibration of the `shor_noise` run.

The run establishes that the syndrome records of the two circuits differ at every nonzero delay on a
second code, of distance 3, with no preparation offset at zero delay. Because those pairs ran on
different qubits, the measured distance at a nonzero delay combines a state term and a layout term,
which are separated under a model in the next section and not on the device. The run does not establish
any of the present models as a quantitative prediction of a device. The `hardware_shor_pinned` section
repeats the round with every circuit on one layout, which separates the two on the device.

---

## shor_device_model

The eight circuits of job `dao2md5r85ps73fen5lg` as `ibm_fez` executed them, 52 cz each, sampled at 16000
shots under a model in which every gate, every scheduled idle and every readout carries its own qubit's
calibration (`gate_level.py`), against what the device measured. Two calibrations: the one read on
2026-09-22 (`results/ibm_fez_calibration.json`, last updated 10:00:30 that day, "today"), and the same
file with the T1, T2 and ancilla readout numbers the job recorded on 2026-09-20 laid over the qubits of
the first circuit ("job day"). Output in `results/shor_device_model.txt`, counts in
`results/hardware_shor_model_today.json` and `results/hardware_shor_model_job_day.json`.

The transpiler placed the circuits on three layouts. The zero-delay pair shares one, data qubits 111, 91,
93, 107, 118, 110, 90, 78, 68; the three nonzero-delay |0_L> circuits share a second, data qubits 91, 111,
113, 87, 78, 90, 110, 108, 97; the three nonzero-delay |1_L> circuits share a third, data qubits 109, 111,
113, 105, 97, 108, 98, 90, 88. At every nonzero delay the measured distance therefore compares two layouts
as well as two states. The model separates them with a twin, the |0_L> circuit rebuilt on the |1_L>
circuit's qubits (`results/hardware_shor_ibm_fez_twins.pkl`, 54 cz against 52): `state` is the distance
between the twin and the |1_L> circuit, two states on one layout, and `layout` is the distance between the
|0_L> circuit and its twin, one state on two layouts.

| delay (us) | measured | model, today | ratio | state | layout | model, job day | ratio | state | layout |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.0108 | 0.0097 | - | - | - | 0.0113 | - | - | - |
| 20 | 0.0988 | 0.0613 | 1.612 | 0.0666 | 0.0648 | 0.0472 | 2.094 | 0.0597 | 0.0530 |
| 50 | 0.2233 | 0.2403 | 0.929 | 0.2315 | 0.0637 | 0.1921 | 1.162 | 0.1916 | 0.0651 |
| 100 | 0.2049 | 0.3164 | 0.647 | 0.3163 | 0.0571 | 0.3235 | 0.633 | 0.3213 | 0.0336 |

Under either calibration the layout term stays between 0.03 and 0.07 at every delay, while the state term
grows with delay: at 20 us the two are the same size, and at 50 and 100 us the state term is 0.19 to 0.32. Neither calibration reproduces the measured magnitudes. The model sits below the
measurement at 20 us and above it at 100 us, and it rises with delay where the measurement falls between
50 and 100 us. Moving from today's calibration to the job-day numbers changes the 20 us model from 0.0613
to 0.0472 and leaves 100 us near 0.32; one data qubit of the first layout had T1 182 us on the job day
and 70 us today.

Trivial-syndrome probability of |0_L>, measured against the model:

| delay (us) | measured | today | job day |
|---|---|---|---|
| 0 | 0.7143 | 0.8184 | 0.8179 |
| 20 | 0.3047 | 0.4290 | 0.4669 |
| 50 | 0.2311 | 0.2734 | 0.2949 |
| 100 | 0.2449 | 0.2602 | 0.2545 |

At zero delay the model puts 18 percent of shots on a nonzero syndrome and the device 29 percent, so the
calibration's gate and readout numbers account for about two thirds of the error the round makes on its
own. Largest difference between model and measurement over the 64 outcomes and both states:

| delay (us) | today | job day |
|---|---|---|
| 0 | 0.1041 | 0.1036 |
| 20 | 0.1243 | 0.1622 |
| 50 | 0.0456 | 0.0759 |
| 100 | 0.0270 | 0.0326 |

The per-outcome agreement is worst at 0 and 20 us and within 0.03 at 100 us, where both distributions are
near their long-delay shape.

---

## hardware_shor_pinned

The Shor-code round of `hardware_shor` again on `ibm_fez`, job `daquif6ekp0c73arbd70`, submitted 2026-09-25,
with all eight circuits on one layout. The layout was chosen once, over 24 transpiler seeds, as data qubits
113, 114, 111, 112, 110, 98, 119, 134, 132 and ancillas 133, 115, 109, 91, 135, 131, and every circuit was
transpiled onto it. The eight circuits the job returned share their qubits, their 47 cz and their
measurement map. Dynamical decoupling and twirling were off. 16000 shots per circuit, 41 quantum seconds on
the open plan. The calibration was read three seconds before submission. The run is
`scripts/shor_pinned_ibm_fez.py` and the analysis `syndrome_leakage/analysis_shor_pinned.py`.

Control: zero delay, where |1_L> differs from |0_L> by a Z on one data qubit before encoding and the two
records must agree.

The records through `estimate.py`, six Z generators, 64 outcomes:

| delay (us) | distance | floor | p | achieved | 95% interval |
|---|---|---|---|---|---|
| 0 | 0.0117 | 0.0144 | 0.8580 | 0.0000 | 0.0000 to 0.0132 |
| 20 | 0.0438 | 0.0238 | 0.0000 | 0.0254 | 0.0121 to 0.0430 |
| 50 | 0.1673 | 0.0282 | 0.0000 | 0.1627 | 0.1489 to 0.1795 |
| 100 | 0.2337 | 0.0295 | 0.0000 | 0.2304 | 0.2171 to 0.2472 |

`achieved` is the distance an observer reaches on held-out shots, with its null subtracted. At zero delay
the records agree within the floor. At 20, 50 and 100 us they differ, with p below 0.0005 at each. The
number of generators that fired gives the same picture, 0.0294, 0.0980 and 0.1224 achieved at the three
nonzero delays.

Against the gate-level model of `shor_device_model`, built from the calibration read before submission:

| delay (us) | measured | model | measured / model | trivial, measured | trivial, model |
|---|---|---|---|---|---|
| 0 | 0.0117 | 0.0088 | | 0.7016 | 0.8394 |
| 20 | 0.0438 | 0.0478 | 0.916 | 0.3792 | 0.4517 |
| 50 | 0.1673 | 0.2483 | 0.674 | 0.2477 | 0.2879 |
| 100 | 0.2337 | 0.3446 | 0.678 | 0.2215 | 0.2594 |

The trivial columns are the probability of the all-zero syndrome for |0_L>. At zero delay the device puts
30 percent of shots on a nonzero syndrome and the model 16 percent, the same gap the 2026-09-20 job showed on
different qubits, 29 against 18. The model has the direction and the rise with delay, and its distance is
too large by a factor of 1.48 at 50 us and 1.47 at 100 us.

Beside the 2026-09-20 job, whose circuits the transpiler placed on three layouts:

| delay (us) | one layout | three layouts |
|---|---|---|
| 0 | 0.0117 | 0.0108 |
| 20 | 0.0438 | 0.0988 |
| 50 | 0.1673 | 0.2233 |
| 100 | 0.2337 | 0.2049 |

The two jobs ran on different qubits with different T1 and T2, so the columns are not a like-for-like
comparison of magnitude. The one-layout distance rises at every step in delay, where the three-layout run
fell between 50 and 100 us.

The run establishes that the Shor code's syndrome record depends on the encoded logical state on this
device: with both states on the same qubits and the same gates, the records differ at every nonzero delay
and agree at zero delay. The two circuits of a pair ran one after the other in one job, so drift within the
job is the one difference between them left uncontrolled. The error the round makes beyond its calibration
is the same fraction on a second layout. One code, one backend, one job.

---

## drop_sign

The sign of the drop of `regions.qubit_drops`: the leak minus the leak with one qubit made noiseless. Each
term of the leak under relaxation carries the damping of every qubit of one Z logical and the survival
1 - gamma of the other qubits of a stabilizer, so a noiseless qubit removes the terms that need its damping and
raises the others, and the sign depends on which wins. Shen and Zhong describe one side of this for unequal
rates: at large exposure the fastest qubits "have already decayed and stop contributing". Exact engine,
amplitude damping, one round, ideal measurement.

Weak damping, 0.001 on every qubit. A qubit that lies on a Z logical of least weight inside the support of a
product of Z generators has a positive drop; a qubit that lies on none has a negative one:

| code | order | on a logical | off | drop > 0 | drop < 0 | disagree | least drop over the leak |
|---|---|---|---|---|---|---|---|
| 3-qubit repetition | 1 | 3 | 0 | 3 | 0 | 0 | 0.3330 |
| [[4,1,2]] | 2 | 4 | 0 | 4 | 0 | 0 | 0.4995 |
| Hamming [[7,1,2]] | 3 | 6 | 1 | 6 | 1 | 0 | -0.0001 |
| Shor [[9,1,3]] | 3 | 9 | 0 | 9 | 0 | 0 | 0.3330 |
| rotated surface [[9,1,3]] | 3 | 9 | 0 | 9 | 0 | 0 | 0.2852 |

Equal damping on every qubit: the damping at which a qubit's drop first changes sign, and that damping as an
idle time over T1:

| code | qubits | damping | idle / T1 |
|---|---|---|---|
| 3-qubit repetition | all three | 0.5000 | 0.6931 |
| [[4,1,2]] | all four | 0.5000 | 0.6931 |
| Hamming [[7,1,2]] | all but qubit 3 | 0.5193 | 0.7325 |
| Shor [[9,1,3]] | all nine | 0.5000 | 0.6931 |
| rotated surface [[9,1,3]] | the four edge qubits | 0.3968 | 0.5055 |
| rotated surface [[9,1,3]] | the four corner qubits | 0.4323 | 0.5662 |
| rotated surface [[9,1,3]] | the centre qubit | 0.7122 | 1.2455 |

Qubit 3 of the Hamming code has a negative drop at every damping.

The Shor code in closed form. Its leak is the product over the three blocks of F = 1 - prod(1 - gamma) -
prod(gamma), and the drop of qubit q is gamma_q (1 - gamma_a - gamma_b) times the F of the other two blocks,
a and b the other two qubits of its block: the sign is set by the two block mates alone. Over 40 draws of nine
rates between 0.02 and 0.95 the engine differs from that form by at most 3.3e-16.

The Z record of the distance 5 surface code, equal damping:

| damping | leak | least drop over the leak | largest | negative |
|---|---|---|---|---|
| 0.45 | 2.7972e-02 | 0.0180 | 0.0499 | 0 |
| 0.60 | 2.3928e-02 | -0.1751 | -0.0720 | 25 |
| 0.75 | 7.1678e-03 | -0.7900 | -0.3151 | 25 |
| 0.90 | 2.2440e-04 | -3.0175 | -1.1267 | 25 |

The Shor code at the T1 of the data qubits of job `daquif6ekp0c73arbd70`, 67 to 166 us, amplitude damping over
the idle alone:

| delay us | damping from | to | leak | least drop over the leak | largest | negative |
|---|---|---|---|---|---|---|
| 20 | 0.1134 | 0.2571 | 0.0715 | 0.1130 | 0.3490 | 0 |
| 50 | 0.2599 | 0.5243 | 0.3339 | -0.0114 | 0.2218 | 1 |
| 100 | 0.4523 | 0.7737 | 0.3902 | -0.3418 | -0.0016 | 9 |
| 200 | 0.7000 | 0.9488 | 0.0907 | -1.6687 | -0.7364 | 9 |
| 400 | 0.9100 | 0.9974 | 0.0016 | -9.4761 | -5.2988 | 9 |

At the job's 100 us delay every one of the nine drops is negative: in this model, making any one data qubit
noiseless raises the leak.

Coverage:

- Five codes of at most nine qubits on the exact engine, and the Z record of one code of 25 qubits.
- One code, the Hamming code, has a qubit off the cover of least-weight logicals; the rule at weak damping
  rests on that qubit and on the 31 qubits that are on a cover.
- Damping over the idle alone in the last table; the gates and the readout of the circuit are not in it.

---

## circuit_parts

The pinned Shor round of `hardware_shor_pinned` under the gate-level model built from the calibration read
before submission, with parts of the circuit made noiseless, from `gate_level.run_round(..., ideal=)`. The two
barriers of the circuit split it into the encoder, the window and the extraction; the readout is the fourth
part. At 100 us the encoder holds 66 gates and 23 idles, the window 15 idles and the extraction 128 gates and
31 idles, on 15 qubits. 8000 shots per circuit, seeds spaced by more than the shot count, the same seeds in
every setting. The sampled counts are in `results/circuit_parts_counts.json` and the tables are computed from
them.

Control: with every part noiseless, every circuit returns the trivial syndrome on every shot.

Distance between the two logical states' records with one part made noiseless. The floor is the mean distance
of two samples of one record at this shot count, for the setting with every part noisy:

| delay us | every part noisy | without encoder | without window | without extraction | without readout | floor |
|---|---|---|---|---|---|---|
| 0 | 0.0116 | 0.0140 | 0.0116 | 0.0061 | 0.0121 | 0.0121 |
| 20 | 0.0426 | 0.0436 | 0.0123 | 0.0499 | 0.0450 | 0.0296 |
| 50 | 0.2450 | 0.2596 | 0.0120 | 0.2704 | 0.2501 | 0.0381 |
| 100 | 0.3525 | 0.3680 | 0.0145 | 0.3999 | 0.3648 | 0.0404 |

With one part noisy alone:

| delay us | every part noisy | only encoder | only window | only extraction | only readout |
|---|---|---|---|---|---|
| 0 | 0.0116 | 0.0034 | 0.0000 | 0.0106 | 0.0023 |
| 20 | 0.0426 | 0.0060 | 0.0546 | 0.0086 | 0.0047 |
| 50 | 0.2450 | 0.0026 | 0.2989 | 0.0136 | 0.0070 |
| 100 | 0.3525 | 0.0050 | 0.4335 | 0.0087 | 0.0031 |

The window carries the distance. With the window noiseless the distance is 0.0120 to 0.0145 at every delay,
at its own floor of 0.0125 with a permutation p of 0.223 to 0.530. The window alone gives more than every part
together. The encoder, the extraction and the readout alone give 0.0023 to 0.0136, each at its own floor; of
their twelve permutation tests two have p under 0.05 (0.048 and 0.025), and a readout flip in this model does
not depend on the state.

The distance of a setting minus the distance with every part noisy, with the 2.5 to 97.5 percent interval of
that difference over 600 resamplings of both records:

| delay us | setting | difference | from | to |
|---|---|---|---|---|
| 50 | without encoder | 0.0146 | -0.0071 | 0.0345 |
| 50 | without window | -0.2330 | -0.2438 | -0.2117 |
| 50 | without extraction | 0.0254 | 0.0039 | 0.0461 |
| 50 | without readout | 0.0051 | -0.0155 | 0.0280 |
| 50 | only window | 0.0539 | 0.0331 | 0.0750 |
| 100 | without encoder | 0.0155 | -0.0043 | 0.0376 |
| 100 | without window | -0.3380 | -0.3504 | -0.3176 |
| 100 | without extraction | 0.0474 | 0.0271 | 0.0667 |
| 100 | without readout | 0.0122 | -0.0089 | 0.0314 |
| 100 | only window | 0.0810 | 0.0604 | 0.1003 |

The distance rises when the extraction is noiseless, by 0.0254 at 50 us and 0.0474 at 100 us, with intervals
that exclude zero. The rises with the encoder or the readout noiseless have intervals that include zero. So the
noise of the extraction lowers the distance the window produces: the difference with a part made noiseless is
signed, as for a qubit in the `regions` section. Shen and Zhong state the direction: additional
state-independent error "randomises both classes together and can reduce" distinguishability. They place the
exposure in the idle time of the syndrome round by its duration, about 2.5 percent of it in the two-qubit
layers; here the parts are switched off in the model one at a time.

The trivial syndrome's probability, every part noisy and the window noiseless:

| delay us | 0_L | 1_L | 0_L, window noiseless | 1_L, window noiseless |
|---|---|---|---|---|
| 0 | 0.8384 | 0.8439 | 0.8384 | 0.8439 |
| 20 | 0.4586 | 0.4551 | 0.8387 | 0.8376 |
| 50 | 0.2806 | 0.2179 | 0.8360 | 0.8365 |
| 100 | 0.2696 | 0.1726 | 0.8431 | 0.8369 |

Coverage:

- The model is split, not the device. The model's distance at 50 and 100 us is 0.2450 and 0.3525, where the
  device gave 0.1673 and 0.2337 in `hardware_shor_pinned`.
- 8000 shots per circuit. A distance at its floor is not resolved, so the encoder, the extraction and the
  readout alone are bounded by their floors, 0.0038 to 0.0095, and not shown to be zero.
- A part is a position in the circuit. A qubit that idles between the gates of the extraction counts as
  extraction, and relaxation during a gate counts with that gate.
- The resampled distances carry their sampling floor a second time, so the intervals are a guide to the sign of
  a difference and not a calibrated test.
- One code, one job's circuits, one calibration.

---

## circuit_parts_surface

The analysis of `circuit_parts` on a second code: the distance 3 rotated surface code, one round of its four Z
generators from `hardware.build_extraction_circuits`, 9 data and 4 ancilla qubits, transpiled for the `FakeFez`
backend of qiskit-ibm-runtime with every circuit on one layout (transpiler seed 15, the fewest two-qubit gates
over seeds 11 to 16, 61 of them), under the model built from that backend's calibration. The data qubits have T1
of 63 to 209 us. At 100 us the encoder holds 83 gates and 24 idles, the window 13 idles and the extraction 142
gates and 32 idles. 8000 shots per circuit; the counts are in `results/circuit_parts_surface_counts.json`.

Control: with every part noiseless, every circuit returns the trivial syndrome on every shot.

Distance between the two logical states' records with one part made noiseless, and the floor of the setting
with every part noisy:

| delay us | every part noisy | without encoder | without window | without extraction | without readout | floor |
|---|---|---|---|---|---|---|
| 0 | 0.0128 | 0.0074 | 0.0128 | 0.0119 | 0.0136 | 0.0120 |
| 20 | 0.0214 | 0.0257 | 0.0105 | 0.0160 | 0.0210 | 0.0195 |
| 50 | 0.0833 | 0.0896 | 0.0142 | 0.0949 | 0.0846 | 0.0221 |
| 100 | 0.1011 | 0.1047 | 0.0105 | 0.1104 | 0.1058 | 0.0224 |

With one part noisy alone:

| delay us | every part noisy | only encoder | only window | only extraction | only readout |
|---|---|---|---|---|---|
| 0 | 0.0128 | 0.0095 | 0.0000 | 0.0080 | 0.0024 |
| 20 | 0.0214 | 0.0066 | 0.0242 | 0.0110 | 0.0044 |
| 50 | 0.0833 | 0.0086 | 0.0996 | 0.0086 | 0.0041 |
| 100 | 0.1011 | 0.0117 | 0.1166 | 0.0085 | 0.0036 |

As on the Shor round, the window carries the distance: with it noiseless the distance is 0.0105 to 0.0142 at
every delay, at its own floor (permutation p 0.198 to 0.670), and the encoder, the extraction and the readout
alone are at their floors, with one of twelve tests at p 0.035. The distance at 20 us, 0.0214, is not resolved
(p 0.292).

The distance of a setting minus the distance with every part noisy, with the 2.5 to 97.5 percent interval over
600 resamplings:

| delay us | setting | difference | from | to |
|---|---|---|---|---|
| 50 | without window | -0.0690 | -0.0811 | -0.0494 |
| 50 | without extraction | 0.0116 | -0.0103 | 0.0305 |
| 50 | only window | 0.0164 | -0.0038 | 0.0354 |
| 100 | without window | -0.0906 | -0.1031 | -0.0706 |
| 100 | without extraction | 0.0092 | -0.0116 | 0.0296 |
| 100 | only window | 0.0155 | -0.0035 | 0.0361 |

The window alone reads above every part together and the distance rises with the extraction noiseless, as on
the Shor round, but here each interval includes zero: at a distance of 0.1 and 8000 shots the sign of those
differences is not resolved.

Coverage:

- A model on a fake backend's calibration; no device run of this code exists here.
- 8000 shots per circuit; the differences between settings other than the window's are inside the sampling
  noise.
- One round, Z generators only, one layout.

---

## hardware

Measured on IBM `ibm_fez`. A 3-qubit repetition code with one round of Z-stabilizer
extraction onto two ancillas, prepared in |0_L> and |1_L>, at four idle delays. Median T1 133.2 us, 4000
shots per circuit. Raw data in `results/hardware_ibm_fez.json`, job metadata in
`results/hardware_ibm_fez_job.json`.

Preparing |1_L> requires X gates that |0_L> does not, so preparation and readout asymmetry contribute to
the measured distance at every delay. The idle-delay sweep separates the two: an amplitude-damping leak
grows with delay and an asymmetry offset is constant. The zero-delay point measures the offset.

| delay (us) | gamma | syndrome TVD | Chernoff | 1-shot attack |
|---|---|---|---|---|
| 0 | 0 | 0.0140 | 0.0014 | 0.4930 |
| 20 | 0.1394 | 0.3412 | 0.1198 | 0.3296 |
| 50 | 0.3130 | 0.5888 | 0.3147 | 0.2055 |
| 100 | 0.5280 | 0.7110 | 0.4480 | 0.1435 |

The offset at zero delay is 0.0140 and the distance rises to 0.7110 at 100 us.

Measured against the simulated ideal channel at matched gamma:

| delay (us) | gamma | measured | simulated | difference |
|---|---|---|---|---|
| 0 | 0 | 0.0140 | 0 | +0.0140 |
| 20 | 0.1394 | 0.3412 | 0.3600 | -0.0187 |
| 50 | 0.3130 | 0.5888 | 0.6451 | -0.0563 |
| 100 | 0.5280 | 0.7110 | 0.7476 | -0.0366 |

Measured values are below simulated at every nonzero delay. The simulated model has no ancilla readout
error and no gate error.

Raw syndrome distributions over (00, 01, 10, 11):

| delay (us) | state | 00 | 01 | 10 | 11 |
|---|---|---|---|---|---|
| 0 | \|0_L> | 0.9505 | 0.0158 | 0.0262 | 0.0075 |
| 0 | \|1_L> | 0.9445 | 0.0297 | 0.0205 | 0.0053 |
| 20 | \|0_L> | 0.9560 | 0.0150 | 0.0262 | 0.0027 |
| 20 | \|1_L> | 0.6148 | 0.1510 | 0.1348 | 0.0995 |
| 50 | \|0_L> | 0.9623 | 0.0135 | 0.0230 | 0.0013 |
| 50 | \|1_L> | 0.3735 | 0.2263 | 0.2052 | 0.1950 |
| 100 | \|0_L> | 0.9600 | 0.0147 | 0.0235 | 0.0018 |
| 100 | \|1_L> | 0.2490 | 0.2495 | 0.2507 | 0.2507 |

The |0_L> distribution is flat across delays. The |1_L> distribution moves from 0.9445 on the trivial
syndrome at zero delay to approximately uniform at 100 us.

### Analysis

Output of `analysis_hardware.py`, in `results/hardware_analysis.txt`:

Sampling error, by resampling both distributions at the shot count:

| delay (us) | measured | bootstrap SE |
|---|---|---|
| 0 | 0.0140 | 0.0031 |
| 20 | 0.3412 | 0.0082 |
| 50 | 0.5888 | 0.0082 |
| 100 | 0.7110 | 0.0074 |

Null floor. Two independent 4000-shot samples drawn from one distribution give a mean distance of 0.0052,
a 95th percentile of 0.0100, and a maximum of 0.0190 over 4000 repetitions. The zero-delay distance of
0.0140 has p = 0.005.

The same records through `estimate.py`, which tests for any difference by permutation and measures what a
fitted likelihood ratio achieves on held-out shots:

| delay (us) | distance | floor | p | achieved | corrected | 95% interval |
|---|---|---|---|---|---|---|
| 0 | 0.0140 | 0.0054 | 0.0065 | 0.0130 | 0.0114 | (0.0000, 0.0440) |
| 20 | 0.3412 | 0.0105 | 0.0000 | 0.3394 | 0.3364 | (0.3102, 0.3685) |
| 50 | 0.5888 | 0.0126 | 0.0000 | 0.5880 | 0.5840 | (0.5629, 0.6131) |
| 100 | 0.7110 | 0.0134 | 0.0000 | 0.7120 | 0.7081 | (0.6902, 0.7337) |

The permutation test puts the zero-delay offset at p = 0.0065, beside the bootstrap's p = 0.005 above. At
100 us the achieved distance of 0.7120 is an attack error of 0.1440 on held-out shots, against the 0.1435
of the table above, which is the closed form at the measured distance.

T1 recovered from the leak. The simulated side of this fit puts one channel on every qubit at a single T1,
and the layout of this job was not recorded. The `disorder` run measures what that costs: three qubits
drawn at random from the backend give leaks spanning more than the model's own value at 20 us, though the
model sits near the median of those draws. The agreement below is therefore a statement about this job,
not a general accuracy of the single-T1 model.

Fitting `measured = a * simulated(gamma(t, T1))` over the three nonzero delays, with `a` a constant
attenuation:

| quantity | value |
|---|---|
| T1 fitted from the syndrome leak | 138.5 us |
| device median T1 from calibration | 133.2 us |
| ratio | 1.040 |
| attenuation a | 0.945 |
| rms residual | 0.0089 |

The fitted T1 is 4.0 percent above the device median. The rms residual is 0.0089 and the bootstrap
standard error is about 0.008. The per-delay ratios of measured to simulated are 0.948, 0.913 and 0.951.

Concerns and limits:

- The transpiler chose the qubits and the saved job record has no layout, so the comparison uses the
  device median T1 across 156 qubits. Qubit-to-qubit T1 variation on this device is wide compared with
  the 4 percent agreement.
- With no layout recorded, it is not established that the two logical states ran on the same qubits.
  The `shor_device_model` section shows the transpiler placing the two states of a pair on different
  qubits in the later job, where the difference between two layouts measured 0.03 to 0.07. That term
  can be present here and cannot be separated from these records.
- Four delay points and two fitted parameters leave two degrees of freedom. The fit is consistent with
  amplitude damping and has limited power to exclude an additional mechanism of similar shape.
- The zero-delay offset is 0.0140 at p = 0.005, against 0.7110 at 100 us.
- One code, one round of extraction, one backend. The measurement covers the amplitude-damping population
  leak. It leaves the coherent channel, the held-memory regime, and other codes untested on hardware.
