# Results

## Setup

### Synthetic

`synthetic_qrng` generates shots from a device model: a source bias, per-qubit readout applied as an
asymmetric assignment matrix, optional cross-qubit copying, and optional drift. Ground truth is known by
construction, so a credited rate can be compared against the entropy the source actually carries. The
estimators run on the resulting bit stream with no knowledge of the model.

### Hardware

IBM `ibm_marrakesh`, job `dablfc5e36ac739fjmmg`, qubits [0, 1], 2000 shots per circuit, submitted through
qiskit-ibm-runtime with the SamplerV2 primitive. Three circuits: prepare |0> and measure, prepare |1> and
measure, and prepare |+> and measure. The first two calibrate the readout on the same device and in the
same job as the third, which is the source being assessed.

---

## demo

`python -m qrng_attest demo --shots 4000 --qubits 3`, in `results/demo.txt`. A synthetic source with a
known model, showing the certificate format.

```
attested min-entropy = 0.9077 bits/bit   (verdict: PASS)
classical estimate: IID track (source passed the IID checks)
extractable = 10766 uniform bits at epsilon = 2^-64
```

---

## nist_vectors

`python scripts/nist_vectors.py`, in `results/nist_vectors.txt`. Known-answer check against the NIST SP
800-90B reference vectors. The reference column is the output NIST publishes for its `ea_non_iid` on the
same vectors, kept in `tests/nist_vectors/*.res`.

`rand1_short`, binary, per estimator:

| estimator | this package | NIST | difference |
|---|---|---|---|
| most_common_value | 0.9610588257 | 0.9610588257 | 1.1e-16 |
| collision | 0.6914641210 | 0.6914641210 | 4.8e-15 |
| markov | 0.9875961045 | 0.9875961045 | 1.1e-15 |
| compression | 0.6117162048 | 0.6117162048 | 1.1e-16 |
| t_tuple | 0.8676244308 | 0.8676244308 | 0 |
| lrs | 0.9626258038 | 0.9626258038 | 0 |
| multimcw | 0.9526180605 | 0.9526180605 | 0 |
| lag | 0.9433337066 | 0.9433337066 | 0 |
| multimmc | 0.9616166678 | 0.9616166678 | 1.1e-16 |
| lz78y | 0.9614462460 | 0.9614462460 | 0 |

`rand4_short` (alphabet 16) and `rand8_short` (alphabet 256) cover the multi-bit path, where a source is
assessed as literal symbols and as its bitstring. Every estimator of both tracks, `H_original`,
`H_bitstring` and the assessed value are compared, 20 values per vector:

| vector | H_original | H_bitstring | assessed | worst difference |
|---|---|---|---|---|
| rand4_short | 3.5674726724 | 0.8038720670 | 3.2154882679 | 2.3e-15 |
| rand8_short | 6.6364412871 | 0.7326117181 | 5.8608937445 | 2.3e-11 |

The 2.3e-11 is on the Lag predictor of `rand8_short`, whose estimate is set by the longest run of correct
predictions; the reference evaluates that bound in extended precision. For the four predictors the
reference also prints the number of correct predictions, the longest run and the number of predictions.
These integers are equal on all 20 predictor blocks of the three vectors.

`python scripts/nist_vectors.py <clone> [vector ...]`, in `results/nist_reference_vectors.txt`, runs the
whole assessment on the vectors of a clone of NIST's repository and compares every min-entropy line of
both tracks, `H_original`, `H_bitstring` and the assessed value with `cpp/selftest/refdata/*.res`. All
eleven vectors of the repository at commit 87c104d0, version 1.1.8:

| vector | samples | alphabet | values compared | assessed | NIST | worst difference |
|---|---|---|---|---|---|---|
| biased-random-bits | 1000000 | 2 | 12 | 0.0177665791 | 0.0177665791 | 8.22e-12 |
| biased-random-bytes | 1000000 | 256 | 20 | 0.2577408710 | 0.2577408710 | 1.86e-10 |
| data.pi | 1165666 | 2 | 12 | 0.5695375935 | 0.5695375935 | 4.33e-15 |
| normal | 1000000 | 179 | 20 | 4.1000957831 | 4.1000957831 | 6.39e-14 |
| rand1_short | 10000 | 2 | 12 | 0.6117162048 | 0.6117162048 | 4.77e-15 |
| rand4_short | 10000 | 16 | 20 | 3.2154882679 | 3.2154882679 | 2.33e-15 |
| rand8_short | 10000 | 256 | 20 | 5.8608937445 | 5.8608937445 | 2.27e-11 |
| ringOsc-nist | 1000000 | 2 | 12 | 0.1264457362 | 0.1264457362 | 3.21e-10 |
| truerand_1bit | 1000000 | 2 | 12 | 0.8296770832 | 0.8296770832 | 2.10e-14 |
| truerand_4bit | 1000000 | 16 | 20 | 3.6025063686 | 3.6025063686 | 3.47e-14 |
| truerand_8bit | 1000000 | 256 | 20 | 7.2338614552 | 7.2338614552 | 1.01e-13 |

180 values, worst difference 3.21e-10 bits. A vector of a million 8-bit samples takes about five minutes.
The reference is NIST's software at one commit; agreement with it is not a NIST validation.

---

## overcredit

`python -m qrng_attest.overcredit`, in `results/overcredit.txt`. A near-deterministic source, P(0) = 0.980
and true min-entropy 0.0291 bits/bit, read through a characterized symmetric readout of fidelity F, 60000
shots per point.

| F | P(0) observed | passes IID | SP 800-90B estimate | attested | estimate minus 0.0291 |
|---|---|---|---|---|---|
| 0.550 | 0.546 | yes | 0.8580 | 0.0000 | +0.8289 |
| 0.650 | 0.641 | yes | 0.6293 | 0.0170 | +0.6002 |
| 0.750 | 0.740 | yes | 0.4247 | 0.0145 | +0.3956 |
| 0.850 | 0.837 | yes | 0.2498 | 0.0187 | +0.2207 |
| 0.950 | 0.932 | yes | 0.0974 | 0.0247 | +0.0682 |
| 0.999 | 0.979 | yes | 0.0282 | 0.0268 | -0.0010 |

Every observed bit is an independent draw at bias `p_obs = F p + (1 - F)(1 - p)`, so the stream is
Bernoulli-IID and passes the IID screen at every point. The SP 800-90B estimate is 0.858 bits/bit at
F = 0.55. The attested value is at or below 0.0291 in every row.

The file also carries the crosstalk and mis-calibration experiments.

---

## hardware_ibm_marrakesh

Measured, in `results/hardware_ibm_marrakesh.txt` with raw values in the accompanying `.json`.

Self-calibration, from the |0> and |1> preparations in the same job:

| qubit | P(read 1 \| prep 0) | P(read 1 \| prep 1) | F0 | F1 |
|---|---|---|---|---|
| 0 | 0.0050 | 0.9945 | 0.9950 | 0.9945 |
| 1 | 0.0045 | 0.9900 | 0.9955 | 0.9900 |

The |+> circuit gives P(1) = 0.5115 and 0.5135 on the two qubits. The certificate:

```
attested min-entropy = 0.8607 bits/bit   (verdict: FLAGGED)
classical estimate: IID track (source passed the IID checks)
extractable = 3316 uniform bits at epsilon = 2^-64
action: q1: readout error sets the attested rate (0.86 bits/bit)
```

The verdict is FLAGGED because the readout error of qubit 1 sets the attested rate.

Scope: one backend, two qubits, 2000 shots, one calibration cycle.
