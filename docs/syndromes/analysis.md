# Analysis of a code

## Analyse one code under one channel

```python
import syndrome_leakage as sl

print(sl.analyze(sl.codes.steane(), "amplitude_damping", gamma=0.2))
```

The report gives the population leak, the phase leak, the analytic order, the measured slope of the exact
simulation, and a Pauli control that reads zero. Channel names come from `channels.LIBRARY`; a Kraus list
may be passed in place of a name.

## Measure the leak directly

```python
from syndrome_leakage.analyze import population_leak
from syndrome_leakage.channels import amplitude_damping

leak, d0, d1 = population_leak(sl.codes.STANDARD["repetition"](), amplitude_damping(0.2))
```

`leak` is the total variation distance between the two state-conditioned syndrome distributions; `d0` and
`d1` are the distributions themselves.

## Run the attack

```python
from syndrome_leakage.eavesdrop import chernoff_exponent, ml_attack, rounds_for_error

C, s_star = chernoff_exponent(d0, d1)
err = ml_attack(d0, d1, rounds=200, trials=200000)
n = rounds_for_error(d0, d1, target=0.01)
```

`ml_attack` simulates the likelihood-ratio test on independent shots and returns its error rate. The
Chernoff rate and the Bhattacharyya bound are reported beside it. The standard error is about
sqrt(e(1-e)/trials), so 4000 trials support two decimals and 200000 support three.

## Use device parameters

```python
from syndrome_leakage.channels import channel_from_t1t2

K = channel_from_t1t2(T1=163.3e-6, T2=77.0e-6, gate_time=1e-6)
```

T1 and T2 in seconds, gate time being the duration of one syndrome-extraction cycle.

The same numbers can be read from a Qiskit backend:

```python
from syndrome_leakage.hardware import channel_from_backend, syndrome_leak_from_backend

K = channel_from_backend(backend, qubit=0)                 # T1, T2 and the sx duration from the backend
r = syndrome_leak_from_backend(backend, code="steane")     # the attack at the worst data-qubit T1
```

`syndrome_leak_from_backend` takes the syndrome-extraction time as four two-qubit gate durations unless
`syndrome_time` is given, and returns the attack point with the T1, the time and the qubits used.

## Repeated extraction

```python
from syndrome_leakage.eavesdrop import repeated_extraction

res = repeated_extraction(code, kraus_override=K, rounds=300, correct=True)
```

Holds one logical state and extracts syndromes repeatedly, applying recovery each round when `correct` is
set. The ensemble is evolved exactly, with no trajectory sampling. Each round reports its own TVD and
Chernoff exponent; summing the exponents gives a partial sum at that round count, which is not a converged
total unless the per-round values have gone to zero.

## Knill-Laflamme matrices

```python
from syndrome_leakage.approx_qec import branch_matrices, branch_spread, syndrome_blocks, population_split

M = branch_matrices(code, amplitude_damping(0.1))      # {Kraus product e: <i_L| A_e^dagger A_e |j_L>}
gap, e = branch_spread(code, amplitude_damping(0.1))   # largest eigenvalue gap and its branch
B = syndrome_blocks(code, amplitude_damping(0.1))      # {syndrome s: Tr(Pi_s E(|j_L><i_L|))}
tvd = population_split(B)
```

`branch_matrices` builds every product of the single-qubit Kraus operators, so its cost is the number of
Kraus operators to the power n; it stops above 65536 products. The blocks from `syndrome_blocks` sum to the
2x2 identity for a trace-preserving channel.

## Closed form for a transversal Z-rotation

```python
from syndrome_leakage.hlc import generator_coefficients, syndrome_probabilities

A = generator_coefficients(code, theta=0.3)                # {(syndrome, g): A_mu,g}, eqs. 45 and 48
p = syndrome_probabilities(code, 0.3, z_expectation=1.0)   # eq. 91 for a state with <Z_L> = 1
```

This follows Hu, Liang and Calderbank (arXiv:2109.13481) for CSS codes with positive stabilizer signs. In
their equation 91 the logical state enters only through <Z_L>. The function enumerates all 2^n Z-strings.

## Logical error after recovery

```python
from syndrome_leakage import protection as pr

pr.protection_row(code, amplitude_damping(0.01))   # {"0_L": ..., "1_L": ..., "+_L": ...}
pr.logical_error(code, kraus, t, phi)              # one logical state
pr.distance(code)                                  # any stabilizer code, weight up to 4
pr.recovery_strings(code)                          # syndrome bits -> recovery Pauli string
```

One round of ideal syndrome measurement, then a correction. The default is the first single-qubit Pauli
with the measured syndrome, or the identity when none has it; `decoder=` takes any callable from syndrome
bits to a Pauli string, and `css_decoder(code)` wraps BP+OSD of `ldpc`, decoding X errors from the Z checks
and Z errors from the X checks. A branch whose correction carries a different syndrome contributes nothing,
since the projector kills it. `outcome(code, error, correction)` scores a correction on syndromes alone,
which works on codes of any size. The logical error is 1 - <psi| rho' |psi>, computed as
1 - sum_s <R_s psi| E(rho) |R_s psi>. It takes any `Code` or strings-only code up to about twelve qubits,
where the n-qubit density matrix fits in memory.

## The leak through a decoder's output

```python
from syndrome_leakage import decoded as dc

dec = dc.min_weight_decoder(code)                  # syndrome bits -> a Pauli of least weight
dc.decoded_leak(code, amplitude_damping(0.05))     # {"syndrome": ..., "correction": ..., "weight": ..., ...}
labels = dc.view_labels(code, dec)                 # one integer per syndrome for each view
dc.frame_spectrum(d0, d1, labels["frame_z"])       # the frame bit's leak for every representative of Z_L
```

An observer may see what a decoder returns instead of the syndrome. Each view is a function of the
syndrome, so its leak is at most the syndrome leak. The views are the correction, the qubits it acts on
(`support`), their number (`weight`), whether it acts (`acted`), and the frame bits `frame_z` and `frame_x`,
which say whether the correction anticommutes with Z_L or X_L. A decoder that returns a correction with
the measured syndrome gives a different correction for every syndrome, and the correction then carries the
whole leak. `view_labels(code, decoder, checks=...)` labels the records of a subset of the generators, and
`transform(labels)` passes a view to the functions of `estimate.py` for measured records. Entry `a` of
`frame_spectrum` is the leak of the frame bit taken against Z_L times the stabilizer with generator mask
`a`; `representative(code, code.zl_str, a)` returns that string and `type_masks(code, "Z")` the masks made
of Z generators alone. Any callable from syndrome bits to a Pauli string serves as the decoder, as in `protection.py`.

## Hold a state through many rounds

```python
from syndrome_leakage import rounds as rd

rows = rd.iterate(code, kraus, rounds=1000, decoder=rd.coset_decoder(code))
T, A = rd.round_map(code, kraus)      # the 4 by 4 logical map and the syndrome map
```

Each row carries the per-round distance, the Chernoff exponent and the running sum. A decoder that
explains every syndrome returns each branch to the code space, so a round is one channel on the logical
qubit: `round_map` returns its transfer matrix on (trace, x, y, z) and the affine map from those
coordinates to the syndrome distribution. The cost lies in building the map; iterating it does not grow
with the round count. `coset_decoder(code)` is complete by construction; the single-qubit table is complete only on
perfect codes, and `is_complete` says which. `eavesdrop.repeated_extraction` is the density-matrix path,
which holds every projector and stops around seven qubits.
