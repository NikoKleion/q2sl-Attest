# Hardware

## Extraction circuits for a code

```python
from syndrome_leakage.hardware import build_extraction_circuits, select_checks

circs, labels = build_extraction_circuits(code, [0.0, 20e-6], checks="z")
```

For each idle delay, one circuit prepares |0_L> and one |1_L>, idles, and measures one round of the
selected generators, one ancilla each, generator j into classical bit j. Z-type generators are measured by
CX into the ancilla and all others through an ancilla in |+>, in the order of `select_checks(code, checks)`.
Preparation is the standard CSS construction: the pivot qubits of Hx in |+>, CX onto the rest of each row,
then X_L for the second state.

| argument | meaning |
|---|---|
| `checks` | `"all"`, `"z"` or `"x"`: the generators measured |
| `prep_errors` | `{(delay, prep): (qubit, "X" or "Y" or "Z")}`, a Pauli after the delay, for checking the syndrome map |
| `flags` | one flag qubit per generator of weight 3 or more, read into the register `"flag"` |
| `faults` | `{check: (gates_done, pauli)}`, a Pauli on that check's ancilla, for testing the flags |
| `final_data` | measure every data qubit in Z at the end, into the register `"data"` |

The flag is coupled to the syndrome ancilla after the first data gate and before the last, the circuit of
Chao and Reichardt (2018), Fig. 2(c). Without a fault it never fires; a single ancilla fault that spreads to
two or more data qubits fires it. This is the distance-3 construction: at distance d, Chao and Reichardt
(2020) use d flags and the conditions of Chamberland and Beverland, which are not built here.
`split_registers(counts, widths)` separates the registers of a result key.

`build_circuits` builds the 3-qubit repetition code with one round of ZZI and IZZ onto two ancillas.
`build_shor_circuits` builds the Shor [[9,1,3]] code; its extraction layout derives from QECops by Jithesh
Mithra (MIT licensed, doi:10.5281/zenodo.19410365, see NOTICE). QECops encodes |0> and |1>, which are
|+_L> and |-_L> in this package's logical basis and carry no population leak, so `build_shor_circuits`
feeds |+> and |-> to the encoder to prepare |0_L> and |1_L>. `checks="z"` measures the six Z pairs only, on
15 qubits; the two X generators carry no population leak.

## Device model

```python
from syndrome_leakage import gate_level as gl

cal = gl.calibration_from_target(backend.target, qubits, pairs)   # or gl.load_calibration(path)
counts = gl.run_round(transpiled_circuits, cal, shots, gl.spaced_seeds(len(transpiled_circuits), shots))
quiet = gl.run_round(transpiled_circuits, cal, shots, seeds, ideal=("window",))   # the delay made noiseless
```

A transpiled round on a model in which every gate, idle and readout carries its own qubit's calibration.
`compress` rebuilds the circuit on the physical qubits it gates. `with_idle_noise` turns every scheduled
delay into a relaxation channel of its own duration on that qubit, so the idle during the delay and the
idles between gates both count. `noise_model` puts depolarizing error and relaxation on every one- and
two-qubit gate and a symmetric flip on every readout. `calibration_from_target` reads the numbers from a
`BackendV2` target, and `save_calibration` writes them with the date so that a run repeats without an
account.

The circuit builders put a barrier before and after the delay, which splits a round into the encoder, the
window and the extraction; the readout is a fourth part. `run_round(..., ideal=)` takes any of `gl.PARTS`
and runs the round with those parts noiseless: their gates carry no error and no relaxation, their idles
are left out, and an ideal readout drops the readout flips. `part_settings()` lists the ten settings of a
drop analysis: every part noisy, each part noiseless in turn, each part noisy alone, and none noisy, which
must return the trivial syndrome on every shot. A part is a position in the circuit, so a qubit that idles
between the gates of the extraction counts as extraction.

Aer draws shots from `seed_simulator`, and runs with close seeds share random numbers;
`spaced_seeds` spaces them by more than the shot count.

## Export a round to Stim

```python
from syndrome_leakage import stim_export as se

text = se.stim_text(code, prep=0, checks="all", gamma=0.05, readout=0.01, flags=True)
```

The text holds the encoder, an idle under the Pauli twirl of amplitude damping, the extraction round, a
final measurement of the data, one DETECTOR per ancilla and per flag, one more per Z-type generator against
the data, and Z_L as OBSERVABLE_INCLUDE(0). A detector is a set of measurements whose xor is deterministic
without noise (Gidney 2021, section 5.6), which holds because the encoder leaves every generator fixed.
Stim carries Pauli noise only, and under a Pauli channel the syndrome leak is zero, so the text serves as an
interchange format for decoders and samplers and as a control. `se.parse(text)` reads the text back
without Stim.

## A delay sweep on a backend

```python
from syndrome_leakage.hardware import build_shor_circuits, syndrome_dist_from_counts, analyse, report

circs, labels = build_shor_circuits([0.0, 20e-6, 50e-6, 100e-6], checks="z")
dist = syndrome_dist_from_counts(counts, shots, n_bits=6)
```

Transpile with `scheduling_method="alap"`, since the circuits carry explicit delays. Turn each result into
a distribution with `syndrome_dist_from_counts`, key the distributions by `(delay, prep)`, and pass them to
`analyse` and `report`, or to `estimate.leak_from_dists`. The zero-delay pair is the control: for the Shor
code |1_L> differs from |0_L> by a Z on one data qubit before encoding, so the two records agree there.

## The pinned run

The two circuits of a pair must run on the same qubits, or the measured distance mixes the state with the
qubits. `scripts/shor_pinned_ibm_fez.py` builds the Shor round with every circuit on one layout.

```bash
python scripts/shor_pinned_ibm_fez.py dry       # submits nothing
python scripts/shor_pinned_ibm_fez.py submit
python -m syndrome_leakage.analysis_shor_pinned
```

`dry` searches 24 transpiler seeds for the layout with the fewest two-qubit gates and the least summed
two-qubit and readout error, transpiles all eight circuits onto it, asserts that they share their qubits,
two-qubit gates and measurement map, estimates the execution time from the scheduled durations and the
repetition delay, and rehearses the round on the device model. The output goes to `build/pinned_dry_run/`.

`submit` loads the circuits the dry run saved and, before the job and again immediately before it, checks
that the account has exactly one instance, that it is on the open plan with free pricing, and that the
allowance covers twice the estimate; any mismatch stops the script. It caps the execution time, turns
dynamical decoupling and twirling off, since pulses during the idle would move the populations the run
measures, and saves the calibration read immediately before submission and the exact circuits beside the
counts. The token is never read: `QiskitRuntimeService()` loads the saved account.

The run of 2026-09-25 is in [results.md](../results.md#hardware_shor_pinned).
