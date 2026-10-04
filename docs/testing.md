# Testing

The suite runs under pytest or, with numpy alone, under its own runner; a test that needs an optional
dependency skips without it.

```bash
python -m pytest tests                  # syndromes and the supporting packages
python tests/run_tests.py
```

Continuous integration runs it with numpy alone on Python 3.11, 3.12 and 3.13, and with every
optional dependency on 3.11 and 3.13, where no test may skip for a missing dependency. One test reaches
qecdb.org and runs only with `Q2SL_NETWORK=1`.

## Files in `tests/`

- `test_syndrome_leakage.py` analytic against exact, the CSS path, Wasserstein, the T1 and T2 channel.
- `test_pauli_boundary.py`, `test_kl_check.py`, `test_eavesdrop.py`, `test_worst_pair.py`,
  `test_coherent_diagonal.py`, `test_scaling.py`, `test_device_grounded.py`, `test_held_memory.py` one file
  per run in [results.md](results.md).
- `test_closed_forms.py` the exact engine against closed forms, including equation 93 of Hu, Liang and
  Calderbank.
- `test_expectations.py` the stabilizer-expectation engine against the projector engine.
- `test_approx_qec.py` branch matrices against equation 39 of Leung et al.
- `test_hlc.py` equation 91 of Hu, Liang and Calderbank against the exact simulation.
- `test_tensor.py` the contraction against the exact engine, the precision at small damping, the width
  limit, and a reproducible plan.
- `test_sampled.py` the reset and Z mixture against the relaxation channel, the records against the exact
  distribution, the Pauli boundary, and the seed.
- `test_zchecks.py` the classical process against the exact engine and the enumeration, T2 independence,
  the integer-program order against the enumeration, the distance, and the decoders.
- `test_estimate.py` the finite-sample tools: calibration of the permutation test, a known leak recovered,
  zero leak reading zero.
- `test_hardware_bridge.py` the device bridge on a stub backend.
- `test_gate_level.py`, `test_stim_export.py`, `test_extraction.py` the device model and its parts made
  noiseless, the Stim text and the extraction circuits.
- `test_multi_logical.py`, `test_rounds.py`, `test_protection.py` codes with k > 1, many rounds, recovery.
- `test_decoded.py` the least weight decoder against every Pauli, the views against the syndrome leak, the
  frame bit over representatives, the Shor closed form, and the pinned run's parity bit.
- `test_regions.py` the record of a set of generators against the whole record summed over the rest, the
  smallest leaking sets, the rule that a set leaks when a product covers a logical, the region engine against
  the exact engine and the distance 5 Z record, a region of the distance 11 code, sets of X and Z generators
  against the exact engine and under T1 and T2, the drops, their sign, and the cover of least-weight logicals.
- `test_load.py`, `test_qecdb.py`, `test_qldpc.py` codes from files, qecdb.org records and qLDPC objects.
- `test_shor_hardware.py` the Shor circuits, the saved ibm_fez records, and the pinned run's records and
  circuits.
- `test_entropy_fusion.py` closed-form scoring against Monte Carlo, the fused posterior against
  enumeration, the Shapley efficiency property.
- `test_reconstruction.py` the toric device's code parameters, matching against the published threshold,
  and the stabilizer Renyi entropy.
- `test_suite.py`, `test_q2sl.py` the suite modules and the command line.
- `test_figures.py` the figures against the tables they draw from.
- `test_results_doc.py` the numbers in [results.md](results.md) against the files in `results/`.
