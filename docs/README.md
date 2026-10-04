# Documentation

q2sl attest computes, in `syndrome_leakage`, the information a stabilizer code's syndrome carries about the
encoded logical state under noise that is not Pauli, its order in the noise strength, the attack an
observer of the syndrome records achieves, and estimates from measured records.

## Pages

| page | contents |
|---|---|
| [installation.md](installation.md) | requirements, optional dependencies, the two test suites |
| [command-line.md](command-line.md) | the `q2sl` command |
| [syndromes/concepts.md](syndromes/concepts.md) | the quantities measured and the conventions |
| [syndromes/analysis.md](syndromes/analysis.md) | a code under a channel: leak, attack, recovery, rounds |
| [syndromes/engines.md](syndromes/engines.md) | the six engines, their methods and their reach |
| [syndromes/codes.md](syndromes/codes.md) | codes from files, qecdb.org and qLDPC |
| [syndromes/records.md](syndromes/records.md) | the leak from measured records |
| [syndromes/hardware.md](syndromes/hardware.md) | extraction circuits, the device model, Stim export, hardware runs |
| [syndromes/limits.md](syndromes/limits.md) | what the models assume |
| [entropy_fusion.md](entropy_fusion.md) | residual min-entropy of a structured value |
| [reconstruction.md](reconstruction.md) | one posterior over a synthetic value from several sources |
| [suite.md](suite.md) | the module registry and SARIF output |
| [results.md](results.md) | every validation run, with its control |
| [testing.md](testing.md) | the test suites and what they check |
| [references.md](references.md) | papers and software |

## Figures

The figures are drawn from the saved output in `results/` by `python -m syndrome_leakage.figures` and are
checked against the tables in [results.md](results.md) by `tests/test_figures.py`.
