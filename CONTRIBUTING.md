# Contributing

Issues and pull requests are welcome.

## Running the tests

```bash
pip install -e ".[dev]"
python -m pytest tests                    # syndromes and the supporting packages
python tests/run_tests.py                 # the same suite without pytest
```

With every optional dependency, which is what continuous integration runs:

```bash
pip install -e ".[full]"
python -m pytest tests -rs                # no test should skip for a missing dependency
```

## Conventions

- Under any Pauli channel the syndrome distribution does not depend on the logical state, so the leak is
  exactly zero. `tests/test_pauli_boundary.py` pins this. A change that makes a Pauli channel leak is a
  bug.
- Every number in a `docs/results.md` table is checked against the saved output in `results/` by
  `tests/test_results_doc.py`. A change that moves a number regenerates the run with
  `python -m syndrome_leakage.experiments <run> --save` and updates the table in the same commit.
- A quantity computed approximately is reported as a bound or an estimate with its interval, and named as
  one. The exact engine refuses rather than approximates past its ceiling.
- Saved output contains no wall clock time or other value that changes between runs of the same code.
- Code or data taken from another project is credited in NOTICE with its licence, and at the point where
  it is used.
- Code carries one-line headers and markers. Methods, derivations and usage belong in `docs/`.
- Hardware runs are never part of the test suite. Their records are saved under `results/` and the tests
  read the saved records.
