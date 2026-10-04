# Installation

Python 3.11 or newer and numpy. Every other dependency is optional.

```bash
pip install -e .           # numpy only
pip install -e ".[full]"   # every optional dependency
```

| extra | packages | used by |
|---|---|---|
| `hardware` | qiskit 2.3 or newer, qiskit-aer, qiskit-ibm-runtime | circuits, backends, the device model, the saved QPY circuits |
| `stim` | stim, pymatching | Stim export, `sampled`, matching decoders |
| `tensor` | cotengra 0.7 or newer | `tensor.l2_leak` |
| `qldpc` | qldpc 0.3 or newer | `load.from_qldpc` |
| `wasserstein` | scipy | the Wasserstein distance, the integer programs of `zchecks` |
| `reconstruction` | torch, ldpc | `reconstruction`, BP+OSD decoding |
| `zx` | pyzx | `suite.zx_fingerprint` |
| `plots` | matplotlib | the figures |
| `dev` | pytest | the test suite |

The package installs the command `q2sl` and the packages `syndrome_leakage`, `entropy_fusion`,
`reconstruction` and `suite`.

Inside a clone, `python q2sl.py` works in place of the installed `q2sl`.

## Tests

The suite runs under pytest or with numpy alone.

```bash
python -m pytest tests                     # syndromes and the supporting packages
python tests/run_tests.py
```

A test that needs an optional dependency skips without it. One test reaches qecdb.org and runs only with
the environment variable `Q2SL_NETWORK=1`.
