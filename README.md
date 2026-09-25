<div align="center">

<img src="assets/logo_readme.png" width="360" alt="q2sl attest">

**Quantum Syndrome and Source Limits**

[![python](https://img.shields.io/badge/python-3.11%2B-30363D?style=flat-square)](https://www.python.org/)
[![tests](https://img.shields.io/badge/tests-320%20passing-30363D?style=flat-square)](docs/testing.md)
[![version](https://img.shields.io/badge/version-1.0-30363D?style=flat-square)](pyproject.toml)
[![license](https://img.shields.io/badge/license-Apache%202.0-30363D?style=flat-square)](LICENSE)

</div>

q2sl attest bounds what an observer can learn from the output of a quantum device, in two settings: the
syndrome record of a stabilizer code, and the output of a quantum random number generator. It is written
in Python on numpy, with Qiskit for circuits, backends and hardware runs.

## Syndromes

Under a Pauli channel the syndrome distribution of a stabilizer code is the same for every encoded logical
state. Under a channel that is not Pauli, such as amplitude damping, it is not, and the syndrome record
carries information about the logical state. `syndrome_leakage` computes that dependence. In this project
"leakage" refers to this information, not to population leaving the computational states.

- The distance between the syndrome distributions of two logical states, exactly up to about 14 qubits.
- The order of that distance in the noise strength, by enumeration and by integer programming, for codes
  up to the [[144,12,12]] bivariate bicycle code.
- The L2 distance by tensor network contraction, a lower bound on the leak, for planar codes up to a
  distance 11 surface code.
- Syndrome records sampled under T1 and T2 relaxation at any code size, and lower bounds read from them.
- The likelihood-ratio test on syndrome records, the number of rounds an observer needs, and the logical
  error after recovery.
- Estimates from measured records, extraction circuits for any CSS code, a gate-level device model, Stim
  export, and codes loaded from files, qecdb.org, qLDPC and ldpc.

## Sources

`qrng_attest`, in [attest/](attest/), assesses the min-entropy of a random source under NIST SP 800-90B,
attests a quantum random number generator from its device physics and readout calibration, and extracts
uniform bits. It is also maintained as the standalone repository `NikoKleion/qrng-attest`.

## Hardware measurement

![the pinned Shor run on ibm_fez](results/figures/hardware_pinned.png)

The Shor [[9,1,3]] code with one round of its six Z generators on IBM `ibm_fez`, job
`daquif6ekp0c73arbd70`, 16000 shots per circuit, every circuit on one layout. The two logical states'
records agree at zero delay and differ at 20, 50 and 100 us, where an observer achieves 0.025, 0.163 and
0.230 on held-out shots, each with p below 0.0005. The dashed line is the gate-level model built from the
calibration read before submission. Full numbers in [docs/results.md](docs/results.md#hardware_shor_pinned).

## Installation

```bash
pip install -e .           # numpy only
pip install -e ".[full]"   # every optional dependency
```

Python 3.11 or newer. The optional dependencies and what uses them are listed in
[docs/installation.md](docs/installation.md).

## Usage

```bash
q2sl assess shor                           # the analyses for one code
q2sl assess file:mycode.npz                # a code from a file
q2sl qecdb n=9 k=1 d=3                     # search qecdb.org for CSS codes
q2sl assess qecdb:674f2504f9caaa7ce7667423
q2sl attest demo                           # entropy attestation of a synthetic source
python -m syndrome_leakage.experiments --save
python -m syndrome_leakage.figures
```

```python
from syndrome_leakage import expectations as ex, zchecks
from syndrome_leakage.channels import amplitude_damping

code = ex.surface_code_3()
leak, d0, d1 = ex.population_leak(code, amplitude_damping(0.1))
order = zchecks.leak_order(code)["order"]
```

## Documentation

- [docs/README.md](docs/README.md) the documentation, by topic
- [docs/results.md](docs/results.md) every validation run, with its control
- [CONTRIBUTING.md](CONTRIBUTING.md) the tests and conventions a change is held to
- [SECURITY.md](SECURITY.md) reporting a vulnerability, or a code reported as protected that leaks

## Support

Maintained by Nikolas Klein. Questions and bug reports go to the issue tracker and are answered on a
best-effort basis; security reports follow [SECURITY.md](SECURITY.md).

## Citation

[CITATION.cff](CITATION.cff) gives the citation.

## License

Apache License 2.0, see [LICENSE](LICENSE). Copyright 2026 Nikolas Klein. Third party material and its
licences are listed in [NOTICE](NOTICE). Qiskit is a trademark of International Business Machines
Corporation; this project is not affiliated with or endorsed by IBM.
