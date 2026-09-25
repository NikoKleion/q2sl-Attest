# suite

A registry of modules over the packages. Each module has a dimension, the input it reads, and a target,
what it reports. Output is SARIF 2.1.0 (OASIS 2020).

- modules: `logical_data_leak`, `syndrome_eavesdropper`, `secret_reconstruction`, `gate_fingerprint`,
  `physical_leakage`, `device_reconstruction`, `entropy_trojan`
- `report.py` a severity score per target from a fixed table
- `zx_fingerprint.py` circuit identity by ZX-calculus (Coecke and Duncan 2011) through PyZX (Kissinger and
  van de Wetering 2020)
- `covert_channel.py` one bit carried by data-dependent noise, read from syndrome shots

```python
import suite

for f in suite.find():
    print(f)
print(suite.report.run_report())
```

```bash
python -m suite
python -m suite find gate-identity
```
