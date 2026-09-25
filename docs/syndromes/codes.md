# Codes

## Load a code from a file

```python
from syndrome_leakage import load

load.load_code("code.npz")                 # keys hx and hz
load.load_code("hx.txt", "hz.alist")       # one file per matrix
load.load_code("generators.txt")           # Pauli strings, one per line
load.from_css_code(obj)                    # anything carrying hx and hz, as bposd's css_code does
```

When the object also carries `lz` and `lx`, its Z logicals are kept exactly, so logical qubit i here is
logical qubit i there. bposd builds the two sets separately and does not pair them, its `lx lz^T` being
full rank and not the identity, so the X logicals are re-paired inside their own span to give
<Z_i, X_j> = delta_ij. `css_strings(Hx, Hz, lx=..., lz=...)` takes the same two matrices directly, and
`use_logicals=False` ignores them.

`read_matrix` takes 0/1 text, `.npy`, alist and MatrixMarket. A text file holding Hx and Hz separated by a
blank line loads on its own. The command line reads the same files:

```bash
q2sl assess file:code.npz
q2sl assess file:hx.txt,hz.alist
```

It reports the size, distance and leak order from the strings, then the leak and the logical error for
codes small enough to simulate, and says which of those it could reach. Past about 14 qubits it prints the
contracted L2 distance and the bound it gives, when cotengra is installed. Where the order search is out of
reach it fits the exponent of the L2 distance over damping 1e-4 to 1e-6 and reports it as the leak order
when the fit lies within 0.05 of an integer; on the d=7 surface code it reads 7.000 in about five seconds.

## Load a code from qecdb.org

```python
from syndrome_leakage import load
rows = load.qecdb_search(n=9, k=1, d=3)        # CSS codes; n, k and d take a value or a (low, high) pair
code = load.from_qecdb(rows[0]["_id"])
```

[qecdb.org](https://qecdb.org) catalogues about 37,000 published stabilizer codes, 15,000 of them CSS, as
Pauli strings. `qecdb_search` lists CSS codes by `n`, `k`, `d` and `family`, and `from_qecdb` fetches one
and builds it. A record's logicals come as symplectic pairs whose X and Z order varies from pair to pair,
so they are sorted by type rather than read by position. A code that is not CSS is refused. Nothing is
bundled: each call reaches qecdb.org, whose JSON interface is undocumented and carries no stated stability
guarantee. From the command line:

```bash
q2sl qecdb n=9 k=1 d=3                   # search; family=surface narrows it
q2sl assess qecdb:674f2504f9caaa7ce7667423
```

## Load a code from qLDPC

```python
from qldpc import codes
from syndrome_leakage import load
code = load.from_qldpc(codes.SurfaceCode(5))     # load.from_css_code(obj) recognises it too
```

qLDPC already pairs its logical operators, so its basis is kept as it is and logical qubit i here is
logical qubit i there. Only qubit CSS codes are read. It needs qLDPC: `pip install -e ".[qldpc]"`.
