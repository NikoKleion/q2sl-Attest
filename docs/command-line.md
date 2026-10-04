# Command line

`q2sl` covers the syndromes and the supporting packages. Inside a clone, `python q2sl.py` works in
place of the installed command.

| command | does |
|---|---|
| `q2sl assess <code>` | the analyses for a code name, `file:<path>`, `qecdb:<id>` or `backend:<name>` |
| `q2sl qecdb n=9 k=1 d=3` | search qecdb.org for CSS codes by `n`, `k`, `d` and `family`, each a value or a range `lo-hi` |
| `q2sl leak` | the syndrome self-test |
| `q2sl eavesdrop` | the likelihood-ratio test on syndrome records |
| `q2sl fusion` | residual min-entropy on synthetic payloads |
| `q2sl reconstruct` | the full reconstruction, with torch and ldpc |
| `q2sl example` | the end to end run in `scenario/one_posterior.py` |
| `q2sl suite` | every suite module, with a SARIF summary |
| `q2sl all` | the numpy demonstrations and the suite |

Code names: `repetition`, `five_qubit`, `steane`, `code_4_1_2`, `hamming_7`, `shor`.

```bash
q2sl assess shor
q2sl assess file:mycode.npz
q2sl assess qecdb:674f2504f9caaa7ce7667423
q2sl assess backend:manila steane
q2sl qecdb n=25 k=1 d=5 family=surface
```

`assess` reports the size, the distance and the leak order, then the leak and the logical error for codes
small enough to simulate. For a CSS code that leaks it names the fewest Z generators whose joint record leaks
and the distance of that record, where scipy is installed. Past about 14 qubits it gives the leak order by
integer programming, where scipy is installed, and the contracted L2 distance with the lower bound it gives,
where cotengra is installed.