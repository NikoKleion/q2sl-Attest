# Sources

The entropy attestation of quantum random number generation is `qrng_attest`, in `attest/`. It has its own
documentation, tests and results.

| document | contents |
|---|---|
| [attest/README.md](../attest/README.md) | overview, the readout and correlation results, installation |
| [attest/USAGE.md](../attest/USAGE.md) | assessment of a file, certification of a backend, extraction, the SP 800-22 battery |
| [attest/RESULTS.md](../attest/RESULTS.md) | the NIST test vectors, the over-credit sweep, the `ibm_marrakesh` run |

The command line reaches it through either command:

```bash
q2sl attest demo                         # a built-in synthetic source
qrng-attest capture capture.bin --bytes  # a captured file
qrng-attest backend fake_manila --simulate
```

## Relation to the syndromes

Both parts bound what an observer can learn from a quantum device's output. For a syndrome record the
bound is on the information the record carries about the encoded logical state; for a random source it is
on the min-entropy the output holds. Both are reported in the direction that does not over-state them: the
syndrome quantities past the exact engines are lower bounds on the leak, and the attested min-entropy is a
lower bound on the source's entropy.

`entropy_fusion` connects the two: it computes the residual min-entropy of a structured value given
several likelihood sources, one of which can be a device read.
