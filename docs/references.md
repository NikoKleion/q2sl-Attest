# References

## external tools

Each is optional and used when installed. `syndrome_leakage.hardware.available()` reports whether qiskit is
present.

- Qiskit and qiskit-ibm-runtime: `syndrome_leakage/hardware.py` reads T1, T2 and gate durations and runs
  the circuits that measure the leak; `reconstruction.device_from_backend` reads error rates.
- ldpc: BP+OSD in `reconstruction/`.
- torch: the pattern net and the neural decoder in `reconstruction/`.
- PyZX: `suite/zx_fingerprint.py`.
- scipy: the exact Wasserstein distance in `syndrome_leakage/wasserstein.py`.

## references

- Gottesman, Stabilizer codes and quantum error correction, arXiv:quant-ph/9705052.
- Knill and Laflamme, Phys. Rev. A 55, 900 (1997): the error-correction conditions.
- Leung, Nielsen, Chuang and Yamamoto, Phys. Rev. A 56, 2567 (1997), arXiv:quant-ph/9704002: approximate
  error-correction conditions and the [[4,1,2]] code under amplitude damping.
- Nielsen and Chuang, Quantum Computation and Quantum Information, section 8.3: amplitude damping and
  relaxation channels.
- Kretschmann, Kribs and Spekkens, arXiv:0711.3438: private and correctable subsystems.
- Cover and Thomas, Elements of Information Theory, section 11.9: the Chernoff information.
- Neyman and Pearson, Phil. Trans. R. Soc. A 231, 289 (1933): the likelihood-ratio test.
- Wagner, Kampermann, Bruss and Kliesch, Quantum 6, 809 (2022), arXiv:2107.14252: stabilizer expectations
  under a Pauli channel in terms of the Pauli error distribution (eq. 40).
- Kobori and Todo, Phys. Rev. A (2025), arXiv:2406.08981: syndrome probabilities from stabilizer
  expectations (eq. 15), and estimating non-Pauli noise from syndrome statistics.
- Blume-Kohout and Young, arXiv:2504.14643: detector error models from syndrome data; polarizations and
  their Walsh-Hadamard inversion (eqs. 7 and 34).
- Hu, Liang and Calderbank, arXiv:2109.13481: syndrome probabilities under a diagonal unitary, and the
  condition for them to be independent of the encoded state.
- Shukla, Browne and Nishio, arXiv:2607.12174: syndrome data as a decoder side channel.
- Kitaev, Ann. Phys. 303, 2 (2003): the toric code.
- Shor, Phys. Rev. A 52, R2493 (1995): the [[9,1,3]] code.
- Gidney, Quantum 5, 497 (2021), arXiv:2103.02202: Stim; the definition of a detector (section 5.6) and
  Pauli-frame noise (section 2.4.1).
- Bravyi, Cross, Gambetta, Maslov, Rall and Yoder, Nature 627, 778 (2024), arXiv:2308.07915: the bivariate
  bicycle codes built in the `tensor` and `sampled` runs.
- Huangfu and Hall, Math. Prog. Comp. 10, 119 (2018): HiGHS, the solver behind `scipy.optimize.milp`.
- Higgott and Gidney, arXiv:2303.15933: sparse blossom, the matching algorithm of PyMatching 2.
- Perlin, qLDPC, <https://github.com/qLDPCOrg/qLDPC>: its `BBCode` builds the same check matrices as the
  runs here, and its randomized distance bounds give 6 and 12.
- Qiskit Aer, `noise.thermal_relaxation_error`: relaxation with T2 <= T1 expressed as a mixture of reset
  and unitary errors, the decomposition `sampled` uses.
- Chao and Reichardt, Phys. Rev. Lett. 121, 050502 (2018), arXiv:1705.02329: the flag circuit, Fig. 2(c).
- Chao and Reichardt, PRX Quantum 1, 010302 (2020), arXiv:1912.09549: flags at arbitrary distance, d of
  them for a distance-d code.
- Roffe, White, Burton and Campbell, Phys. Rev. Research 2, 043423 (2020), arXiv:2005.07016: BP+OSD and the
  bposd package, whose `css_code.compute_logicals` sets the logical-operator convention read here.
- Mithra, QECops v2.5 (2026), doi:10.5281/zenodo.19410365, MIT licensed,
  <https://github.com/JitheshMithra/QECops>: the Shor-code syndrome-extraction circuit layout ported
  in `hardware.build_shor_circuits`, and the quenched per-qubit noise model that the `disorder` run
  uses. See [NOTICE](../NOTICE).
- Panteleev and Kalachev, Quantum 5, 585 (2021): BP+OSD decoding.
- Roffe, White, Burton and Campbell, Phys. Rev. Research 2, 043423 (2020): BP+OSD and the ldpc package.
- Leone, Oliviero and Hamma, Phys. Rev. Lett. 128, 050402 (2022): stabilizer Renyi entropy.
- Wang, Harrington and Preskill, Ann. Phys. 303, 31 (2003): the toric code threshold under bit flips.
- Coecke and Duncan, New J. Phys. 13, 043016 (2011): ZX-calculus.
- Kissinger and van de Wetering, EPTCS 318, 229 (2020): PyZX.
- Shapley, A value for n-person games, Contributions to the Theory of Games II (1953).
- Luhn, US Patent 2,950,048 (1960): the Luhn check digit.
- OASIS, Static Analysis Results Interchange Format (SARIF) Version 2.1.0 (2020).
