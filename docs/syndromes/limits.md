# Limits

Every run except the one in `hardware.py` is a simulation. Device parameters enter as numbers, so a result
at given T1 and T2 is a statement about that noise model.

`Code` builds all 2^(n-k) syndrome projectors, each 2^n by 2^n, so codes past about eleven qubits cannot
be constructed; memory sets the limit. A code with k > 1 is read one logical qubit at a time, and the
reading depends on the logical basis. `expectations.syndrome_dist` needs 2^m Pauli expectations of cost
n 2^n each and reaches about fourteen qubits, for pure logical states and one channel on every qubit.

The analytic leak order is compared with the exact simulation on codes up to nine qubits, and with the
exponent fitted from the contraction on surface d=5 and the hypergraph product of rep4, where both engines
give 5 and 4. The contraction returns the L2 distance and a lower bound on the leak, not the leak, reads
CSS codes only, and refuses past its width limit. Sampling covers relaxation with T2 <= T1 under ideal
stabilizer measurement, and each statistic of its records bounds the leak from below only. `zchecks` reads
the Z checks alone, under relaxation alone, with ideal stabilizer measurement.

A drop from `qubit_drops` ranks qubits at weak damping only: past a damping of about 0.4 on the codes run,
making one qubit noiseless raises the leak for some qubits, and at 0.6 for every qubit of every code run. A
qubit on no Z logical of least weight has a negative drop at every damping.

`regions` returns the leak of a part of the record, which is at most the leak of the whole record and can be
zero where the whole is not. `z_region_leak`, `smallest_leaking_set` and `region_order` read Z generators of
CSS codes under relaxation with ideal measurement. `css_region_leak` reads X and Z generators of a CSS code
under relaxation alone, a decay probability below 1 and a coherence factor per qubit; under any other channel
a set of generators needs the exact engine. The integer program proves its answer where the solver finishes,
108 s for the distance 11 surface code, and otherwise reports that it stopped.

Leak measures compare the two logical basis states. `worst_case_pair` searches a grid of 58 logical
states for a more distinguishable pair.

The test uses the two state-conditioned distributions, so the code, the channel and the noise strength
are taken as known.
`rounds_for_error` is an asymptotic figure from the Chernoff exponent. At small round counts `ml_attack`
gives the error directly.

Summing per-round exponents gives a partial sum at that round count. It is a total only if the per-round
values have gone to zero, which the `rounds` run measures: the order-3 codes stop moving by a hundred
rounds and the repetition code is still adding at three thousand.

`rounds.iterate` needs a decoder that explains every syndrome, so that each branch returns to the code
space and the round becomes a map on the logical qubit. Readout error breaks that, since a correction
chosen from a flipped syndrome leaves the state outside the code space, so this path takes ideal
measurement only.

`protection.py` corrects single-qubit Paulis after one round with ideal syndrome measurement.

Aer draws its shots from `seed_simulator`, and two runs whose seeds are close share random numbers: on the
Shor circuits, seeds 5000 and 5001 give identical counts and 5037 gives a distance of 0.0007, against 0.023
for seeds a thousand apart. Comparing two states means two runs, so give them seeds spaced by more than the
shot count, or the measured distance carries none of the sampling noise a device would show.

The channel library applies one single-qubit channel to each qubit independently, the same one everywhere
or one per qubit. Correlated noise, crosstalk and leakage out of the computational subspace are not
modelled. Measurement error enters as a flip on each syndrome bit through `readout=`, which covers the
readout of the ancillas; error on the data qubits during extraction is not part of it, and the `readout`
run shows what that leaves out.
