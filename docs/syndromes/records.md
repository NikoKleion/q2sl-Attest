# Measured records

The leak from two measured records: the plug-in distance and its floor, a permutation test, and the
distance an observer achieves on held-out shots.

```python
from syndrome_leakage import estimate as es

r = es.leak_from_counts(counts0, counts1, n_outcomes=64)   # or leak_from_dists for distributions
es.null_floor(p, shots)                                    # the reading at zero leak
es.shots_for_leak(p0, p1, factor=2.0)                      # shots that put the leak above that floor
print(es.report([("20 us", r)]))
```

Counts may be keyed by bit strings, with or without the spaces Qiskit puts between registers, by hex
strings, or by integers, and a per-shot list or array goes through `dist_from_memory`.

`r["tvd"]` is the plug-in distance, which sits above `r["null"]` even when the two states are identical, so
it is not an estimate of the leak on its own. `r["p_value"]` is a permutation test for any difference
between the records. `r["distance"]` fits the likelihood ratio on half of each record and scores it on the
other half, and `r["distance_corrected"]` takes off the reading the same procedure gives on records drawn
from the pooled counts. Passing `transform=` applies a statistic to each distribution first, such as
`hardware.syndrome_weight_dist`.

The plug-in distance between two sampled distributions is biased upward, and when the true distance is zero
the bias is the whole reading, so it is reported beside the floor of the same shot count. The achieved
distance splits each record, fits the likelihood ratio on one half and scores it on the other, which
removes that bias; it is a lower bound, since the rule is fitted on finite data.
