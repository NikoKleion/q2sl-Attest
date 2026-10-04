# entropy_fusion

Residual min-entropy, in bits, of a synthetic structured value given the likelihoods of several experts.
The score is the sum over positions of -log2 of the largest posterior probability.

- `core.py` Payload and Expert base classes, fusion, scoring, Monte Carlo evaluation, and Shapley values
  (Shapley 1953) over expert subsets.
- `experts.py` field format, learned pattern, simulated hardware read, observed read, revealed positions.
- `payloads.py` a 16 digit card number with a Luhn check digit (Luhn 1960), a Florida phone number, a device
  error string.
- `fast.py` closed-form evaluate and shapley for payloads without a checksum.
- `adaptive.py` fusion with the prior weighted by an estimated trust in [0, 1].
- `fisher.py` the correlation under p of the log-likelihood vectors of two experts, and a redundancy weight
  of 1 minus the largest such correlation with an included expert.
- `factorgraph.py` group factors for correlated positions, such as an area code.
- `registry.py` payload and expert registries.

```python
import entropy_fusion as ef

pay = ef.make_payload("phone_number")
experts = [ef.make_expert("field_format"), ef.make_expert("pattern", pay), ef.make_expert("revealed")]
print(ef.evaluate(pay, experts, reveals=2))
```

Synthetic payloads only.

Checked against independent answers. The fused posterior of a five-digit Luhn number under two experts
equals the enumeration of all 10^5 instances to 1e-12, and the Shapley values add up to the total drop in
residual bits, their defining property.
