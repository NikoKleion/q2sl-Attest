# the Stim text: the grammar, the twirl, and every detector deterministic on an independent simulator
from _harness import needs
from syndrome_leakage import expectations as ex
from syndrome_leakage import hardware as hw
from syndrome_leakage import stim_export as se
from syndrome_leakage.codes import STANDARD

KNOWN = {"R", "H", "X", "CX", "CY", "CZ", "M", "TICK", "DETECTOR", "OBSERVABLE_INCLUDE", "PAULI_CHANNEL_1"}


def _annotations(text):
    n_meas, detectors, observable = 0, [], None
    for name, args, targets in se.parse(text):
        assert name in KNOWN, name
        if name in ("CX", "CY", "CZ"):
            assert len(targets) % 2 == 0, (name, targets)
        if name == "M":
            n_meas += len(targets)
        if name in ("DETECTOR", "OBSERVABLE_INCLUDE"):
            recs = [int(t[4:-1]) for t in targets]
            assert all(-n_meas <= r < 0 for r in recs), (recs, n_meas)
            if name == "DETECTOR":
                detectors.append([n_meas + r for r in recs])
            else:
                observable = [n_meas + r for r in recs]
    return n_meas, detectors, observable


def test_twirled_damping_is_a_pauli_channel():
    for g in (0.0, 0.01, 0.3, 1.0):
        px, py, pz = se.twirled_damping(g)
        assert px == py == g / 4
        assert 0 <= pz and px + py + pz <= 1 + 1e-12
    px, _py, pz = se.twirled_damping(0.1)
    assert abs(pz - 0.0006583509747) < 1e-12


def test_text_follows_the_grammar():
    for code in (STANDARD["repetition"](), STANDARD["steane"](), hw.shor_code(), ex.surface_code_3()):
        for flags in (False, True):
            text = se.stim_text(code, prep=1, flags=flags, gamma=0.05, readout=0.01)
            n_meas, detectors, observable = _annotations(text)
            z_checks = len(hw.select_checks(code, "z"))
            flagged = sum(1 for g in code.stab_strings if sum(c != "I" for c in g) >= 3) if flags else 0
            assert n_meas == len(code.stab_strings) + flagged + code.n
            assert len(detectors) == len(code.stab_strings) + flagged + z_checks
            assert observable and len(observable) == sum(ch == "Z" for ch in code.zl_str)


def test_every_detector_is_deterministic_without_noise():
    needs("qiskit")
    needs("qiskit_aer")
    from qiskit import transpile
    from qiskit_aer import AerSimulator
    sim = AerSimulator(method="statevector")
    for code in (STANDARD["steane"](), ex.surface_code_3()):
        m, n = len(code.stab_strings), code.n
        flagged = [j for j, g in enumerate(code.stab_strings) if sum(c != "I" for c in g) >= 3]
        for prep in (0, 1):
            n_meas, detectors, observable = _annotations(se.stim_text(code, prep=prep, flags=True))
            qc = hw.build_extraction_circuits(code, [0.0], flags=True, final_data=True)[0][prep]
            counts = sim.run(transpile(qc, sim), shots=64, seed_simulator=31_000_001 + prep).result().get_counts()
            seen, obs = [set() for _ in detectors], set()
            for (syn, flag, data), _n in hw.split_registers(counts, (m, m, n)).items():
                record = [(syn >> j) & 1 for j in range(m)] + [(flag >> j) & 1 for j in flagged]
                record += [(data >> q) & 1 for q in range(n)]
                assert len(record) == n_meas
                for i, det in enumerate(detectors):
                    seen[i].add(sum(record[r] for r in det) % 2)
                obs.add(sum(record[r] for r in observable) % 2)
            assert all(len(s) == 1 for s in seen), (code.name, prep)
            assert obs == {prep}, (code.name, prep, obs)


def test_stim_reads_the_text_and_pymatching_decodes_it():
    needs("stim")
    needs("pymatching")
    import pymatching
    import stim
    for code in (STANDARD["repetition"](), STANDARD["steane"](), hw.shor_code(), ex.surface_code_3()):
        for flags in (False, True):
            for prep in (0, 1):
                c = stim.Circuit(se.stim_text(code, prep=prep, flags=flags))
                det, obs = c.compile_detector_sampler(seed=5).sample(64, separate_observables=True)
                assert det.sum() == 0 and obs.sum() == 0, (code.name, flags, prep)
        noisy = stim.Circuit(se.stim_text(code, prep=0, gamma=0.05, readout=0.01))
        dem = noisy.detector_error_model(decompose_errors=True)
        det, obs = noisy.compile_detector_sampler(seed=11).sample(20000, separate_observables=True)
        pred = pymatching.Matching.from_detector_error_model(dem).decode_batch(det)
        assert (pred != obs).any(axis=1).mean() < obs.mean(), code.name


def test_a_pauli_twirl_shows_no_leak_in_stim():
    # the Pauli boundary on an independent simulator: the two logical states give one distribution
    needs("stim")
    import numpy as np
    import stim
    code = hw.shor_code()

    def dist(prep, seed):
        c = stim.Circuit(se.stim_text(code, prep=prep, gamma=0.2))
        det = c.compile_detector_sampler(seed=seed).sample(200000)
        keys = (det @ (1 << np.arange(det.shape[1], dtype=np.int64))).astype(np.int64)
        vals, counts = np.unique(keys, return_counts=True)
        return dict(zip(vals.tolist(), (counts / len(keys)).tolist()))

    def tvd(a, b):
        return 0.5 * sum(abs(a.get(k, 0) - b.get(k, 0)) for k in set(a) | set(b))

    across = tvd(dist(0, 1_000_003), dist(1, 2_000_006))
    floor = tvd(dist(0, 21_000_009), dist(0, 28_000_012))
    assert across < 2.5 * floor, (across, floor)
    assert across < 0.03
