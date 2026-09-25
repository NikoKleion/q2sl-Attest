# the figures are drawn from the same files as the results tables, and agree with them
import os

import numpy as np

from _harness import needs

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# docs/results.md, shor_device_model, at 0, 20, 50 and 100 us
MEASURED = [0.0108, 0.0988, 0.2233, 0.2049]
TODAY = [0.0097, 0.0613, 0.2403, 0.3164]
JOB_DAY = [0.0113, 0.0472, 0.1921, 0.3235]
STATE = [0.0666, 0.2315, 0.3163]
LAYOUT = [0.0648, 0.0637, 0.0571]


def _figures():
    needs("matplotlib")
    from syndrome_leakage import figures as fg
    return fg


def test_hardware_distances_match_the_results_table():
    fg = _figures()
    measured, today, job_day = fg._load()
    for store, want in ((measured, MEASURED), (today, TODAY), (job_day, JOB_DAY)):
        got = [fg._distance(store, d) for d in fg.DELAYS]
        for a, b in zip(got, want):
            assert abs(a - b) < 5e-5, (got, want)


def test_state_and_layout_split_matches_the_results_table():
    fg = _figures()
    _measured, today, _job = fg._load()
    from syndrome_leakage.estimate import dist_from_counts
    for i, d in enumerate([x for x in fg.DELAYS if x]):
        twin = dist_from_counts(today[f"{d}_0_on_1"], fg.OUTCOMES)[0]
        state = fg._tvd(twin, fg._dist(today, d, 1))
        layout = fg._tvd(fg._dist(today, d, 0), twin)
        assert abs(state - STATE[i]) < 5e-5, (d, state)
        assert abs(layout - LAYOUT[i]) < 5e-5, (d, layout)


def test_the_pinned_figure_matches_the_results_table():
    # docs/results.md, hardware_shor_pinned: measured, model and achieved at 0, 20, 50 and 100 us
    needs("matplotlib")
    import tempfile
    from syndrome_leakage import figures
    measured, model, achieved = figures.pinned_figure(os.path.join(tempfile.mkdtemp(), "p.png"))
    assert np.allclose(measured, [0.0117, 0.0438, 0.1673, 0.2337], atol=5e-5)
    assert np.allclose(model, [0.0088, 0.0478, 0.2483, 0.3446], atol=5e-5)
    assert np.allclose(achieved, [0.0000, 0.0254, 0.1627, 0.2304], atol=5e-5)


def test_the_precision_figure_matches_the_results_table():
    needs("matplotlib")
    import tempfile
    from syndrome_leakage import figures
    ratio = figures.precision_figure(os.path.join(tempfile.mkdtemp(), "p.png"))
    assert np.allclose(ratio, [1.114906, 1.117721, 1.118003, 1.118031, 1.118034, 1.118034], atol=5e-7)


def test_the_readme_figures_are_present():
    figdir = os.path.join(ROOT, "results", "figures")
    for name in ("fig1_hardware_shor.png", "fig2_leak_order.png", "hardware_pinned.png", "precision.png",
                 "records.png"):
        path = os.path.join(figdir, name)
        assert os.path.exists(path) and os.path.getsize(path) > 10_000, path
