# tests for the top-level CLI and the assess dispatch.
import q2sl
from suite.core import Finding
from suite.modules import DeviceReconstruction


def test_cli_help_runs():
    q2sl.main(["help"])


def test_assess_help_runs():
    q2sl.main(["assess"])


def test_assess_code_leaking_and_protected():
    q2sl._assess_code("repetition")
    q2sl._assess_code("steane")


def test_reconstruction_module_optout_contract():
    # returns a Finding when torch and ldpc are present, else None; never crashes
    r = DeviceReconstruction().run(samples=4)
    assert r is None or isinstance(r, Finding)


def test_find_never_returns_none():
    import suite
    assert all(f is not None for f in suite.find())
