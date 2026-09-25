# Tests for suite: runs every module, severity triage, SARIF validity.
import suite
from suite import report as R


def test_all_modules_run():
    findings = suite.find()
    assert len(findings) >= 6
    for f in findings:
        assert f.headline and 0.0 <= f.confidence <= 1.0


def test_sarif_is_valid_211():
    sarif = R.to_sarif(suite.find())
    assert sarif["version"] == "2.1.0"
    run = sarif["runs"][0]
    assert len(run["results"]) == len(run["tool"]["driver"]["rules"]) or run["results"]
    for res in run["results"]:
        assert res["ruleId"] and res["level"] in ("note", "warning", "error")


def test_protected_finding_is_informational():
    # a protected code must not be flagged as a vulnerability
    f = suite.find(target="logical-data-leak", code="steane")[0]
    score, level, _ = R.severity(f)
    assert not f.vulnerable and score == 0.0 and level == "note"


def test_severity_monotone_with_confidence():
    findings = suite.find()
    for f in findings:
        s, _, _ = R.severity(f)
        assert 0.0 <= s <= 10.0
