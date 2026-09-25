# suite.report: Findings to SARIF 2.1.0 with a severity score
import json

from . import core

_RUBRIC = {
    "entropy-trojan":          (9.0, "C:H/I:H", "RNG output with bias or correlation"),
    "secret-reconstruction":   (8.0, "C:H",     "few guesses remain for the structured value"),
    "device-reconstruction":   (8.2, "C:H",     "the structured value is recovered from a device-trained decoder and patterns"),
    "logical-data-leak":       (7.0, "C:H",     "the syndrome distribution depends on the logical state"),
    "gate-identity":           (6.0, "C:L",     "the operation is identifiable from the syndrome"),
    "leakage-events":          (5.0, "A:M/I:M", "leakage out of the computational subspace"),
}


def severity(finding):
    # score 0 to 10, SARIF level, severity vector
    base, vec, _ = _RUBRIC.get(finding.target, (5.0, "C:L", ""))
    if not finding.vulnerable:
        return 0.0, "note", "informational (protection confirmed)"
    score = round(min(10.0, base * finding.confidence), 1)
    level = "error" if score >= 7 else "warning" if score >= 4 else "note"
    return score, level, f"CVSS-like {vec} score {score}/10 (bespoke rubric)"


def to_sarif(findings, tool_version="0.1"):
    rules, seen = [], set()
    for f in findings:
        if f.target not in seen:
            seen.add(f.target)
            base, vec, why = _RUBRIC.get(f.target, (5.0, "C:L", ""))
            rules.append({"id": f.target, "name": f.target,
                          "shortDescription": {"text": f.target},
                          "fullDescription": {"text": why},
                          "properties": {"cvss_like": vec, "impact_base": base}})
    results = []
    for f in findings:
        score, level, vec = severity(f)
        results.append({
            "ruleId": f.target,
            "level": level,
            "message": {"text": f.headline},
            "locations": [{"physicalLocation": {
                "artifactLocation": {"uri": f"simulated://{f.dimension}/{f.module}"}}}],
            "properties": {"dimension": f.dimension, "module": f.module, "metric": f.metric,
                           "confidence": round(f.confidence, 3), "vulnerable": f.vulnerable,
                           "severity_score": score, "severity_vector": vec, "detail": f.detail},
        })
    return {
        "$schema": "https://json.schemastore.org/sarif-2.1.0.json",
        "version": "2.1.0",
        "runs": [{
            "tool": {"driver": {"name": "q2sl suite", "version": tool_version,
                                "informationUri": "https://quantumvillage.org", "rules": rules}},
            "results": results,
        }],
    }


def run_report(target=None, dimension=None, indent=2, **kw):
    findings = core.find(target=target, dimension=dimension, **kw)
    return json.dumps(to_sarif(findings), indent=indent)
