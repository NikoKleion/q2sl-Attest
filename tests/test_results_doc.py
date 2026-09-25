# numbers in the docs/results.md tables against the saved files in results/
import os
import re
from decimal import Decimal

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESDIR = os.path.join(ROOT, "results")
NUMBER = re.compile(r"(?<![\w.\[])[-+]?\d+\.\d+(?:e[-+]?\d+)?|(?<![\w.\[])\d+e[-+]?\d+|(?<![\w.\[,])\d+(?![\w.\]])")
FILES = {"hardware": ["hardware_ibm_fez.txt", "hardware_analysis.txt"],
         "hardware_shor": ["hardware_shor_analysis.txt"],
         "hardware_shor_pinned": ["hardware_shor_pinned_analysis.txt"]}


def sections():
    doc = open(os.path.join(ROOT, "docs", "results.md"), encoding="utf-8").read()
    out = {}
    for part in doc.split("\n## ")[1:]:
        name, body = part.split("\n", 1)
        out[name.strip()] = body
    return out


def table_numbers(body):
    for line in body.split("\n"):
        if not line.startswith("|") or set(line) <= set("|-: "):
            continue
        cells = line.strip("|").split("|")
        if all(not NUMBER.search(c) for c in cells[1:]):
            continue
        for cell in cells[1:]:
            yield from NUMBER.findall(cell.split("+/-")[0])


def test_table_numbers_are_in_the_saved_output():
    checked = 0
    for name, body in sections().items():
        files = FILES.get(name, [f"{name}.txt"])
        if not all(os.path.exists(os.path.join(RESDIR, f)) for f in files):
            continue
        text = "".join(open(os.path.join(RESDIR, f), encoding="utf-8").read() for f in files)
        saved = [float(v) for v in NUMBER.findall(text)]
        for token in table_numbers(body):
            d = Decimal(token)
            tol = 0.5 * 10.0 ** d.as_tuple().exponent
            assert any(abs(v - float(d)) <= tol * (1 + 1e-9) for v in saved), \
                f"docs/results.md section {name} shows {token}, not found in {files}"
            checked += 1
    assert checked > 50


def test_every_run_has_a_section():
    names = {f[:-4] for f in os.listdir(RESDIR) if f.endswith(".txt")}
    names -= {"hardware_ibm_fez", "hardware_analysis", "hardware_shor_analysis", "hardware_shor_pinned_analysis"}
    assert names <= set(sections())
