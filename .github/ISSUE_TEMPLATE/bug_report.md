---
name: Bug report
about: A command or function that fails, or returns a result you believe is wrong
labels: bug
---

**Version**
Output of `python -c "import syndrome_leakage; print(syndrome_leakage.__version__)"`:

**Python version, platform, and installed extras**

**What you ran**
The command, or the smallest piece of code that reproduces it.

**What it returned**

**What you expected, and why**
For a leak value, include the code (a registry name, or its Hx and Hz) and the channel with its strength.

A code reported as protected when it leaks, or a leak reported smaller than the true one, belongs in
SECURITY.md rather than a public issue.
