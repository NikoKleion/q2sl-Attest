# Security

## Scope of reports

Two kinds of report belong here rather than in a public issue.

- A vulnerability in the code: for example, a file loader that executes content from the file it reads.
- A result that understates a leak: a code reported as protected when its syndrome depends on the
  logical state, or a leak reported smaller than the true one. A protected result may be relied on as a
  property of a code, so an error in that direction is treated as a security issue.

A result that overstates a leak, a crash, or a wrong number in a direction that cannot mislead someone
about protection can go in a public issue.

## How to report

Email nikomaxman@gmail.com with the subject line `q2sl attest security`. Include the version
(`python -c "import syndrome_leakage; print(syndrome_leakage.__version__)"`), the code (a registry name,
or its Hx and Hz), the channel and its strength, and the command or call that produced the result.

## Supported versions

| Version | Supported |
|---|---|
| 1.0.x | yes |
