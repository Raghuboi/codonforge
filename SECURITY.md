# Security Policy

CodonForge is a local command-line research tool. It does not run a server, process network requests, or execute untrusted plugins.

## Supported versions

| Version | Supported |
| --- | --- |
| 0.1.x | Yes |

## Reporting a vulnerability

Please open a private security advisory on GitHub or contact the maintainer through the repository profile.

Do not include secrets, private datasets, or unpublished biological data in public issues.

## Scope

Security-relevant issues include:

- unsafe file handling
- panic/crash on malformed input when an error should be returned
- dependency vulnerabilities
- incorrect handling of user-provided output paths
- accidental inclusion of secrets or private data in examples/tests

Scientific correctness issues should be filed as regular GitHub issues unless they create a security or safety risk.
