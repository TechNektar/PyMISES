# Security Policy

## Supported Versions

We currently support the following versions of PyMISES with security updates:

| Version | Supported          |
| ------- | ------------------ |
| 0.1.x   | :white_check_mark: |

## Reporting a Vulnerability

If you discover a security vulnerability in PyMISES, please report it responsibly:

### Private Disclosure

1. **Do NOT** create a public GitHub issue for security vulnerabilities
2. Send an email to: info@technektar.com
3. Include the following information:
   - Description of the vulnerability
   - Steps to reproduce the issue
   - Potential impact
   - Suggested fix (if available)

### What to Expect

- **Response Time**: We aim to respond within 48 hours
- **Investigation**: We will investigate and assess the vulnerability
- **Fix Timeline**: Critical vulnerabilities will be addressed within 7 days
- **Disclosure**: We will coordinate with you on public disclosure timing

### Scope

Security issues in the following areas are in scope:
- Code execution vulnerabilities
- Data validation issues
- Authentication/authorization problems
- Dependency vulnerabilities

### Out of Scope

The following are generally out of scope:
- Issues in third-party dependencies (please report to those projects)
- Social engineering attacks
- Physical security issues

## Security Best Practices

When using PyMISES:

1. **Keep dependencies updated**: Regularly update NumPy, SciPy, and other dependencies
2. **Validate inputs**: Always validate user inputs, especially geometry files
3. **Sandbox execution**: Run PyMISES in isolated environments for untrusted data
4. **Monitor logs**: Check for unusual computational behavior or errors

## Acknowledgments

We appreciate the security research community and will acknowledge researchers who help improve PyMISES security (with their permission).

---

This security policy is effective as of May 26, 2024.
