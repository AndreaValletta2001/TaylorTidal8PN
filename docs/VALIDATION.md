# Validation

Install the package and run its automated tests:

```bash
python -m pip install .
python -m unittest discover -s tests -v
python examples/generate_waveform.py
```

The unit tests check selected phase coefficients, spin tails, magnetic tidal
terms, the 7PN electric quadrupole, rotational tidal Love numbers (including
their zero default), input limits and public waveform interfaces.

These checks test implementation consistency. They do not independently
validate the accuracy of the phenomenological tidal relations outside their
calibration range, the completeness of an 8PN expansion, or waveform-based
Bayesian inference results. On GitHub, the included CI workflow additionally
runs the tests against multiple Python versions.
