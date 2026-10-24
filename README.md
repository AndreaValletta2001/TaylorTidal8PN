# TaylorTidal8PN

Frequency-domain TaylorF2 gravitational-wave models for quasi-circular,
aligned-spin compact binaries with **selected** tidal contributions through
8PN order. The Python package is named `tln_waveforms` and uses Bilby for
tidal-parameter conversion.

## Installation

Requires Python 3.10 or newer. In the repository root:

```bash
python -m pip install .
python -m unittest discover -s tests -v
python examples/generate_waveform.py
```

The example saves `waveform_example.npz` with frequency, plus, and cross arrays.
It does not require a detector PSD or an inference dataset.

## Example

```python
import numpy as np
from tln_waveforms import waveform_8pn

frequencies = np.linspace(20, 1000, 4096)
h = waveform_8pn(
    frequencies,
    chirp_mass=1.2, eta=0.24,
    chi_1=0.0, chi_2=0.0,
    lambda_tilde=400.0, delta_lambda_tilde=0.0,
    phase=0.0, luminosity_distance=100.0, theta_jn=0.5,
    reference_frequency=100.0,
    minimum_frequency=20.0, maximum_frequency=1000.0,
    magnetic_fluid="irrotational",
)
h_plus, h_cross = h["plus"], h["cross"]
```

The public interfaces use detector-frame chirp mass in solar masses,
luminosity distance in Mpc, frequency in Hz, angles in radians, and
aligned, dimensionless component spins. The low-level `pn_phase` instead
accepts total mass in kg.

## Available models

| Function or module | Included tidal contributions |
| --- | --- |
| `waveform_6p5pn` | Electric quadrupole at 5PN/6PN; magnetic quadrupole at 6PN; tidal tail and spin-tidal terms at 6.5PN |
| `waveform_7pn` | 6.5PN model plus the electric quadrupole and octupole terms at 7PN |
| `waveform_7pn_nooct` | 7PN model without the electric octupole |
| `waveform_8pn` | 7PN model plus selected quadratic and dynamical 8PN terms |
| `waveform_8pn_quadratic` | Quadratic 8PN term only, on top of the 7PN model |
| `waveform_8pn_dynamic` | Dynamical 8PN term only, on top of the 7PN model |

Additional named modules provide an 8PN model without dynamical running,
an agnostic quadratic 8PN coefficient, and an agnostic total 8PN coefficient
without running. These export `taylor_f2_APR_relative_binning` from their
respective modules in `tln_waveforms`.

## Conventions and limitations

- The magnetic response is irrotational by default; a static-fluid QUR is
  also available. These prescriptions should not be conflated with those
  for rotational Love numbers.
- `include_dynamic_log=False` suppresses the explicit 8PN logarithm only;
  the separate no-running model sets the running beta coefficients to zero
  also in the finite contribution.
- The waveform's `minimum_frequency` and `maximum_frequency` impose a hard
  analysis band. The public strain interface does not automatically choose
  a Roche cutoff or apply a taper; these must be supplied consistently by
  the calling analysis.

For conventions, equations, and references see
[phase coefficients](docs/COEFFICIENTS.md), [tidal relations](docs/QUR.md),
and [validation](docs/VALIDATION.md).

## Citation and license

Citation information and a distribution license will be finalized with the
authors before this repository is made public.
