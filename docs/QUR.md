# Tidal relations and normalizations

Natural logarithms are used throughout. The default magnetic-fluid state
for binary inspirals is irrotational.

## Magnetic quadrupolar response

The dimensionless quantities are `Lambda2=lambda2/m^5` and
`Sigma2=sigma2/m^5` in geometrized units. For the fitted magnetic response,
`ln(abs(Sigma2)) = sum_n a_n ln(Lambda2)^n`, with coefficients:

| Fluid state | Sign of Sigma2 | a0 | a1 | a2 | a3 | a4 | a5 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Irrotational (default) | negative | -2.03 | 0.487 | 0.00969 | 0.00103 | -0.0000937 | 0.00000224 |
| Static | positive | -2.66 | 0.786 | -0.01 | 0.00128 | -0.0000637 | 0.00000118 |

Source: [Jiménez Forteza et al., Eq. (22), Table I](https://arxiv.org/abs/1807.08016).
Zero electric deformability is mapped to zero magnetic deformability;
negative or nonfinite inputs are rejected. The fits need not remain accurate
at very high `Lambda2`, beyond their calibrated regime.

## Electric octupole

The model uses the power-law approximation

```text
Lambda3 = 0.23214267180912393 * Lambda2**1.345740158802329
```

Its input data were digitized and fitted with a power law; the data and
reproducible fitting script are provided in
[`examples/fit_octupole_qur.py`](../examples/fit_octupole_qur.py) (requires SciPy).
These numerical coefficients are **not** the quartic-polynomial coefficients
published in [Yagi (2014)](https://arxiv.org/abs/1311.0872), which is a
reference for multipolar Love relations. The function
`lambda3_over_lambda_fit` is a compatibility alias that returns `Lambda3`,
not the ratio `Lambda3/Lambda2`.

## Spin-induced quadrupole

```text
kappa = exp(0.194 + 0.0936*z + 0.0474*z**2 - 0.00421*z**3 + 0.000123*z**4)
z = ln(Lambda2)
```

For `Lambda2=0` the implementation uses `kappa=1`.
Source: [Yagi and Yunes, Eq. (54), Table I](https://arxiv.org/abs/1303.1528).

## Quadratic 8PN tidal response

The positive quantity `q=-G*lambda222/(G*m)^7` is modeled as

```text
z = ln(Lambda2/4)
q = exp(3.77 + 1.51*z + 0.005*z**2)
```

The factor `4` and the sign in the definition of `q` matter.
Source: [Pani et al.](https://arxiv.org/abs/2512.14663).

## Dynamical 8PN tidal response

The fit uses `z=ln(Lambda2)` and either of the two prescriptions:

```text
mu_rs=1:    ln(D) = 3.3840 + 1.1621*z + 0.019180*z**2
mu_rs=0.01: ln(D) = 3.9171 + 1.0707*z + 0.023100*z**2
beta2 = -(64/45) * (1 + 321*Lambda2/112)
```

The phase engine uses `mu_rs=1` unless changed explicitly.
Source: [Apostolidis et al.](https://arxiv.org/abs/2606.19446).
The no-running models set `beta2=0` in both the finite and logarithmic parts;
`include_dynamic_log=False` removes only the explicit logarithm.

## Rotational tidal Love numbers at 6.5PN

Eight optional inputs are accepted:
`lambda1_23`, `sigma1_23`, `lambda1_32`, `sigma1_32` and their body-2
counterparts. They are zero by default; there is no automatic rotational QUR.
The inputs are dimensionless coefficients normalized by each star's mass:

```text
lambda1_23 = lambda23^(1) / m1**4
lambda1_32 = lambda32^(1) / m1**4
```

The corresponding `sigma` parameters have the same mass normalization.
Thus a term of the form `lambda23^(i)/M**4` in the 6.5PN phase is evaluated
as `X_i**4 * lambda_i_23`, where `X_i=m_i/M`.
The implemented phase formula follows [Jiménez Forteza et al., Eq. (10)](https://arxiv.org/abs/1807.08016).
The static-fluid versus irrotational-fluid distinction for rotational
responses is not resolved by a QUR in this implementation.

## Other source-specific helpers

The APR mass–deformability approximation in the phase modules is an
equation-of-state-dependent helper and not an EoS-universal relation.
