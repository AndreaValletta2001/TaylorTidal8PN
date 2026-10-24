# Phase coefficients and references

## Convention

Let `X_i=m_i/M`, `eta=X_1*X_2`, `delta=X_1-X_2>=0`,
`chi_s=(chi_1+chi_2)/2`, `chi_a=(chi_1-chi_2)/2`, and
`v=(pi*G*M*f/c^3)^(1/3)`. The complex strain uses `exp(i*phi)`.
The implemented intrinsic phase is

```text
phi_PP = -pi/4 - 3/(128*eta) * v^-5 * sum_n[(p_n+s_n)*v^n]
phi_T  = 3/(128*eta) * [T10*v^5 + T12*v^7 + T13*v^8 + T16*v^11]
         + C7*v^9
phi = phi_PP - phi_T
```

The public strain additionally subtracts the phase at the reference frequency
and applies its orbital-phase offset. With this phase convention, a negative
leading tidal coefficient produces a positive contribution to `phi`. The
relative sign follows from the adopted Fourier convention.
`pn_phase_order` is the integer index n (twice the PN order); its default -1
includes all implemented point-particle/spin terms. `tidal_pn_order` is in PN
units and selects the 6.5, 7 or 8PN models.

## Point-particle and spin sectors

The source retains the aligned-spin orbit, spin-spin,
spin-induced quadrupole and cubic-spin coefficients.
The cubic-spin terms retain their inherited black-hole-like prescription;
no independent neutron-star spin-induced octupole QUR is introduced.
Reference expressions are in [Mishra et al., Eq. (6)](https://arxiv.org/pdf/1601.05588)
and the [LALSuite PN coefficient source](https://lscsoft.docs.ligo.org/lalsuite/dev/lalsimulation/_l_a_l_sim_inspiral_p_n_coefficients_8c_source.html).
The spin-induced quadrupole relation is documented in [QUR.md](QUR.md).

The following spin-tail terms are implemented.
With `a_A=X_1*chi_1`, `a_B=X_2*chi_2`,
`a_Q^2=kappa_1*a_A^2+2*a_A*a_B+kappa_2*a_B^2`:

```text
s7_SS_tail = -pi * [400*a_Q^2 + (15/2)*(a_A-a_B)^2]

s8_SO_tail = pi*(1-3*ln(v)) * {
    (233915/168 - 99185*eta/252)*delta*chi_a
  + (233915/168 - 3970375*eta/2268 + 19655*eta^2/189)*chi_s
}
```

These are [Nagar et al., Eqs. (47)–(48)](https://arxiv.org/pdf/1812.07923);
the second expression is equivalently Mishra et al., Eq. (6c). The source uses
the equivalent mass-weighted spin basis for the latter. Tests independently
isolate both additions and verify the default includes the 4PN spin-orbit tail.

## Linear tidal sectors

The following are coefficients inside the common `3/(128*eta)` convention:

```text
T10 = -(39/2)*lambda_tilde
T12 = -(3115/64)*lambda_tilde
      + (6595/364)*delta*delta_lambda_tilde + sigma_tilde
T13 = lambda_hat + sigma_hat + gamma_hat + (39*pi/2)*lambda_tilde
```

For each body i and its companion j, sum the following contributions:

```text
sigma_tilde_i = [6920/7 - 20740/(21*X_i)] * X_i^5 * Sigma_i
lambda_hat_i = {
    [593/4 - 1105/(8*X_i) + 567*X_i/8 - 81*X_i^2]*chi_j
  + [-6607/8 + 6639*X_i/8 - 81*X_i^2]*chi_i
} * X_i^5 * Lambda_i
sigma_hat_i = {
    [-9865/3 + 4933/(3*X_i) + 1644*X_i]*chi_j - chi_i
} * X_i^5 * Sigma_i
```

These magnetic and spin-tidal coefficients follow
[Jiménez Forteza et al., Eqs. (6)–(8)](https://arxiv.org/pdf/1807.08016).
Both magnetic terms use the same signed Sigma QUR; its full coefficient table
is in [QUR.md](QUR.md). Rotational tidal Love numbers are zero by default.
The 6.5PN gamma_hat expression follows Jimenez Forteza et al., Eq. (10).
Optional RTLNs are normalized by m_i**4; see QUR.md.

The implemented 7PN electric-sector direct phase
coefficient is

```text
C7 = coeffLT(m1,m2)*lambda_tilde + coeffDLT(m1,m2)*delta_lambda_tilde
     - (125/12)*(X_1^5*Lambda3_1 + X_2^5*Lambda3_2)
```

The exact rational mass functions are included in each engine under `coeffLT`
and `coeffDLT`. The final line is omitted in the no-octupole model. These are
already coefficients of `v^9`: no additional `3/(128*eta)` is applied to C7.

## Selected 8PN sectors

Define the positive dimensionless magnitude
`q_i=-G*lambda222_i/(G*m_i)^7` using the fit in QUR.md. Its input is
`ln(Lambda_i/4)`, not `ln(Lambda_i)`.

```text
rho_plus  = 2*(q_1*X_1^6*X_2 + q_2*X_2^6*X_1)
rho_minus = 2*(q_1*X_1^6*X_2 - q_2*X_2^6*X_1)
Q8 = 19*rho_plus - 17*delta*rho_minus
T16_quadratic = -(45/11)*Q8
```

The positive fitted magnitude and the explicit minus sign correspond to
[Pani et al., Eqs. (77), (85), (96), footnote 13](https://arxiv.org/html/2512.14663v2).
For dynamics define `D_i=Lambda_dyn_i` from the mu_rs=1 fit in QUR.md,
`beta_i=-(64/45)*(1+321*Lambda_i/112)`, and `x=v^2`:

```text
F8 = sum_i X_i^5 * [(8*X_i+147*eta)*D_i
                   + (38*X_i/11+1253*eta/44)*beta_i]
B8_i = X_i^5 * (8*X_i+147*eta)*beta_i
T16_dynamic = -(15/11) * [F8 + sum_i B8_i*ln(1/(2*gamma*x*X_i))]
```

Here `gamma=dynamic_log_gamma` defaults to 1. These are
[Apostolidis et al., Eqs. (60)–(61), (66), (70)–(71)](https://arxiv.org/html/2606.19446v2).
The no-running models set beta_i=0 in both finite and logarithmic terms.
`include_dynamic_log=False` in the main engine removes only the logarithm.
The default selected 8PN coefficient is the sum of the quadratic and dynamical
parts. The free `lambda_8pn_quadratic` replaces Q8; the free
`lambda_8pn_total` replaces `3*Q8+F8(beta=0)`, with coefficient `-(15/11)`.

The model name `complete` refers to these selected sectors. It does not assert
PN completeness: for example, there is no 7.5PN tidal-tail term. Empirical QURs
are extrapolated in the large-Lambda subsolar regime, as discussed in QUR.md.
