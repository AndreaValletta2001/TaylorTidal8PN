# Phase coefficients and references

## Convention

Let $`X_i=m_i/M`$, $`\eta=X_1X_2`$, $`\delta=X_1-X_2\geq 0`$,
$`\chi_s=(\chi_1+\chi_2)/2`$, $`\chi_a=(\chi_1-\chi_2)/2`$, and
$`v=(\pi G M f/c^3)^{1/3}`$. The complex strain uses $`e^{i\phi}`$.
The implemented intrinsic phase is

```math
\begin{aligned}
\phi_{\rm PP} &= -\frac{\pi}{4} - \frac{3}{128\eta}\, v^{-5} \sum_n \left(p_n+s_n\right) v^n , \\
\phi_{\rm T} &= \frac{3}{128\eta} \left[ T_{10}\, v^5 + T_{12}\, v^7 + T_{13}\, v^8 + T_{16}\, v^{11} \right] + C_7\, v^9 , \\
\phi &= \phi_{\rm PP} - \phi_{\rm T} .
\end{aligned}
```

The public strain additionally subtracts the phase at the reference frequency
and applies its orbital-phase offset. With this phase convention, a negative
leading tidal coefficient produces a positive contribution to $`\phi`$. The
relative sign follows from the adopted Fourier convention.
`pn_phase_order` is the integer index $`n`$ (twice the PN order); its default -1
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
With $`a_A=X_1\chi_1`$, $`a_B=X_2\chi_2`$ and
$`a_Q^2=\kappa_1 a_A^2+2a_Aa_B+\kappa_2 a_B^2`$:

```math
\begin{aligned}
s_7^{\rm SS,\,tail} &= -\pi \left[ 400\, a_Q^2 + \frac{15}{2} \left(a_A-a_B\right)^2 \right] , \\
s_8^{\rm SO,\,tail} &= \pi \left(1-3\ln v\right) \left[ \left( \frac{233915}{168} - \frac{99185}{252}\eta \right) \delta\chi_a
  + \left( \frac{233915}{168} - \frac{3970375}{2268}\eta + \frac{19655}{189}\eta^2 \right) \chi_s \right] .
\end{aligned}
```

These are [Nagar et al., Eqs. (47)–(48)](https://arxiv.org/pdf/1812.07923);
the second expression is equivalently Mishra et al., Eq. (6c). The source uses
the equivalent mass-weighted spin basis for the latter. Tests independently
isolate both additions and verify the default includes the 4PN spin-orbit tail.

## Linear tidal sectors

The following are coefficients inside the common $`3/(128\eta)`$ convention:

```math
\begin{aligned}
T_{10} &= -\frac{39}{2}\tilde\Lambda , \\
T_{12} &= -\frac{3115}{64}\tilde\Lambda + \frac{6595}{364}\,\delta\,\delta\tilde\Lambda + \tilde\Sigma , \\
T_{13} &= \hat\Lambda + \hat\Sigma + \hat\Gamma + \frac{39\pi}{2}\tilde\Lambda .
\end{aligned}
```

For each body $`i`$ and its companion $`j`$, sum the following contributions:

```math
\begin{aligned}
\tilde\Sigma_i &= \left[ \frac{6920}{7} - \frac{20740}{21X_i} \right] X_i^5\, \Sigma_i , \\
\hat\Lambda_i &= \left\{ \left[ \frac{593}{4} - \frac{1105}{8X_i} + \frac{567}{8}X_i - 81X_i^2 \right] \chi_j
  + \left[ -\frac{6607}{8} + \frac{6639}{8}X_i - 81X_i^2 \right] \chi_i \right\} X_i^5\, \Lambda_i , \\
\hat\Sigma_i &= \left\{ \left[ -\frac{9865}{3} + \frac{4933}{3X_i} + 1644X_i \right] \chi_j - \chi_i \right\} X_i^5\, \Sigma_i .
\end{aligned}
```

These magnetic and spin-tidal coefficients follow
[Jiménez Forteza et al., Eqs. (6)–(8)](https://arxiv.org/pdf/1807.08016).
Both magnetic terms use the same signed $`\Sigma`$ QUR; its full coefficient table
is in [QUR.md](QUR.md). Rotational tidal Love numbers are zero by default.
The 6.5PN $`\hat\Gamma`$ expression follows Jimenez Forteza et al., Eq. (10).
Optional RTLNs are normalized by $`m_i^4`$; see QUR.md.

The implemented 7PN electric-sector direct phase coefficient is the mass
quadrupole and mass octupole part of the SPA phase in
[Henry, Faye & Blanchet, Sec. VI](https://arxiv.org/abs/2005.13367), including
the corrections of their Erratum (arXiv v4). With
$`\mu_\pm = \left(X_2X_1^4\Lambda_1 \pm X_1X_2^4\Lambda_2\right)/2`$:

```math
\begin{aligned}
C_7 = -\frac{9}{16\eta^2} \Bigg[ &\left( \frac{136190135}{27433728} + \frac{978554825}{4572288}\eta - \frac{281935}{6048}\eta^2 + \frac{5}{3}\eta^3 \right) \mu_+ \\
 &+ \delta \left( \frac{136190135}{27433728} + \frac{213905}{2592}\eta + \frac{1585}{1296}\eta^2 \right) \mu_- \Bigg]
 - \frac{125}{12} \left( X_1^5\Lambda_{3,1} + X_2^5\Lambda_{3,2} \right) .
\end{aligned}
```

`coeffLT` and `coeffDLT` return the equivalent coefficients of `lambda_tilde`
and `delta_lambda_tilde`. The octupole term is omitted in the no-octupole model.
The 7PN magnetic term and the 7.5PN tails of the same reference are not
included. These are already coefficients of $`v^9`$: no additional
$`3/(128\eta)`$ is applied to $`C_7`$.

## Selected 8PN sectors

Define the positive dimensionless magnitude
$`q_i=-G\lambda_{222}^{(i)}/(Gm_i)^7`$ using the fit in QUR.md. Its input is
$`\ln(\Lambda_i/4)`$, not $`\ln\Lambda_i`$.

```math
\begin{aligned}
\rho_\pm &= 2\left( q_1 X_1^6 X_2 \pm q_2 X_2^6 X_1 \right) , \\
Q_8 &= 19\rho_+ - 17\,\delta\,\rho_- , \\
T_{16}^{\rm quad} &= -\frac{45}{11}\, Q_8 .
\end{aligned}
```

The positive fitted magnitude and the explicit minus sign correspond to
[Pani et al., Eqs. (77), (85), (96), footnote 13](https://arxiv.org/html/2512.14663v2).
For dynamics define $`D_i=\Lambda^{\rm dyn}_i`$ from the $`\mu r_s=1`$ fit in QUR.md,
$`\beta_i=-\frac{64}{45}\left(1+\frac{321}{112}\Lambda_i\right)`$, and $`x=v^2`$:

```math
\begin{aligned}
F_8 &= \sum_i X_i^5 \left[ \left(8X_i+147\eta\right) D_i + \left( \frac{38}{11}X_i + \frac{1253}{44}\eta \right) \beta_i \right] , \\
B_8^{(i)} &= X_i^5 \left(8X_i+147\eta\right) \beta_i , \\
T_{16}^{\rm dyn} &= -\frac{15}{11} \left[ F_8 + \sum_i B_8^{(i)} \ln\frac{1}{2\gamma x X_i} \right] .
\end{aligned}
```

Here $`\gamma`$ is `dynamic_log_gamma`, which defaults to 1. These are
[Apostolidis et al., Eqs. (60)–(61), (66), (70)–(71)](https://arxiv.org/html/2606.19446v2).
The no-running models set $`\beta_i=0`$ in both finite and logarithmic terms.
`include_dynamic_log=False` in the main engine removes only the logarithm.
The default selected 8PN coefficient is the sum of the quadratic and dynamical
parts. The free `lambda_8pn_quadratic` replaces $`Q_8`$; the free
`lambda_8pn_total` replaces $`3Q_8+F_8(\beta=0)`$, with coefficient $`-15/11`$.
