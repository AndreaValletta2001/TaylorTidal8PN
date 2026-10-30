# Tidal relations and normalizations

Natural logarithms are used throughout. The default magnetic-fluid state
for binary inspirals is irrotational.

## Magnetic quadrupolar response

The dimensionless quantities are $`\Lambda_2=\lambda_2/m^5`$ and
$`\Sigma_2=\sigma_2/m^5`$ in geometrized units. For the fitted magnetic response,

```math
\ln\left|\Sigma_2\right| = \sum_{n=0}^{5} a_n \left(\ln\Lambda_2\right)^n ,
```

with coefficients:

| Fluid state | Sign of $`\Sigma_2`$ | $`a_0`$ | $`a_1`$ | $`a_2`$ | $`a_3`$ | $`a_4`$ | $`a_5`$ |
|:---|:---:|---:|---:|---:|---:|---:|---:|
| Irrotational (default) | Negative | -2.03 | 0.487 | 0.00969 | 0.00103 | -0.0000937 | 0.00000224 |
| Static | Positive | -2.66 | 0.786 | -0.01 | 0.00128 | -0.0000637 | 0.00000118 |

Source: [Jiménez Forteza et al., Eq. (22), Table I](https://arxiv.org/abs/1807.08016).
Zero electric deformability is mapped to zero magnetic deformability;
negative or nonfinite inputs are rejected.


## Electric octupole

The electric octupolar tidal deformability is determined from the quadrupolar
tidal deformability using the multipolar tidal Love quasi-universal relation
of Yagi (2014).

The relation is expressed as a fourth-order polynomial in the natural
logarithm of the dimensionless quadrupolar deformability:

```math
\ln \Lambda_3 =
-1.15
+ 1.18 \ln \Lambda_2
+ 0.0251 (\ln \Lambda_2)^2
- 0.00131 (\ln \Lambda_2)^3
+ 0.0000252 (\ln \Lambda_2)^4.
```

Equivalently, Lambda_3 is obtained by exponentiating the polynomial.

The relation is applied independently to each component of the binary.
For Lambda_2 = 0, the implementation sets Lambda_3 = 0, corresponding
to the black-hole limit. Negative tidal deformabilities are rejected.

Reference: Yagi (2014), multipole tidal Love relations,
arXiv:1311.0872.

## Spin-induced quadrupole

```math
\kappa = \exp\left( 0.194 + 0.0936\,z + 0.0474\,z^2 - 0.00421\,z^3 + 0.000123\,z^4 \right) ,
\qquad z = \ln\Lambda_2 .
```

For $`\Lambda_2=0`$ the implementation uses $`\kappa=1`$. The fit diverges as
$`\Lambda_2\to 0`$ ($`\kappa\sim 10^{10}`$ at $`\Lambda_2=10^{-6}`$), so for
$`0<\Lambda_2<1`$, below any neutron star, it is evaluated at $`\Lambda_2=1`$.
Source: [Yagi and Yunes, Eq. (54), Table I](https://arxiv.org/abs/1303.1528).

## Quadratic 8PN tidal response

The positive quantity $`q=-G\lambda_{222}/(Gm)^7`$ is modeled as

```math
q = \exp\left( 3.77 + 1.51\,z + 0.005\,z^2 \right) , \qquad z = \ln\frac{\Lambda_2}{4} .
```

The factor 4 and the sign in the definition of $`q`$ matter.
Source: [Pani et al.](https://arxiv.org/abs/2512.14663).

## Dynamical 8PN tidal response

The fit uses $`z=\ln\Lambda_2`$ and either of the two prescriptions:

```math
\begin{aligned}
\mu r_s = 1 :&\quad \ln D = 3.3840 + 1.1621\,z + 0.019180\,z^2 , \\
\mu r_s = 0.01 :&\quad \ln D = 3.9171 + 1.0707\,z + 0.023100\,z^2 , \\
&\quad \beta_2 = -\frac{64}{45} \left( 1 + \frac{321}{112}\Lambda_2 \right) .
\end{aligned}
```

The phase engines always use $`\mu r_s=1`$; the $`\mu r_s=0.01`$ fit is available
only through `dynamic_love_from_lambda0`.
Source: [Apostolidis et al.](https://arxiv.org/abs/2606.19446).
The no-running models set $`\beta_2=0`$ in both the finite and logarithmic parts;
`include_dynamic_log=False` removes only the explicit logarithm.

## Rotational tidal Love numbers at 6.5PN

Eight optional inputs are accepted:
`lambda1_23`, `sigma1_23`, `lambda1_32`, `sigma1_32` and their body-2
counterparts. They are zero by default; there is no automatic rotational QUR.
The inputs are dimensionless coefficients normalized by each star's mass:

```math
\texttt{lambda1\_23} = \frac{\lambda_{23}^{(1)}}{m_1^4} , \qquad
\texttt{lambda1\_32} = \frac{\lambda_{32}^{(1)}}{m_1^4} .
```

The corresponding `sigma` parameters have the same mass normalization.
Thus a term of the form $`\lambda_{23}^{(i)}/M^4`$ in the 6.5PN phase is evaluated
as $`X_i^4`$ times `lambda_i_23`, where $`X_i=m_i/M`$.
The implemented phase formula follows [Jiménez Forteza et al., Eq. (10)](https://arxiv.org/abs/1807.08016).
The static-fluid versus irrotational-fluid distinction for rotational
responses is not resolved by a QUR in this implementation.

## Other source-specific helpers

The APR mass–deformability approximation in the phase modules is an
equation-of-state-dependent helper and not an EoS-universal relation.
