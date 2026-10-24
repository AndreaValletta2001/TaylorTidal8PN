"""TaylorF2 with selected tidal sectors through 8PN.

The complex strain uses exp(i*phi), with phi = phi_PP - phi_T.
The 7PN electric coefficients are direct v^9 phase coefficients, so their
common 3/(128*eta) prefactor is compensated in the phase assembly.
Masses supplied to pn_phase are in kg; public waveform masses are in solar masses.
The hard analysis band is applied using minimum_frequency/maximum_frequency.
"""
#Python libraries
import numpy as np
from bilby.gw.conversion import lambda_tilde_delta_lambda_tilde_to_lambda_1_lambda_2 as lambda_conv
from .magnetic import sigma_from_lambda_fit
from .octupole import lambda3_from_lambda_fit
#------------------ Constants ------------------------------------------

##Euler constant 
gamma_e = 0.5772156649
## solar mass
solar_m = 1.9884099021470415E30 #Kg
## Megaparsec
m_parsec = 3.085677581491367E22 #m
## speed of light
light_s = 299792458. #m/s
## Newton's constant
g_newton = 6.67430E-11 #m^3/kg/s
#EOS constant
C = 7.3e4*0.5**4.7

########################### FUNCTION DEFINITIONS ###############################
########################### FUNCTION DEFINITIONS ###############################

# Roche frequency
def f_RO(mass_1,mass_2):
    f_ro = -26.9 - 35.5*mass_1 - 3.02*mass_1**2 + 1690*mass_2 - 575*mass_2**2
    return f_ro

# Tapering with a generic cutoff frequency f_cut.
def tapering(f, M, f_RO, delta_gamma):
    # Schwarzschild-ISCO gravitational-wave frequency. M is in kg.
    fisco = light_s**3 / (np.pi * 6.0**1.5 * g_newton * M)
    gamma = f_RO / fisco

    inner_frac = (1.0 + np.exp(-gamma / delta_gamma)) / (
        1.0
        + np.exp(
            (np.asarray(f) / fisco - gamma) / delta_gamma
        )
    )

    return inner_frac

def compute_lambda_tilde(eta,lambda_s,lambda_a):
    #5pN leading tidal coefficient from http://www.arxiv.org/abs/0911.3535 
    delta = np.sqrt(1-4*eta)
    return 16/13*((1+7*eta-31*eta**2)*lambda_s+\
                  delta*(1+9*eta-11*eta**2)*lambda_a)

def compute_delta_lambda_tilde(eta,lambda_s,lambda_a):
    #6pN tidal coefficient from https://arxiv.org/abs/1410.8866
    delta = np.sqrt(1-4*eta)
    return delta*(1-13272/1319*eta+8944/1319*eta**2)*lambda_s+\
        (1-15910/1319*eta+32850/1319*eta**2+3380/1319*eta**3)*lambda_a

def compute_sigma_tilde(eta1, eta2,sigma_1,sigma_2):
    return (6920/7-20740/(21*eta1))*eta1**5*sigma_1+(6920/7-20740/(21*eta2))*eta2**5*sigma_2

def compute_lambda_hat(eta1,eta2,chi_1,chi_2,lambda_1,lambda_2):
    t_a = (593/4-1105/(8*eta1)+567*eta1/8-81*eta1**2)*chi_2
    t_b = (-6607/8+6639*eta1/8-81*eta1**2)*chi_1
    t_c = (593/4-1105/(8*eta2)+567*eta2/8-81*eta2**2)*chi_1
    t_d = (-6607/8+6639*eta2/8-81*eta2**2)*chi_2
    return (t_a+t_b)*eta1**5*lambda_1+(t_c+t_d)*eta2**5*lambda_2

def compute_sigma_hat(eta1,eta2,chi_1,chi_2,sigma_1,sigma_2):
    t_a = (-9865/3+4933/(3*eta1)+1644*eta1)*chi_2-chi_1
    t_b = (-9865/3+4933/(3*eta2)+1644*eta2)*chi_1-chi_2
    return t_a*eta1**5*sigma_1+t_b*eta2**5*sigma_2

def compute_kappa_hat(lambda_tilde):
    return (39*np.pi/2)*lambda_tilde

def compute_gamma_hat(chi_1, chi_2, eta1, eta2, eta, m_tot,
                      lambda1_23, sigma1_23, lambda1_32, sigma1_32,
                      lambda2_23, sigma2_23, lambda2_32, sigma2_32):
    """Rotational tidal contribution at 6.5PN.

    Reference: Jimenez Forteza et al. (2018), Eq. (10).

    The eight RTLN inputs are dimensionless and normalized by
    the fourth power of their respective component masses (G=c=1).

    The factors X_i**4 convert the component-mass normalization
    into the total-mass normalization of the published expression.

    All RTLNs default to zero. No rotational QUR is assumed.
    """
    values = (
        lambda1_23, sigma1_23, lambda1_32, sigma1_32,
        lambda2_23, sigma2_23, lambda2_32, sigma2_32
    )
    if not all(np.isfinite(value) for value in values):
        raise ValueError("Rotational tidal Love numbers must be finite")

    X1, X2 = eta1, eta2

    gamma1 = chi_1 * X1**4 * (
        (856*X1 - 816*X1**2) * lambda1_23
        - (833*X1/3 - 278*X1**2) * sigma1_23
        - eta * (272*lambda1_32 - 204*sigma1_32)
    )

    gamma2 = chi_2 * X2**4 * (
        (856*X2 - 816*X2**2) * lambda2_23
        - (833*X2/3 - 278*X2**2) * sigma2_23
        - eta * (272*lambda2_32 - 204*sigma2_32)
    )

    return gamma1 + gamma2


def mass_components(chirp_mass, eta):
    delta = np.sqrt(1-4*eta)
    mass_1=chirp_mass/(2*eta**(3/5))*(1+delta)
    mass_2=chirp_mass/(2*eta**(3/5))*(1-delta)
    return mass_1, mass_2

def lambda_from_mass_APR(mass):
    
    #fit from mathematica
    a = -5.095773636110624
    b = -0.4072809551261860
    c = 3.258321049591981
    
    #tide
    tide=10**(a * np.log10(mass) + b * (np.log10(mass))**2 + c)
    
    return tide

#fit for light NSs _> extract const assuming tilt=4.7
def log_const_fit(lambda_i, m):
    
    #see eq. (11) of https://arxiv.org/abs/2402.18656
    return np.log10(lambda_i*(m/0.5)/(7.3*10**4))

#fit for light NSs _> extract constand tilt given masses and tides
def const_tilt_fit(lambda_1, lambda_2, mass_1, mass_2):
    
    #EOS tilt
    tilt = (np.log10(lambda_1/lambda_2)) / (np.log10(mass_2/mass_1))

    #EOS constant
    const = np.sqrt(lambda_1 * lambda_2 * (mass_1 * mass_2 / (0.5)**2)**tilt) / 7.3e4

    return const, tilt

#I-Love-Q
def quadrupole_from_loveq(lambda_i):
    if lambda_i!=0:

        # Quadrupoles for a NS. The fit diverges as Lambda -> 0 (kappa ~ 1e10
        # at Lambda = 1e-6), so it is evaluated at Lambda >= 1, below any NS.
        z = np.log(max(lambda_i, 1.0))
        kappa=np.exp(0.194+0.0936*z+0.0474*z**2-0.00421*z**3+0.000123*z**4)
    else: 
        # Quadrupoles for a Kerr BH
        kappa=1
    return kappa

def compute_electric_quadrupole_7pn(mass_1, mass_2, lambda_1, lambda_2):
    """7PN (NNL) electric-quadrupole tidal coefficient of v^9 in the phase.

    Henry, Faye & Blanchet, Phys. Rev. D 102, 044033 (2020), SPA phase of
    Sec. VI, including the corrections of the Erratum (arXiv:2005.13367v4).
    With mu_A = m_A^5 Lambda_A, mu_pm = (m2/m1 mu_1 +- m1/m2 mu_2)/(2 M^5).
    This is a direct phase coefficient: no 3/(128*eta) prefactor is applied.
    """
    M = mass_1 + mass_2
    X1 = mass_1 / M
    X2 = mass_2 / M
    nu = X1 * X2
    Delta = X1 - X2

    mu_plus = 0.5 * (X2 * X1**4 * lambda_1 + X1 * X2**4 * lambda_2)
    mu_minus = 0.5 * (X2 * X1**4 * lambda_1 - X1 * X2**4 * lambda_2)

    c_plus = (136190135/27433728 + 978554825/4572288*nu
              - 281935/6048*nu**2 + 5/3*nu**3)
    c_minus = 136190135/27433728 + 213905/2592*nu + 1585/1296*nu**2

    return -9 / (16 * nu**2) * (c_plus * mu_plus + Delta * c_minus * mu_minus)

def coeffLT(m1,m2):
    """Coefficient of lambda_tilde in compute_electric_quadrupole_7pn."""
    return compute_electric_quadrupole_7pn(m1, m2, *lambda_conv(1.0, 0.0, m1, m2))

def coeffDLT(m1,m2):
    """Coefficient of delta_lambda_tilde in compute_electric_quadrupole_7pn."""
    return compute_electric_quadrupole_7pn(m1, m2, *lambda_conv(0.0, 1.0, m1, m2))



# ==================================
# Static quadratic Love numbers
# ==================================
# These functions implement the 8PN static quadratic-Love coefficient
# internally, keeping the public waveform signature identical to TaylorF2_APR_7PN.
# The sampled tidal variables remain only lambda_tilde and delta_lambda_tilde.

def lambda222bar_from_lambda2bar(lambda2bar):
    """
    Quasi-universal relation for the static quadratic Love number.

        Y = 3.77 + 1.51 X + 0.005 X^2
        X = log[ G lambda_2_EFT / (G m)^5 ] = log(Lambda_2/4)
        Y = log[ -G lambda_222 / (G m)^7 ]

    Natural logarithms are assumed.  The input is the usual dimensionless
    quadrupolar tidal deformability Lambda_2. The return value is the
    positive magnitude -G*lambda_222/(G*m)^7; see Pani et al., Eq. (77).
    """
    lambda2bar = np.asarray(lambda2bar, dtype=float)
    good = np.isfinite(lambda2bar) & (lambda2bar > 0.0)

    out = np.zeros_like(lambda2bar, dtype=float)
    if not np.any(good):
        return out.item() if out.shape == () else out

    X = np.log(lambda2bar[good] / 4.0)
    Y = 3.77 + 1.51 * X + 0.005 * X**2
    out[good] = np.exp(Y)
    return out.item() if out.shape == () else out


def compute_quadratic_love_rho_tilde(mass_1, mass_2, lambda222_1, lambda222_2):
    """
    Effective quadratic-Love combinations.

    lambda222_i stores the positive magnitude -G*lambda_222^(i)/(G*m_i)^7.
    The returned combinations use these magnitudes; the phase coefficient
    below supplies the physical minus sign (Pani et al., Eqs. 85 and 96).
    """
    M = mass_1 + mass_2
    X1 = mass_1 / M
    X2 = mass_2 / M

    rho_tilde_plus = 2.0 * (lambda222_1 * X1**6 * X2 + lambda222_2 * X2**6 * X1)
    rho_tilde_minus = 2.0 * (lambda222_1 * X1**6 * X2 - lambda222_2 * X2**6 * X1)

    return rho_tilde_plus, rho_tilde_minus


def compute_quadratic_love_8PN_term(mass_1, mass_2, delta, lambda222_1, lambda222_2):
    """
    Coefficient to be added to tidal_phase_terms[4].

    Convention used here:
        phi_T = 3/(128 eta) v^(-5) [ ... + coeff_8PN v^16 ].
    """
    rho_tilde_plus, rho_tilde_minus = compute_quadratic_love_rho_tilde(
        mass_1,
        mass_2,
        lambda222_1,
        lambda222_2,
    )

    return -(45.0 / 11.0) * (19.0 * rho_tilde_plus - 17.0 * delta * rho_tilde_minus)


# ==================================
# Dynamical quadrupolar Love numbers
# ==================================
# These functions implement the 8PN dynamical-tide coefficient internally,
# keeping the public waveform signature identical to TaylorF2_APR_7PN.
# The sampled tidal variables remain only lambda_tilde and delta_lambda_tilde.

def dynamic_love_from_lambda0(lambda0, mu_rs=1.0):
    """
    Quasi-universal relation for the dynamical Love coefficient.

        X = log Lambda_0
        Y = log Lambda_2^dyn

    Fits:
        mu_rs = 1.0:  Y = 3.3840 + 1.1621 X + 0.019180 X^2
        mu_rs = 0.01: Y = 3.9171 + 1.0707 X + 0.023100 X^2

    Natural logarithms are assumed.  Lambda_0 is the usual dimensionless
    static quadrupolar tidal deformability.
    """
    lambda0 = np.asarray(lambda0, dtype=float)
    good = np.isfinite(lambda0) & (lambda0 > 0.0)

    out = np.zeros_like(lambda0, dtype=float)
    if not np.any(good):
        return out.item() if out.shape == () else out

    X = np.log(lambda0[good])
    if np.isclose(mu_rs, 1.0):
        Y = 3.3840 + 1.1621 * X + 0.019180 * X**2
    elif np.isclose(mu_rs, 0.01):
        Y = 3.9171 + 1.0707 * X + 0.023100 * X**2
    else:
        raise ValueError("Available dynamical-Love fits are mu_rs=1.0 and mu_rs=0.01.")

    out[good] = np.exp(Y)
    return out.item() if out.shape == () else out


def beta2_dynamic_running(lambda0):
    """
    Universal logarithmic running coefficient beta_2.
    lambda0 is the usual dimensionless static quadrupolar Lambda.
    """
    return -(64.0 / 45.0) * (1.0 + (321.0 / 112.0) * lambda0)


def compute_dynamic_love_8PN_terms(v, mass_1, mass_2, eta, lambda0_1, lambda0_2,
                                   lambda_dyn_1, lambda_dyn_2,
                                   include_dynamic_log=True,
                                   dynamic_log_gamma=1.0):
    """
    Coefficient to be added to tidal_phase_terms[4].

    Convention used here:
        phi_T = 3/(128 eta) v^(-5) [ ... + coeff_8PN v^16 ].

    The sign is chosen to match the existing tidal convention in this file,
    where the leading 5PN term is -39/2 * lambda_tilde.
    """
    M = mass_1 + mass_2
    X1 = mass_1 / M
    X2 = mass_2 / M

    beta2_1 = beta2_dynamic_running(lambda0_1)
    beta2_2 = beta2_dynamic_running(lambda0_2)

    Lambda8PN_finite = (
        X1**5 * ((8.0 * X1 + 147.0 * eta) * lambda_dyn_1
                 + ((38.0 / 11.0) * X1 + (1253.0 / 44.0) * eta) * beta2_1)
        + X2**5 * ((8.0 * X2 + 147.0 * eta) * lambda_dyn_2
                   + ((38.0 / 11.0) * X2 + (1253.0 / 44.0) * eta) * beta2_2)
    )

    coeff_finite = -(15.0 / 11.0) * Lambda8PN_finite

    if not include_dynamic_log:
        return coeff_finite

    B8PN_1 = X1**5 * (8.0 * X1 + 147.0 * eta) * beta2_1
    B8PN_2 = X2**5 * (8.0 * X2 + 147.0 * eta) * beta2_2

    x = np.maximum(v**2, 1e-300)
    gamma_run = float(dynamic_log_gamma)
    if gamma_run <= 0.0:
        raise ValueError("dynamic_log_gamma must be positive.")

    # r_s^(i) = 2 G m_i, hence GM/(x gamma r_s_i) = 1/(2 gamma x X_i).
    log_arg_1 = 1.0 / (2.0 * gamma_run * x * X1)
    log_arg_2 = 1.0 / (2.0 * gamma_run * x * X2)

    coeff_log = -(15.0 / 11.0) * (
        B8PN_1 * np.log(log_arg_1) + B8PN_2 * np.log(log_arg_2)
    )

    return coeff_finite + coeff_log


########################### TAYLOR-F2 PHASE #############################

# pn_phase_order=-1 includes the implemented point-particle/spin terms.
# Nonspinning terms stop at 3.5PN; the spin-orbit tail reaches 4PN.
# tidal_pn_order separately selects the implemented tidal sectors.
@np.errstate(divide='ignore', invalid='ignore')
def pn_phase(
    frequency_array,
    total_mass,
    eta,
    chi_s,
    chi_a,
    lambda_tilde,
    delta_lambda_tilde,
    lambda_8pn_quadratic=None,
    lambda1_23=0.0,
    sigma1_23=0.0,
    lambda1_32=0.0,
    sigma1_32=0.0,
    lambda2_23=0.0,
    sigma2_23=0.0,
    lambda2_32=0.0,
    sigma2_32=0.0,
    pn_phase_order=-1,
    include_dynamic_log=True,
    dynamic_log_gamma=1.0,
    tidal_pn_order=8.0,
    include_octupole_7pn=True,
    include_quadratic_8pn=True,
    include_dynamic_8pn=True,
    magnetic_fluid="irrotational",
    include_magnetic_6pn=True,
    include_magnetic_spin_tidal=True,
):

    delta = np.sqrt(1-4*eta)
    v = (g_newton/(light_s**3)*np.pi*frequency_array*total_mass)**(1/3.)
    
    # Initialize coefficient arrays.    
    pp_phase_terms = np.zeros((9, len(v)))
    spin_phase_terms = np.zeros((9, len(v)))
    tidal_phase_terms = np.zeros((5, len(v)))
    
    # Masses computation
    M_ad = total_mass/solar_m
    mass_1, mass_2 = mass_components(M_ad*eta**(3/5),eta)

    # Tides computation
    lambda_1, lambda_2 = lambda_conv(lambda_tilde, delta_lambda_tilde, mass_1, mass_2)

    # Magnetic tidal coefficients reconstructed through quasi-universal relations.
    # QUR magnetic TLNs enter the 6PN term and the spin-tidal 6.5PN term.
    sigma_1 = sigma_from_lambda_fit(lambda_1, fluid=magnetic_fluid)
    sigma_2 = sigma_from_lambda_fit(lambda_2, fluid=magnetic_fluid)

    # Static quadratic Love numbers and dynamical Love coefficients derived
    # internally from the static quadrupolar TLNs.
    lambda222_1 = lambda222bar_from_lambda2bar(lambda_1)
    lambda222_2 = lambda222bar_from_lambda2bar(lambda_2)
    lambda_dyn_1 = dynamic_love_from_lambda0(lambda_1, mu_rs=1.0)
    lambda_dyn_2 = dynamic_love_from_lambda0(lambda_2, mu_rs=1.0)

    L31 = lambda3_from_lambda_fit(lambda_1)
    L32 = lambda3_from_lambda_fit(lambda_2)

    X1 = mass_1 / (mass_1 + mass_2)
    X2 = mass_2 / (mass_1 + mass_2)

    # Mass-fraction notation used by the 6.5PN spin-tidal sector.
    eta1 = X1
    eta2 = X2
    m_tot = mass_1 + mass_2

    # I-Love-Q computation
    kappa_1 = quadrupole_from_loveq(lambda_1)
    kappa_2 = quadrupole_from_loveq(lambda_2)
    kappa_s = (kappa_1+kappa_2)/2
    kappa_a = (kappa_1-kappa_2)/2
    #kappa_s=1.0
    #kappa_a=0.0
    

    #Non-zero point particle phase coefficients 
    pp_phase_terms[0] = 1.
    pp_phase_terms[2] = 3715/756.+55*eta/9.
    pp_phase_terms[3] = -16*np.pi     
    pp_phase_terms[4] = 15293365/508032.+27145*eta/504.+3085*eta**2/72. 
    pp_phase_terms[5] = (1.+np.log(v**3)) *\
             ((38645*np.pi/756.)-(65*np.pi*eta/9.)) 
    pp_phase_terms[6] = 11583231236531/4694215680.-6848*gamma_e/21.-640*np.pi**2/3.+(-15737765635/3048192. +
             2255*np.pi**2/12.)*eta+76055*eta**2/1728.-127825*eta**3/1296.-6848/63. *\
             np.log(64*v**3)
    pp_phase_terms[7] = 77096675*np.pi/254016.+378515*np.pi*eta/1512.-74045*np.pi*eta**2/756.
    
    #Non-zero spin phase coefficients 
    spin_phase_terms[3] = 113/3.*delta*chi_a+(113/3.-76*eta/3.)*chi_s
    
    spin_phase_terms[4] = -5/8.*chi_s**2*(1+156*eta+80*delta*kappa_a +
             80*(1-2*eta)*kappa_s)+chi_a**2*(-5/8.-50*delta*kappa_a-50*kappa_s+100*eta*(1+kappa_s)) -\
             5/4.*chi_a*chi_s*(delta+80*(1-2*eta)*kappa_a+80*delta*kappa_s)
    spin_phase_terms[5] = (1.+np.log(v**3))*\
             (delta*(-732985/2268.-(140*eta/9.))*chi_a + (-732985/2268.+24260*eta/81.+340*eta**2/9.)*chi_s)
    spin_phase_terms[6] = np.pi*(2270/3.*delta*chi_a+(2270/3.-520*eta)*chi_s)+chi_s**2*((-1344475/2016.) +
             829705/504*eta+3415/9.*eta**2+delta*(26015/28.-1495/6.*eta)*kappa_a +
             (26015/28.-44255/21.*eta-240*eta**2)*kappa_s)+chi_a**2*(-1344475/2016.+267815/252.*eta -
             240*eta**2+delta*(26015/28.-1495/6.*eta)*kappa_a+(26015/28.-44255/21.*eta-240*eta**2)*kappa_s)\
             + chi_a*chi_s*((26015/14.-88510/21.*eta-480*eta**2)*kappa_a+delta*(-1344475/1008. +
             745/18.*eta+(26015/14.-1495/3.*eta)*kappa_s))
    spin_phase_terms[7] = (-25150083775/3048192.+26804935/6048.*eta-1985/48.*eta**2)*delta*chi_a\
             + (-25150083775/3048192.+10566655595/762048.*eta-1042165/3024.*eta**2+5345/36.*eta**3)*chi_s
    
    #Cubic spin terms
    spin_phase_terms[7]+= +(14585/24.-2380*eta)*delta*chi_a**3\
     +(14585/24.-475/6.*eta+100/3.*eta**2)*chi_s**3+(14585/8.-215/2.*eta)*delta*chi_a*chi_s**2+\
     (14585/8.-7270*eta+80*eta**2)*chi_a**2*chi_s
     
    # 3.5PN spin-spin tail (Nagar et al. 2019, Eq. 47).
    # aQ2 includes the spin-spin cross term. kappa_i are spin-induced
    # quadrupole coefficients; for black holes, kappa_i=1.
    aA_tail = X1 * (chi_s + chi_a)
    aB_tail = X2 * (chi_s - chi_a)
    aQ2_tail = (kappa_1 * aA_tail**2 + 2.0*aA_tail*aB_tail
                + kappa_2*aB_tail**2)
    a12_tail = aA_tail - aB_tail
    spin_phase_terms[7] += -np.pi * (400.0 * aQ2_tail + 7.5*a12_tail**2)

    # 4PN spin-orbit tail (Nagar et al. 2019, Phys. Rev. D 99, 044007, Eq. 48).
    # a0 = X1*chi_1 + X2*chi_2; a12 = X1*chi_1 - X2*chi_2.
    # ln(pi*M*f) = 3*ln(v), with M in geometrized units.
    a0_4pn = chi_s + delta*chi_a
    a12_4pn = delta*chi_s + chi_a
    spin_phase_terms[8] = np.pi * (1.0 - 3.0*np.log(v)) * (
        (2388425.0/2268.0 - 9925.0*eta/27.0)*a0_4pn
        + (1538855.0/4536.0 - 19655.0*eta/756.0)*delta*a12_4pn
    )

    # Electric tides (5PN) 
    tidal_phase_terms[0] = -39/2.*lambda_tilde

    # Electric and magnetic quadrupolar tides (6PN).
    sigma_tilde = compute_sigma_tilde(X1, X2, sigma_1, sigma_2) if include_magnetic_6pn else 0.0
    tidal_phase_terms[1] = -3115/64.*lambda_tilde + 6595/364.*delta*delta_lambda_tilde + sigma_tilde

    # Spin-tidal couplings + tidal tail (6.5PN)
    chi_1 = chi_s + chi_a
    chi_2 = chi_s - chi_a

    lambda_hat = compute_lambda_hat(X1, X2, chi_1, chi_2, lambda_1, lambda_2)
    sigma_hat = compute_sigma_hat(X1, X2, chi_1, chi_2, sigma_1, sigma_2) if include_magnetic_spin_tidal else 0.0

    # Rotational tidal Love numbers default to zero.
    # All eight rotational Love numbers default to zero.
    gamma_hat = compute_gamma_hat(
        chi_1,
        chi_2,
        eta1,
        eta2,
        eta,
        m_tot,
        lambda1_23,
        sigma1_23,
        lambda1_32,
        sigma1_32,
        lambda2_23,
        sigma2_23,
        lambda2_32,
        sigma2_32,
    )

    kappa_hat = compute_kappa_hat(lambda_tilde)

    tidal_phase_terms[2] = lambda_hat + sigma_hat + gamma_hat + kappa_hat

    # Electric quadrupole and optional electric octupole (7PN).
    # The electric quadrupole and octupole coefficients are direct coefficients
    # of v^9 in the phase. The common prefactor is compensated below.
    if float(tidal_pn_order) >= 7.0:
        tidal_phase_terms[3] = compute_electric_quadrupole_7pn(
            mass_1, mass_2, lambda_1, lambda_2
        )
        if include_octupole_7pn:
            oct_term = -(X1**5 * L31 + X2**5 * L32)
            tidal_phase_terms[3] += (125.0 / 12.0) * oct_term

    # Selected 8PN sectors.
    if float(tidal_pn_order) >= 8.0:
        if include_quadratic_8pn:
            if lambda_8pn_quadratic is None:
                tidal_phase_terms[4] += compute_quadratic_love_8PN_term(
                    mass_1,
                    mass_2,
                    delta,
                    lambda222_1,
                    lambda222_2,
                )
            else:
                tidal_phase_terms[4] += (
                    -(45.0 / 11.0)
                    * float(lambda_8pn_quadratic)
                )
        if include_dynamic_8pn:
            tidal_phase_terms[4] += compute_dynamic_love_8PN_terms(
                v,
                mass_1,
                mass_2,
                eta,
                lambda_1,
                lambda_2,
                lambda_dyn_1,
                lambda_dyn_2,
                include_dynamic_log=include_dynamic_log,
                dynamic_log_gamma=dynamic_log_gamma,
            )
    
    #Build point-particle phase
    phi_PP=0
    if pn_phase_order==-1:
        dim = len(pp_phase_terms)
    else:
        dim = int(pn_phase_order)+1
    for i in range(dim):
        phi_PP -= 3/(128.*eta)*v**(-5)*(pp_phase_terms[i]+spin_phase_terms[i])*v**i
    phi_PP = -np.pi/4. + phi_PP
        
    # Build tidal phase
    phi_T = 3/(128.*eta)*v**(-5)*(tidal_phase_terms[0]*v**10+tidal_phase_terms[1]*v**12+tidal_phase_terms[2]*v**13+(128.*eta)/3*tidal_phase_terms[3]*v**14+tidal_phase_terms[4]*v**16)

    phi = phi_PP - phi_T
    
    return phi

    ########################### WAVEFORM ###############################

@np.errstate(divide='ignore', invalid='ignore')
def Taylor_f2_APR(frequency_array, chirp_mass, eta, chi_1, chi_2,\
                 lambda_tilde, delta_lambda_tilde,
                 phase, luminosity_distance, theta_jn, tilt_1=0., tilt_2=0., phi_12=0., phi_jl=0., **kwargs):

    # Reference frequency, minimum frequency, pN phase order
    waveform_kwargs = dict(reference_frequency=50.0, minimum_frequency = 5, pn_phase_order=-1)
    waveform_kwargs.update(kwargs)
    reference_frequency = waveform_kwargs['reference_frequency']
    pn_phase_order = waveform_kwargs['pn_phase_order']
    minimum_frequency = waveform_kwargs['minimum_frequency']
    maximum_frequency = waveform_kwargs.get('maximum_frequency', np.inf)
    magnetic_fluid = waveform_kwargs.get('magnetic_fluid', 'irrotational')
    include_magnetic_6pn = waveform_kwargs.get('include_magnetic_6pn', True)
    include_magnetic_spin_tidal = waveform_kwargs.get('include_magnetic_spin_tidal', True)
    include_dynamic_log = waveform_kwargs.get('include_dynamic_log', True)
    dynamic_log_gamma = waveform_kwargs.get('dynamic_log_gamma', 1.0)
    tidal_pn_order = float(waveform_kwargs.get('tidal_pn_order', 8.0))
    include_octupole_7pn = bool(
        waveform_kwargs.get('include_octupole_7pn', True)
    )
    include_quadratic_8pn = bool(
        waveform_kwargs.get('include_quadratic_8pn', True)
    )
    include_dynamic_8pn = bool(
        waveform_kwargs.get('include_dynamic_8pn', True)
    )

    lambda_8pn_quadratic = waveform_kwargs.get(
        "lambda_8pn_quadratic",
        None,
    )

    # Optional rotational tidal Love arguments for the 6.5PN gamma_hat term.
    lambda1_23 = waveform_kwargs.get('lambda1_23', 0.0)
    sigma1_23 = waveform_kwargs.get('sigma1_23', 0.0)
    lambda1_32 = waveform_kwargs.get('lambda1_32', 0.0)
    sigma1_32 = waveform_kwargs.get('sigma1_32', 0.0)
    lambda2_23 = waveform_kwargs.get('lambda2_23', 0.0)
    sigma2_23 = waveform_kwargs.get('sigma2_23', 0.0)
    lambda2_32 = waveform_kwargs.get('lambda2_32', 0.0)
    sigma2_32 = waveform_kwargs.get('sigma2_32', 0.0)

    # Mass, beta, quadrupole and tidal def.
    mass_1, mass_2 = mass_components(chirp_mass,eta)

    # Conversion to SI units                                      
    l_dist = luminosity_distance*m_parsec
    m_tot = mass_1 + mass_2                        
    total_mass = m_tot*solar_m

    # Symmetric and antisymmetric parameters
    chi_a = (chi_1-chi_2)/2.
    chi_s = (chi_1+chi_2)/2.
    
    # Electric quadrupolar tidal terms (5PN and 6PN)
    #lambda_1, lambda_2 = lambda_from_mass_APR(mass_1), lambda_from_mass_APR(mass_2)
    #lambda_s           = (lambda_1 + lambda_2)/2
    #lambda_a           = (lambda_1 - lambda_2)/2
    #lambda_tilde       = compute_lambda_tilde(eta,lambda_s,lambda_a)
    #delta_lambda_tilde = compute_delta_lambda_tilde(eta,lambda_s,lambda_a)

    # Magnetic tidals (6PN)
    #sigma_tilde        = sigma_tilde

    # Spin-tidal couplings (6.5PN)
    #lambda_hat         = compute_lambda_hat(eta1, eta2, chi_1, chi_2, lambda_1, lambda_2)
    #sigma_hat          = sigma_hat
    #gamma_hat          = compute_gamma_hat(chi_1, chi_2, eta1, eta2, eta, m_tot, lambda1_23, sigma1_23, lambda1_32, sigma1_32, lambda2_23, sigma2_23, lambda2_32, sigma2_32)
    
    # Tidal tail term (6.5PN)
    #kappa_hat          = compute_kappa_hat(lambda_tilde)
    
    # Rescaled frequency
    v = (g_newton/(light_s**3)*np.pi*frequency_array*total_mass)**(1/3.)

    # Phase
    phi = pn_phase(
        frequency_array,
        total_mass,
        eta,
        chi_s,
        chi_a,
        lambda_tilde,
        delta_lambda_tilde,
        lambda_8pn_quadratic=lambda_8pn_quadratic,
        lambda1_23=lambda1_23,
        sigma1_23=sigma1_23,
        lambda1_32=lambda1_32,
        sigma1_32=sigma1_32,
        lambda2_23=lambda2_23,
        sigma2_23=sigma2_23,
        lambda2_32=lambda2_32,
        sigma2_32=sigma2_32,
        pn_phase_order=pn_phase_order,
        magnetic_fluid=magnetic_fluid,
        include_magnetic_6pn=include_magnetic_6pn,
        include_magnetic_spin_tidal=include_magnetic_spin_tidal,
        include_dynamic_log=include_dynamic_log,
        dynamic_log_gamma=dynamic_log_gamma,
        tidal_pn_order=tidal_pn_order,
        include_octupole_7pn=include_octupole_7pn,
        include_quadratic_8pn=include_quadratic_8pn,
        include_dynamic_8pn=include_dynamic_8pn,
    )
    
    # Apply the reference-frequency and orbital-phase shifts.
    t_c=0
    phase_c = -pn_phase(
        np.array([reference_frequency]),
        total_mass,
        eta,
        chi_s,
        chi_a,
        lambda_tilde,
        delta_lambda_tilde,
        lambda_8pn_quadratic=lambda_8pn_quadratic,
        lambda1_23=lambda1_23,
        sigma1_23=sigma1_23,
        lambda1_32=lambda1_32,
        sigma1_32=sigma1_32,
        lambda2_23=lambda2_23,
        sigma2_23=sigma2_23,
        lambda2_32=lambda2_32,
        sigma2_32=sigma2_32,
        pn_phase_order=pn_phase_order,
        magnetic_fluid=magnetic_fluid,
        include_magnetic_6pn=include_magnetic_6pn,
        include_magnetic_spin_tidal=include_magnetic_spin_tidal,
        include_dynamic_log=include_dynamic_log,
        dynamic_log_gamma=dynamic_log_gamma,
        tidal_pn_order=tidal_pn_order,
        include_octupole_7pn=include_octupole_7pn,
        include_quadratic_8pn=include_quadratic_8pn,
        include_dynamic_8pn=include_dynamic_8pn,
    ).astype(float) - 3/4.*np.pi + 2*phase
    phi-= (2*np.pi*frequency_array*t_c-phase_c)
    
    # Amplitude
    amplitude = g_newton/(light_s**3)*total_mass*np.sqrt(5*np.pi/96.) *\
                (g_newton/(light_s**2)*total_mass/l_dist)*np.sqrt(eta)*(v)**(-7/2.)
    
    
    # Construct frequency-domain plus and cross strain
    ampl_plus = amplitude*(1+np.cos(theta_jn)**2)
    ampl_cross = amplitude*(2*np.cos(theta_jn))

    wave_plus = ampl_plus*(np.cos(phi)+1j*np.sin(phi))
    wave_cross = ampl_cross*(np.sin(phi)-1j*np.cos(phi))
    
    # Apply the caller-supplied hard analysis band.
    # The inference pipeline applies its smooth taper separately.
    outside_band = (
        (np.asarray(frequency_array) < minimum_frequency)
        | (np.asarray(frequency_array) > maximum_frequency)
    )
    wave_plus[outside_band] = 0.0 + 1j * 0.0
    wave_cross[outside_band] = 0.0 + 1j * 0.0

    # Returns:A dictionary with plus and cross polarization strain
    waveform = {'plus': wave_plus,
                'cross': wave_cross}

    return waveform#, amplitude, phi

def taylor_f2_APR_relative_binning(frequency_array, chirp_mass, eta, chi_1, chi_2,\
                 lambda_tilde, delta_lambda_tilde, 
                 phase, luminosity_distance, theta_jn, tilt_1=0., tilt_2=0., phi_12=0., phi_jl=0., **kwargs):
                 

    bin_edges = kwargs.get("frequency_bin_edges", None)
        
    # =========================
    # NO RELATIVE BINNING
    # =========================
    if bin_edges is None:
        return Taylor_f2_APR(frequency_array, chirp_mass, eta, chi_1, chi_2,\
                 lambda_tilde, delta_lambda_tilde, 
                 phase, luminosity_distance, theta_jn, tilt_1=0, tilt_2=0, phi_12=0, phi_jl=0, **kwargs)
                 
    # Modified for relative binning, returns the waveform on a set of frequencies corresponding to frequency-bin edges. 
    # =========================
    # RELATIVE BINNING MODE
    # =========================
    f_bin_edges = bin_edges
    
    return Taylor_f2_APR(f_bin_edges, chirp_mass, eta, chi_1, chi_2,\
                 lambda_tilde, delta_lambda_tilde,
                 phase, luminosity_distance, theta_jn, tilt_1=0., tilt_2=0., phi_12=0., phi_jl=0., **kwargs)


# =============================================================================
# Named source models for the model-comparison runs
# =============================================================================

def _relative_binning_with_model(
    frequency_array,
    chirp_mass,
    eta,
    chi_1,
    chi_2,
    lambda_tilde,
    delta_lambda_tilde,
    phase,
    luminosity_distance,
    theta_jn,
    *,
    tidal_pn_order,
    include_octupole_7pn,
    include_quadratic_8pn,
    include_dynamic_8pn,
    tilt_1=0.0,
    tilt_2=0.0,
    phi_12=0.0,
    phi_jl=0.0,
    **kwargs,
):
    model_kwargs = dict(kwargs)
    model_kwargs.update(
        tidal_pn_order=float(tidal_pn_order),
        include_octupole_7pn=bool(include_octupole_7pn),
        include_quadratic_8pn=bool(include_quadratic_8pn),
        include_dynamic_8pn=bool(include_dynamic_8pn),
    )
    return taylor_f2_APR_relative_binning(
        frequency_array,
        chirp_mass,
        eta,
        chi_1,
        chi_2,
        lambda_tilde,
        delta_lambda_tilde,
        phase,
        luminosity_distance,
        theta_jn,
        tilt_1=tilt_1,
        tilt_2=tilt_2,
        phi_12=phi_12,
        phi_jl=phi_jl,
        **model_kwargs,
    )


def taylor_f2_APR_8PN_complete_relative_binning(
    frequency_array,
    chirp_mass,
    eta,
    chi_1,
    chi_2,
    lambda_tilde,
    delta_lambda_tilde,
    phase,
    luminosity_distance,
    theta_jn,
    tilt_1=0.0,
    tilt_2=0.0,
    phi_12=0.0,
    phi_jl=0.0,
    **kwargs,
):
    return _relative_binning_with_model(
        frequency_array,
        chirp_mass,
        eta,
        chi_1,
        chi_2,
        lambda_tilde,
        delta_lambda_tilde,
        phase,
        luminosity_distance,
        theta_jn,
        tidal_pn_order=8.0,
        include_octupole_7pn=True,
        include_quadratic_8pn=True,
        include_dynamic_8pn=True,
        tilt_1=tilt_1,
        tilt_2=tilt_2,
        phi_12=phi_12,
        phi_jl=phi_jl,
        **kwargs,
    )


def taylor_f2_APR_8PN_quadratic_relative_binning(
    frequency_array,
    chirp_mass,
    eta,
    chi_1,
    chi_2,
    lambda_tilde,
    delta_lambda_tilde,
    phase,
    luminosity_distance,
    theta_jn,
    tilt_1=0.0,
    tilt_2=0.0,
    phi_12=0.0,
    phi_jl=0.0,
    **kwargs,
):
    return _relative_binning_with_model(
        frequency_array,
        chirp_mass,
        eta,
        chi_1,
        chi_2,
        lambda_tilde,
        delta_lambda_tilde,
        phase,
        luminosity_distance,
        theta_jn,
        tidal_pn_order=8.0,
        include_octupole_7pn=True,
        include_quadratic_8pn=True,
        include_dynamic_8pn=False,
        tilt_1=tilt_1,
        tilt_2=tilt_2,
        phi_12=phi_12,
        phi_jl=phi_jl,
        **kwargs,
    )


def taylor_f2_APR_8PN_dynamic_relative_binning(
    frequency_array,
    chirp_mass,
    eta,
    chi_1,
    chi_2,
    lambda_tilde,
    delta_lambda_tilde,
    phase,
    luminosity_distance,
    theta_jn,
    tilt_1=0.0,
    tilt_2=0.0,
    phi_12=0.0,
    phi_jl=0.0,
    **kwargs,
):
    return _relative_binning_with_model(
        frequency_array,
        chirp_mass,
        eta,
        chi_1,
        chi_2,
        lambda_tilde,
        delta_lambda_tilde,
        phase,
        luminosity_distance,
        theta_jn,
        tidal_pn_order=8.0,
        include_octupole_7pn=True,
        include_quadratic_8pn=False,
        include_dynamic_8pn=True,
        tilt_1=tilt_1,
        tilt_2=tilt_2,
        phi_12=phi_12,
        phi_jl=phi_jl,
        **kwargs,
    )


def taylor_f2_APR_7PN_relative_binning(
    frequency_array,
    chirp_mass,
    eta,
    chi_1,
    chi_2,
    lambda_tilde,
    delta_lambda_tilde,
    phase,
    luminosity_distance,
    theta_jn,
    tilt_1=0.0,
    tilt_2=0.0,
    phi_12=0.0,
    phi_jl=0.0,
    **kwargs,
):
    return _relative_binning_with_model(
        frequency_array,
        chirp_mass,
        eta,
        chi_1,
        chi_2,
        lambda_tilde,
        delta_lambda_tilde,
        phase,
        luminosity_distance,
        theta_jn,
        tidal_pn_order=7.0,
        include_octupole_7pn=True,
        include_quadratic_8pn=False,
        include_dynamic_8pn=False,
        tilt_1=tilt_1,
        tilt_2=tilt_2,
        phi_12=phi_12,
        phi_jl=phi_jl,
        **kwargs,
    )


def taylor_f2_APR_7PN_nooct_relative_binning(
    frequency_array,
    chirp_mass,
    eta,
    chi_1,
    chi_2,
    lambda_tilde,
    delta_lambda_tilde,
    phase,
    luminosity_distance,
    theta_jn,
    tilt_1=0.0,
    tilt_2=0.0,
    phi_12=0.0,
    phi_jl=0.0,
    **kwargs,
):
    return _relative_binning_with_model(
        frequency_array,
        chirp_mass,
        eta,
        chi_1,
        chi_2,
        lambda_tilde,
        delta_lambda_tilde,
        phase,
        luminosity_distance,
        theta_jn,
        tidal_pn_order=7.0,
        include_octupole_7pn=False,
        include_quadratic_8pn=False,
        include_dynamic_8pn=False,
        tilt_1=tilt_1,
        tilt_2=tilt_2,
        phi_12=phi_12,
        phi_jl=phi_jl,
        **kwargs,
    )


def taylor_f2_APR_6p5PN_relative_binning(
    frequency_array,
    chirp_mass,
    eta,
    chi_1,
    chi_2,
    lambda_tilde,
    delta_lambda_tilde,
    phase,
    luminosity_distance,
    theta_jn,
    tilt_1=0.0,
    tilt_2=0.0,
    phi_12=0.0,
    phi_jl=0.0,
    **kwargs,
):
    return _relative_binning_with_model(
        frequency_array,
        chirp_mass,
        eta,
        chi_1,
        chi_2,
        lambda_tilde,
        delta_lambda_tilde,
        phase,
        luminosity_distance,
        theta_jn,
        tidal_pn_order=6.5,
        include_octupole_7pn=False,
        include_quadratic_8pn=False,
        include_dynamic_8pn=False,
        tilt_1=tilt_1,
        tilt_2=tilt_2,
        phi_12=phi_12,
        phi_jl=phi_jl,
        **kwargs,
    )


def taylor_f2_APR_8PN_complete_agnostic_quadratic_relative_binning(
    frequency_array,
    chirp_mass,
    eta,
    chi_1,
    chi_2,
    lambda_tilde,
    delta_lambda_tilde,
    lambda_8pn_quadratic,
    phase,
    luminosity_distance,
    theta_jn,
    tilt_1=0.0,
    tilt_2=0.0,
    phi_12=0.0,
    phi_jl=0.0,
    **kwargs,
):
    """
    Complete 8PN model in which the effective static quadratic
    coefficient

        lambda_8pn_quadratic
            = 19 rho_tilde_plus
              - 17 delta rho_tilde_minus

    is sampled independently.

    All other sectors, including the dynamical 8PN response,
    retain their quasi-universal relations.
    """

    return _relative_binning_with_model(
        frequency_array,
        chirp_mass,
        eta,
        chi_1,
        chi_2,
        lambda_tilde,
        delta_lambda_tilde,
        phase,
        luminosity_distance,
        theta_jn,
        tidal_pn_order=8.0,
        include_octupole_7pn=True,
        include_quadratic_8pn=True,
        include_dynamic_8pn=True,
        lambda_8pn_quadratic=lambda_8pn_quadratic,
        tilt_1=tilt_1,
        tilt_2=tilt_2,
        phi_12=phi_12,
        phi_jl=phi_jl,
        **kwargs,
    )
