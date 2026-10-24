"""Independent sign, spin-tail and public-input regression tests."""
import importlib
import unittest
import numpy as np

ENGINES = ['TaylorF2_APR_tidal_orders_fRoche',
           'TaylorF2_APR_tidal_orders_agnostic_quadratic_fRoche',
           'TaylorF2_APR_tidal_orders_agnostic_total_norunning_fRoche']

class PhaseCorrections(unittest.TestCase):
    def test_electric_tidal_sign(self):
        # At zero spin, disabling magnetism and truncating below 7PN leaves
        # the published electric 5PN, 6PN and 6.5PN terms only.
        for name in ENGINES:
            with self.subTest(engine=name):
                w=importlib.import_module('tln_waveforms.'+name)
                f=np.array([20.,100.,400.,900.]);eta=.24;lt=500.;dlt=10.
                args=dict(frequency_array=f,total_mass=2.8*w.solar_m,eta=eta,
                          chi_s=0.,chi_a=0.,tidal_pn_order=6.5,
                          include_magnetic_6pn=False,include_magnetic_spin_tidal=False)
                a=w.pn_phase(**args,lambda_tilde=lt,delta_lambda_tilde=dlt)
                b=w.pn_phase(**args,lambda_tilde=0.,delta_lambda_tilde=0.)
                v=(np.pi*w.g_newton*2.8*w.solar_m*f/w.light_s**3)**(1/3)
                psi=3/(128*eta)*(-39/2*lt*v**5+(-3115/64*lt+6595/364*np.sqrt(1-4*eta)*dlt)*v**7+39*np.pi/2*lt*v**8)
                np.testing.assert_allclose(a-b,-psi,rtol=1e-8,atol=2e-10)

    def test_four_pn_spin_orbit_tail(self):
        # Independent symmetric/antisymmetric spin basis, Mishra Eq. (6c).
        for name in ENGINES:
            with self.subTest(engine=name):
                w=importlib.import_module('tln_waveforms.'+name)
                f=np.array([20.,100.,400.,900.]);eta=.24;s=.05;a=.08
                args=dict(frequency_array=f,total_mass=2.8*w.solar_m,eta=eta,
                          chi_s=s,chi_a=a,lambda_tilde=500.,delta_lambda_tilde=10.)
                v=(np.pi*w.g_newton*2.8*w.solar_m*f/w.light_s**3)**(1/3)
                c=np.pi*((233915/168-99185/252*eta)*np.sqrt(1-4*eta)*a+(233915/168-3970375/2268*eta+19655/189*eta**2)*s)*(1-3*np.log(v))
                np.testing.assert_allclose(w.pn_phase(**args,pn_phase_order=8)-w.pn_phase(**args,pn_phase_order=7),-3/(128*eta)*c*v**3,rtol=2e-8,atol=2e-10)
                np.testing.assert_allclose(w.pn_phase(**args),w.pn_phase(**args,pn_phase_order=8),rtol=0,atol=0)

    def test_three_half_pn_spin_spin_tail(self):
        # The even-in-spin part of the 3.5PN increment isolates this tail:
        # linear and cubic terms cancel, and the zero-spin term is removed.
        for name in ENGINES:
            with self.subTest(engine=name):
                w=importlib.import_module('tln_waveforms.'+name)
                f=np.array([30.,200.,800.]);eta=.24;lt=500.;dlt=10.;s=.05;a=.08
                args=dict(frequency_array=f,total_mass=2.8*w.solar_m,eta=eta,
                          lambda_tilde=lt,delta_lambda_tilde=dlt)
                def increment(cs,ca):
                    return w.pn_phase(**args,chi_s=cs,chi_a=ca,pn_phase_order=7)-w.pn_phase(**args,chi_s=cs,chi_a=ca,pn_phase_order=6)
                actual=(increment(s,a)+increment(-s,-a))/2-increment(0.,0.)
                X=(1+np.sqrt(1-4*eta))/2;Y=1-X
                l1,l2=w.lambda_conv(lt,dlt,2.8*X,2.8*Y)
                def kappa(L):
                    x=np.log(L)
                    return np.exp(.194+.0936*x+.0474*x*x-.00421*x**3+.000123*x**4)
                A=X*(s+a);B=Y*(s-a)
                c=-np.pi*(400*(kappa(l1)*A*A+2*A*B+kappa(l2)*B*B)+7.5*(A-B)**2)
                v=(np.pi*w.g_newton*2.8*w.solar_m*f/w.light_s**3)**(1/3)
                np.testing.assert_allclose(actual,-3/(128*eta)*c*v*v,rtol=1e-8,atol=5e-10)

    def test_seven_pn_electric_quadrupole(self):
        # Henry, Faye & Blanchet (2020) with Erratum, arXiv:2005.13367v4:
        # general SPA phase of Sec. VI and, separately, its equal-mass form.
        for name in ENGINES:
            with self.subTest(engine=name):
                w=importlib.import_module('tln_waveforms.'+name)
                f=np.array([30.,200.,800.])
                for eta,lt,dlt in [(.24,500.,10.),(.2,300.,-40.)]:
                    args=dict(frequency_array=f,total_mass=2.8*w.solar_m,eta=eta,chi_s=0.,chi_a=0.,
                              lambda_tilde=lt,delta_lambda_tilde=dlt,include_octupole_7pn=False)
                    actual=w.pn_phase(**args,tidal_pn_order=7)-w.pn_phase(**args,tidal_pn_order=6.5)
                    X=(1+np.sqrt(1-4*eta))/2;Y=1-X
                    l1,l2=w.lambda_conv(lt,dlt,2.8*X,2.8*Y)
                    mp=(Y*X**4*l1+X*Y**4*l2)/2;mm=(Y*X**4*l1-X*Y**4*l2)/2
                    c=-9/(16*eta**2)*((136190135/27433728+978554825/4572288*eta-281935/6048*eta**2+5/3*eta**3)*mp
                                      +(X-Y)*(136190135/27433728+213905/2592*eta+1585/1296*eta**2)*mm)
                    v=(np.pi*w.g_newton*2.8*w.solar_m*f/w.light_s**3)**(1/3)
                    np.testing.assert_allclose(actual,-c*v**9,rtol=1e-8,atol=1e-10)
                # Equal bodies: psi = -(117/2) mu v^5 [... + 29323235/3429216 mu v^4], mu = Lambda/32.
                np.testing.assert_allclose(w.coeffLT(1.4,1.4),-117/64*29323235/3429216,rtol=1e-12)

    def test_spin_induced_quadrupole_small_lambda(self):
        # Below Lambda=1 the fit is frozen at its Lambda=1 value; zero is a black hole.
        for name in ENGINES:
            with self.subTest(engine=name):
                w=importlib.import_module('tln_waveforms.'+name)
                for L in [1e-8,1e-4,.5,1.]:
                    self.assertAlmostEqual(w.quadrupole_from_loveq(L),np.exp(.194),places=12)
                self.assertEqual(w.quadrupole_from_loveq(0.),1)

    def test_numpy_error_state_is_untouched(self):
        for name in ENGINES:
            with self.subTest(engine=name):
                w=importlib.import_module('tln_waveforms.'+name)
                with np.errstate(divide='warn',invalid='warn'):
                    before=np.geterr()
                    w.Taylor_f2_APR(np.array([0.,20.,100.]),chirp_mass=1.2,eta=.24,chi_1=0.,chi_2=0.,
                                    lambda_tilde=400.,delta_lambda_tilde=0.,phase=0.,luminosity_distance=100.,
                                    theta_jn=.5,minimum_frequency=20.)
                    self.assertEqual(np.geterr(),before)


    def test_rotational_tidal_love_numbers(self):
        for name in ENGINES:
            with self.subTest(engine=name):
                w = importlib.import_module("tln_waveforms." + name)

                X1, X2, eta = 0.6, 0.4, 0.24
                chi1, chi2 = 0.1, -0.05
                args = (chi1, chi2, X1, X2, eta, 2.8)

                # The default contribution must vanish identically.
                self.assertEqual(
                    w.compute_gamma_hat(*args, *([0.0]*8)),
                    0.0
                )

                # Check nonzero inputs against Eq. (10).
                values = (10., -3., 2., 1., 5., -2., 4., 3.)

                expected = (
                    chi1 * X1**4 * (
                        (856*X1-816*X1**2)*values[0]
                        -(833*X1/3-278*X1**2)*values[1]
                        -eta*(272*values[2]-204*values[3])
                    )
                    + chi2 * X2**4 * (
                        (856*X2-816*X2**2)*values[4]
                        -(833*X2/3-278*X2**2)*values[5]
                        -eta*(272*values[6]-204*values[7])
                    )
                )

                np.testing.assert_allclose(
                    w.compute_gamma_hat(*args, *values),
                    expected,
                    rtol=1e-13
                )

                f = np.array([40., 100., 300.])
                common = dict(
                    frequency_array=f,
                    total_mass=2.8*w.solar_m,
                    eta=eta,
                    chi_s=(chi1+chi2)/2,
                    chi_a=(chi1-chi2)/2,
                    lambda_tilde=400.,
                    delta_lambda_tilde=0.,
                )

                phi_zero = w.pn_phase(**common)
                phi_rtln = w.pn_phase(
                    **common,
                    lambda1_23=values[0],
                    sigma1_23=values[1],
                    lambda1_32=values[2],
                    sigma1_32=values[3],
                    lambda2_23=values[4],
                    sigma2_23=values[5],
                    lambda2_32=values[6],
                    sigma2_32=values[7],
                )

                v = (
                    np.pi*w.g_newton*2.8*w.solar_m*f/w.light_s**3
                )**(1/3)

                np.testing.assert_allclose(
                    phi_rtln-phi_zero,
                    -3/(128*eta)*expected*v**8,
                    rtol=1e-6,
                    atol=1e-9,
                )


if __name__ == '__main__':
    unittest.main()
