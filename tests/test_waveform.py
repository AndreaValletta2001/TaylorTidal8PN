import importlib
import math
import numpy as np
import unittest
from tln_waveforms import sigma_from_lambda_fit
from tln_waveforms import TaylorF2_APR_tidal_orders_fRoche as wf

ENGINES=['TaylorF2_APR_tidal_orders_fRoche','TaylorF2_APR_tidal_orders_agnostic_quadratic_fRoche','TaylorF2_APR_tidal_orders_agnostic_total_norunning_fRoche']
PARAMS=dict(chirp_mass=1.2,eta=.24,chi_1=.03,chi_2=-.01,lambda_tilde=400.,delta_lambda_tilde=0.,phase=0.,luminosity_distance=100.,theta_jn=.5,reference_frequency=100.,minimum_frequency=20.,maximum_frequency=1000.)

def test_published_magnetic_fit(L):
    x=math.log(L)
    expected=-math.exp(-2.03+.487*x+.00969*x*x+.00103*x**3-9.37e-5*x**4+2.24e-6*x**5)
    np.testing.assert_allclose(sigma_from_lambda_fit(L),expected,rtol=2e-14)
    assert sigma_from_lambda_fit(L,'static')>0

def test_invalid_love_numbers(bad):
    with unittest.TestCase().assertRaises(ValueError):sigma_from_lambda_fit(bad)

def test_magnetic_phase_terms_and_waveform_options(engine):
    mod=importlib.import_module('tln_waveforms.'+engine)
    f=np.array([20.,100.,300.,700.]);eta=.24;mt=2.8
    args=dict(frequency_array=f,total_mass=mt*mod.solar_m,eta=eta,chi_s=.01,chi_a=.02,lambda_tilde=400.,delta_lambda_tilde=0.)
    m1,m2=mod.mass_components(mt*eta**.6,eta);x1,x2=m1/mt,m2/mt
    l1,l2=mod.lambda_conv(400.,0.,m1,m2)
    s1,s2=sigma_from_lambda_fit(l1),sigma_from_lambda_fit(l2)
    v=(mod.g_newton/mod.light_s**3*np.pi*f*mt*mod.solar_m)**(1/3)
    full=mod.pn_phase(**args)
    no6=mod.pn_phase(**args,include_magnetic_6pn=False)
    coeff=(6920/7-20740/(21*x1))*x1**5*s1+(6920/7-20740/(21*x2))*x2**5*s2
    np.testing.assert_allclose(full-no6,-3/(128*eta)*coeff*v**7,rtol=1e-7,atol=2e-10)
    c1,c2=.03,-.01
    hat=((-9865/3+4933/(3*x1)+1644*x1)*c2-c1)*x1**5*s1+((-9865/3+4933/(3*x2)+1644*x2)*c1-c2)*x2**5*s2
    nohat=mod.pn_phase(**args,include_magnetic_spin_tidal=False)
    np.testing.assert_allclose(full-nohat,-3/(128*eta)*hat*v**8,rtol=1e-6,atol=2e-10)
    ir=mod.Taylor_f2_APR(f,**PARAMS)
    st=mod.Taylor_f2_APR(f,**PARAMS,magnetic_fluid='static')
    assert not np.allclose(ir['plus']/st['plus'],1,rtol=1e-9,atol=1e-9)
    for pol in ir:
        np.testing.assert_allclose(abs(ir[pol]),abs(st[pol]),rtol=1e-13,atol=1e-40)
        assert np.isfinite(ir[pol]).all()

def test_all_public_models_and_relative_binning(name):
    mod=importlib.import_module('tln_waveforms.TaylorF2_APR_'+name+'_fRoche')
    f=np.array([10.,20.,100.,300.,700.,1100.])
    params=dict(PARAMS)
    if 'agnostic_quadratic' in name: params['lambda_8pn_quadratic']=1000.
    if 'agnostic_total' in name: params['lambda_8pn_total']=1000.
    h=mod.taylor_f2_APR_relative_binning(f,**params)
    b=mod.taylor_f2_APR_relative_binning(f,**params,frequency_bin_edges=f[1:5])
    for pol in h:
        assert np.isfinite(h[pol]).all()
        assert h[pol][0]==h[pol][-1]==0
        np.testing.assert_array_equal(b[pol],h[pol][1:5])

def test_no_tide_and_equal_mass_limits():
    assert sigma_from_lambda_fit(0)==0
    h=wf.Taylor_f2_APR(np.array([20.,100.,300.]),**(PARAMS|dict(eta=.25,lambda_tilde=0.,delta_lambda_tilde=0.)))
    assert all(np.isfinite(x).all() for x in h.values())

class WaveformChecks(unittest.TestCase):
    def test_polynomial(self):
        for L in [1.,10.,100.,300.,1000.,1e4,1e5]:
            with self.subTest(L=L): test_published_magnetic_fit(L)
    def test_invalid(self):
        for bad in [-1.,float('nan'),float('inf')]:
            with self.subTest(bad=bad): test_invalid_love_numbers(bad)
    def test_magnetic_engines(self):
        for engine in ENGINES:
            with self.subTest(engine=engine): test_magnetic_phase_terms_and_waveform_options(engine)
    def test_public_models(self):
        for name in ['6p5PN','7PN','7PN_nooct','8PN_complete','8PN_quadratic','8PN_dynamic','8PN_complete_norunning','8PN_complete_agnostic_quadratic','8PN_complete_agnostic_total_norunning']:
            with self.subTest(model=name): test_all_public_models_and_relative_binning(name)
    def test_limits(self):
        test_no_tide_and_equal_mass_limits()

if __name__ == '__main__': unittest.main()
