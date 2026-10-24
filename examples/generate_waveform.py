"""Generate an 8PN waveform and save portable arrays; run from any directory."""
import numpy as np
from tln_waveforms import waveform_8pn
f=np.linspace(20.,1000.,4096)
h=waveform_8pn(f,chirp_mass=1.2,eta=.24,chi_1=0.,chi_2=0.,lambda_tilde=400.,delta_lambda_tilde=0.,phase=0.,luminosity_distance=100.,theta_jn=.5,reference_frequency=100.,minimum_frequency=20.,maximum_frequency=1000.,magnetic_fluid='irrotational')
np.savez('waveform_example.npz',frequency=f,plus=h['plus'],cross=h['cross'])
print('Saved waveform_example.npz')
