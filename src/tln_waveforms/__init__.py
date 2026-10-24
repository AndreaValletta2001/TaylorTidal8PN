"""TaylorF2 tidal models: selected tidal sectors through 8PN."""
__version__ = "0.1.0rc2"
from .TaylorF2_APR_tidal_orders_fRoche import (
    Taylor_f2_APR, pn_phase,
    taylor_f2_APR_8PN_complete_relative_binning as waveform_8pn,
    taylor_f2_APR_8PN_quadratic_relative_binning as waveform_8pn_quadratic,
    taylor_f2_APR_8PN_dynamic_relative_binning as waveform_8pn_dynamic,
    taylor_f2_APR_7PN_relative_binning as waveform_7pn,
    taylor_f2_APR_7PN_nooct_relative_binning as waveform_7pn_nooct,
    taylor_f2_APR_6p5PN_relative_binning as waveform_6p5pn,
)
from .magnetic import sigma_from_lambda_fit
