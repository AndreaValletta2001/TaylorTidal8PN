"""
Complete 8PN waveform with an independently sampled effective
static quadratic-Love coefficient.
"""

from .TaylorF2_APR_tidal_orders_agnostic_quadratic_fRoche import *
from .TaylorF2_APR_tidal_orders_agnostic_quadratic_fRoche import (
    taylor_f2_APR_8PN_complete_agnostic_quadratic_relative_binning
    as taylor_f2_APR_relative_binning,
)

MODEL_NAME = (
    "8PN complete with agnostic static quadratic Love"
)
