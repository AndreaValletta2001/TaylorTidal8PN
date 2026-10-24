"""
8PN recovery with one independent parameter representing
the total quadratic+dynamical 8PN coefficient, with beta_2=0.
"""

from .TaylorF2_APR_tidal_orders_agnostic_total_norunning_fRoche import *

from .TaylorF2_APR_tidal_orders_agnostic_total_norunning_fRoche import (
    taylor_f2_APR_8PN_complete_agnostic_total_relative_binning
    as taylor_f2_APR_relative_binning,
)

MODEL_NAME = (
    "8PN agnostic total quadratic+dynamical, beta2=0"
)
