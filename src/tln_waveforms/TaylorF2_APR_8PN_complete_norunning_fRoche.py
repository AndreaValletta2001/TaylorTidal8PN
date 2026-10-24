"""8PN complete QUR model with beta_2 running set identically to zero."""

from .TaylorF2_APR_tidal_orders_agnostic_total_norunning_fRoche import *

from .TaylorF2_APR_tidal_orders_agnostic_total_norunning_fRoche import (
    taylor_f2_APR_8PN_complete_relative_binning
    as taylor_f2_APR_relative_binning,
)

MODEL_NAME = (
    "8PN complete QUR, quadratic + dynamical, beta2=0"
)
