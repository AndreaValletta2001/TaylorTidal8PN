"""Thin wrapper selecting the 8PN quadratic Love only tidal model."""

from .TaylorF2_APR_tidal_orders_fRoche import *
from .TaylorF2_APR_tidal_orders_fRoche import (
    taylor_f2_APR_8PN_quadratic_relative_binning as taylor_f2_APR_relative_binning,
)

MODEL_NAME = '8PN quadratic Love only'
