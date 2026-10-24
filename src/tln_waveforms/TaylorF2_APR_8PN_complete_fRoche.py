"""Thin wrapper selecting the 8PN complete (quadratic + dynamical) tidal model."""

from .TaylorF2_APR_tidal_orders_fRoche import *
from .TaylorF2_APR_tidal_orders_fRoche import (
    taylor_f2_APR_8PN_complete_relative_binning as taylor_f2_APR_relative_binning,
)

MODEL_NAME = '8PN complete (quadratic + dynamical)'
