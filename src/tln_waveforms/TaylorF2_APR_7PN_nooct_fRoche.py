"""Thin wrapper selecting the 7PN without electric octupole tidal model."""

from .TaylorF2_APR_tidal_orders_fRoche import *
from .TaylorF2_APR_tidal_orders_fRoche import (
    taylor_f2_APR_7PN_nooct_relative_binning as taylor_f2_APR_relative_binning,
)

MODEL_NAME = '7PN without electric octupole'
