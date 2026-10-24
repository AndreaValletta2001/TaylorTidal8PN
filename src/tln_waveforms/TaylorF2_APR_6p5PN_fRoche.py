"""Thin wrapper selecting the 6.5PN tidal model."""

from .TaylorF2_APR_tidal_orders_fRoche import *
from .TaylorF2_APR_tidal_orders_fRoche import (
    taylor_f2_APR_6p5PN_relative_binning as taylor_f2_APR_relative_binning,
)

MODEL_NAME = '6.5PN'
