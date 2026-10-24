"""Electric-octupole power-law fit to digitized data.
See docs/QUR.md and examples/fit_octupole_qur.py for its inputs.
"""
import numpy as np
A = 0.23214267180912393
B = 1.345740158802329

def lambda3_from_lambda_fit(L):
    values = np.asarray(L, dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError("Lambda_2 must be finite and nonnegative")
    result = A * values**B
    return result.item() if result.shape == () else result

lambda3_over_lambda_fit = lambda3_from_lambda_fit  # Compatibility alias: returns Lambda_3.
