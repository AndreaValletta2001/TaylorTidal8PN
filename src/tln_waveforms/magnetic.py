"""Magnetic quadrupolar QUR: arXiv:1807.08016v3, Eq. (22), Table I.

Sigma = sigma_2 / m**5 in G=c=1; natural logarithms.
Irrotational is the inspiral default. Static is an explicit comparison branch.
Zero maps to zero as the no-tide convention (the fit itself requires Lambda>0).
No refitting or file/plot side effects occur at import.
"""
import numpy as np

COEFFICIENTS = {
    "irrotational": (-2.03, .487, .00969, .00103, -9.37e-5, 2.24e-6),
    "static": (-2.66, .786, -.01, .00128, -6.37e-5, 1.18e-6),
}

def sigma_from_lambda_fit(L, fluid="irrotational"):
    """Return dimensionless Sigma_2 for scalar or array nonnegative Lambda_2."""
    if fluid not in COEFFICIENTS:
        raise ValueError("fluid must be 'irrotational' or 'static'")
    values = np.asarray(L, dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError("Lambda_2 must be finite and nonnegative")
    positive = values > 0
    # Use a safe logarithm at zero and enforce the no-tide convention below.
    x = np.log(np.where(positive, values, 1.0))
    y = np.polynomial.polynomial.polyval(x, COEFFICIENTS[fluid])
    sign = -1.0 if fluid == "irrotational" else 1.0
    with np.errstate(over="raise", invalid="raise"):
        result = np.where(positive, sign * np.exp(y), 0.0)
    return result.item() if result.shape == () else result

def sigma_over_lambda_fit(L, fluid="irrotational"):
    """Return Sigma_2/Lambda_2; define the no-tide value at zero as zero."""
    values = np.asarray(L, dtype=float)
    sigma = sigma_from_lambda_fit(values, fluid=fluid)
    out = np.zeros_like(values)
    np.divide(sigma, values, out=out, where=values > 0)
    return out.item() if out.shape == () else out
