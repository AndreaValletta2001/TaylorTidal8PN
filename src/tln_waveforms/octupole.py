"""Electric octupolar tidal deformability from the multipolar
tidal Love quasi-universal relation of Yagi (2014).

The relation is a fourth-order polynomial in ln(Lambda_2).
"""

import numpy as np


def lambda3_from_lambda_fit(L):
    """Compute Lambda_3 from Lambda_2 using Yagi's QUR.

    Lambda_2 = 0 is mapped to Lambda_3 = 0 (black-hole limit).
    Negative or non-finite values are rejected.
    """
    values = np.asarray(L, dtype=float)

    if np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError("Lambda_2 must be finite and nonnegative")

    result = np.zeros_like(values, dtype=float)
    mask = values > 0

    x = np.log(values[mask])

    result[mask] = np.exp(
        -1.15
        + 1.18 * x
        + 0.0251 * x**2
        - 0.00131 * x**3
        + 0.0000252 * x**4
    )

    return result.item() if result.shape == () else result


# Compatibility alias: returns Lambda_3.
lambda3_over_lambda_fit = lambda3_from_lambda_fit
