"""Core gamma-model utilities extracted from the master notebook."""

import numpy as np
from scipy.special import gamma


def model_gamma(t, R_peak, t_peak, t_onset, t_mean, epsilon=0.1):
    """
    Compute a shifted gamma-like emission curve.

    Args:
        t (array): Time values.
        R_peak (float): Peak emission amplitude.
        t_peak (float): Time of peak emission.
        t_onset (float): Onset time of emission.
        t_mean (float): Mean-related shape control parameter.
        epsilon (float): Smoothing factor for onset ramp.

    Returns:
        array: Predicted emission values.
    """
    alpha = (t_mean - t_onset) / (t_mean - t_peak)
    beta = t_mean - t_peak
    t_adj = t - t_onset

    core = ((np.clip(t_adj, 1e-10, None) / (t_peak - t_onset)) ** (alpha - 1)) * np.exp(
        -(t_adj - (t_peak - t_onset)) / beta
    )

    ramp = 1 / (1 + np.exp(-(t - t_onset) / epsilon))

    model_pred = R_peak * ramp * core

    model_pred = np.nan_to_num(model_pred, nan=0.0, posinf=0.0, neginf=0.0)
    return model_pred


def shape(t_onset, t_peak, t_mean):
    """
    Compute the shape ratio of the emission curve.

    Args:
        t_onset (float): Onset time.
        t_peak (float): Peak time.
        t_mean (float): Mean time.

    Returns:
        float: Shape ratio.
    """
    return (t_mean - t_peak) / (t_mean - t_onset)


def duration(t_onset, t_peak, t_mean):
    """
    Compute the duration between onset and mean time.

    Args:
        t_onset (float): Onset time.
        t_peak (float): Peak time.
        t_mean (float): Mean time.

    Returns:
        float: Duration value.
    """
    return t_mean - t_onset


def total_integral_model_gamma(R_peak, t_peak, t_onset, t_mean):
    """
    Compute the total integral of the gamma-based emission model.

    Args:
        R_peak (float): Peak emission amplitude.
        t_peak (float): Peak time.
        t_onset (float): Onset time.
        t_mean (float): Mean time.

    Returns:
        float: Total integrated emission.
    """
    alpha = (t_mean - t_onset) / (t_mean - t_peak)
    beta = t_mean - t_peak
    t0 = t_peak - t_onset
    A = R_peak * (t0 ** (1 - alpha)) * np.exp(t0 / beta)
    total = A * (beta**alpha) * gamma(alpha)

    return total
