"""Calabrese and Cardy central charge estimators.

periodic: S(l) = (c/3) ln[(L/pi) sin(pi l / L)] + const
open, block at one end: S(l) = (c/6) ln[(2L/pi) sin(pi l / L)] + const
"""
import numpy as np


def chord_periodic(L, l):
    return np.log((L / np.pi) * np.sin(np.pi * np.asarray(l) / L))


def chord_open(L, l):
    return np.log((2 * L / np.pi) * np.sin(np.pi * np.asarray(l) / L))


def _factor(boundary):
    if boundary == "periodic":
        return 3.0
    if boundary == "open":
        return 6.0
    raise ValueError("boundary must be 'open' or 'periodic'")


def fit_c(ls, S, L, boundary):
    """Least squares fit of S against the chord log. Returns (c, intercept)."""
    f = _factor(boundary)
    chord = chord_periodic if boundary == "periodic" else chord_open
    x = chord(L, np.asarray(ls, dtype=float))
    slope, intercept = np.polyfit(x, np.asarray(S, dtype=float), 1)
    return float(f * slope), float(intercept)


def c_from_two_sizes(S_half_1, L1, S_half_2, L2, boundary):
    """c from the half chain entropy at two sizes."""
    return float(_factor(boundary) * (S_half_2 - S_half_1) / np.log(L2 / L1))
