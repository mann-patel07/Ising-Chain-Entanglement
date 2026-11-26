import numpy as np

from central_charge import chord_periodic, chord_open, fit_c, c_from_two_sizes


def test_fit_periodic():
    L = 40
    ls = np.arange(2, 39)
    S = (0.5 / 3) * chord_periodic(L, ls) + 0.3
    c, b = fit_c(ls, S, L, "periodic")
    assert abs(c - 0.5) < 1e-12 and abs(b - 0.3) < 1e-12


def test_open_fit_and_two_sizes():
    L = 64
    ls = np.arange(4, 61)
    S = (0.5 / 6) * chord_open(L, ls) + 0.3
    c, _ = fit_c(ls, S, L, "open")
    assert abs(c - 0.5) < 1e-12
    S1 = (0.5 / 6) * chord_open(100, 50) + 0.3
    S2 = (0.5 / 6) * chord_open(200, 100) + 0.3
    assert abs(c_from_two_sizes(S1, 100, S2, 200, "open") - 0.5) < 1e-12
