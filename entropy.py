"""Entanglement entropy of pure states stored as vectors of length 2^L (site 0 is the most significant bit)."""
import numpy as np


def schmidt_values(psi, L, l):
    """All squared Schmidt values for the cut after l sites, sorted from largest to smallest."""
    s = np.linalg.svd(psi.reshape(2**l, 2**(L - l)), compute_uv=False)
    return np.sort(s**2)[::-1]


def entanglement_entropy(psi, L, l):
    """Von Neumann entropy (natural log) of sites 0..l-1 of the pure state psi (length 2^L)."""
    if l == 0 or l == L:
        return 0.0
    s = np.linalg.svd(psi.reshape(2**l, 2**(L - l)), compute_uv=False)
    p = s**2
    p = p[p > 1e-16]
    p = p / p.sum()
    return float(-np.sum(p * np.log(p)))


def entropy_profile(psi, L):
    """Array of S(l) for l = 0, 1, ..., L."""
    return np.array([entanglement_entropy(psi, L, l) for l in range(L + 1)])
