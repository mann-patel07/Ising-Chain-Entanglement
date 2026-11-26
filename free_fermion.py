"""Exact solution of the open transverse field Ising chain through Majorana fermions.

Use only for the open boundary and h >= 1. For h < 1 the open chain has a nearly
zero edge mode and the ground state covariance is not well defined numerically.
"""
import numpy as np


def binary_entropy(p):
    """-p ln p - (1-p) ln(1-p), with 0 ln 0 = 0, without warnings."""
    p = np.asarray(p, dtype=float)
    out = np.zeros_like(p)
    m = (p > 0.0) & (p < 1.0)
    q = p[m]
    out[m] = -q * np.log(q) - (1.0 - q) * np.log1p(-q)
    return out


def majorana_matrix(L, h, J=1.0):
    A = np.zeros((2 * L, 2 * L))
    for j in range(L):
        A[2 * j, 2 * j + 1] = 2 * h
        A[2 * j + 1, 2 * j] = -2 * h
    for j in range(L - 1):
        A[2 * j + 1, 2 * j + 2] = 2 * J
        A[2 * j + 2, 2 * j + 1] = -2 * J
    return A


def ground_state_covariance(A):
    """Returns (Gamma, E0). Gamma[m,n] = (i/2) <[c_m, c_n]> in the ground state."""
    M = -A @ A
    w, W = np.linalg.eigh(M)
    w = np.clip(w, 0.0, None)
    # The rounding floor of sqrt(w.min()) is about 1e-8, so 1e-6 never depends on rounding.
    if np.sqrt(w.min()) < 1e-6:
        raise ValueError("near zero mode, free fermion ground state is not well defined here")
    inv_sqrt = (W / np.sqrt(w)) @ W.T
    Gamma = -A @ inv_sqrt
    E0 = -0.25 * np.sum(np.sqrt(w))
    return Gamma, E0


def block_entropy(Gamma, l):
    """Von Neumann entropy (natural log) of the first l sites (sites 0 to l-1)."""
    if l == 0 or 2 * l == Gamma.shape[0]:
        return 0.0
    G = Gamma[:2 * l, :2 * l]
    nu = np.abs(np.linalg.eigvalsh(1j * G))
    nu = np.clip(nu, 0.0, 1.0)
    p = (1.0 + nu) / 2.0
    return 0.5 * float(np.sum(binary_entropy(p)))


def entropy_profile_ff(Gamma, L):
    """S(l) for l = 0..L.

    The chain is symmetric under reflection and the ground state is unique, so S(l) = S(L - l)
    for the ground state. Compute only l <= L/2 and mirror.
    """
    S = np.zeros(L + 1)
    for l in range(1, L // 2 + 1):
        S[l] = block_entropy(Gamma, l)
        S[L - l] = S[l]
    return S
