"""Exact diagonalization of the transverse field Ising chain, H = -J sum X X - h sum Z.

Basis convention: a basis state is an integer n in [0, 2^L). The bit of site i is
b_i = (n >> (L - 1 - i)) & 1, so site 0 is the most significant bit.
"""
import numpy as np
import scipy.sparse
import scipy.sparse.linalg


def sector_states(L, parity):
    """Sorted int64 array of all n in [0, 2^L) whose number of 1 bits has the given parity (0 even, 1 odd)."""
    n = np.arange(2**L, dtype=np.int64)
    count = np.zeros_like(n)
    for k in range(L):
        count += (n >> k) & 1
    return n[(count % 2) == parity]


def build_hamiltonian(L, h, boundary, parity, J=1.0):
    """Sparse CSR Hamiltonian restricted to one parity sector. Returns (H, states)."""
    if boundary not in ("open", "periodic"):
        raise ValueError("boundary must be 'open' or 'periodic', got %r" % (boundary,))
    states = sector_states(L, parity)
    D = len(states)
    index = -np.ones(2**L, dtype=np.int64)
    index[states] = np.arange(D)
    shifts = (L - 1 - np.arange(L))
    bits = (states[:, None] >> shifts[None, :]) & 1
    diag = -h * np.sum(1 - 2 * bits, axis=1).astype(float)
    bonds = [(i, i + 1) for i in range(L - 1)]
    if boundary == "periodic":
        bonds.append((L - 1, 0))
    rows = [np.arange(D)]
    cols = [np.arange(D)]
    vals = [diag]
    for (i, j) in bonds:
        mask = (1 << (L - 1 - i)) | (1 << (L - 1 - j))
        rows.append(np.arange(D))
        cols.append(index[states ^ mask])
        vals.append(np.full(D, -J))
    H = scipy.sparse.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(D, D)
    ).tocsr()
    return H, states


def _lowest(L, h, boundary, parity, J=1.0):
    H, states = build_hamiltonian(L, h, boundary, parity, J)
    D = len(states)
    if D <= 4096:
        w, v = np.linalg.eigh(H.toarray())
        return float(w[0]), v[:, 0], states
    w, v = scipy.sparse.linalg.eigsh(
        H, k=1, which="SA", tol=1e-12, v0=np.random.default_rng(12345).standard_normal(D)
    )
    return float(w[0]), v[:, 0], states


def ground_state(L, h, boundary, parity=0, J=1.0):
    """Returns (E0, psi_full). psi_full has length 2^L, is zero outside the sector and has norm 1."""
    E0, vec, states = _lowest(L, h, boundary, parity, J)
    psi = np.zeros(2**L)
    psi[states] = vec
    psi /= np.linalg.norm(psi)
    if psi[np.argmax(np.abs(psi))] < 0:
        psi = -psi
    return E0, psi


def sector_energies(L, h, boundary, J=1.0):
    """Returns (E_even, E_odd), the lowest energy in each parity sector."""
    return (_lowest(L, h, boundary, 0, J)[0], _lowest(L, h, boundary, 1, J)[0])
