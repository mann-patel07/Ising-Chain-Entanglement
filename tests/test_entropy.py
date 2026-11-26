import numpy as np
import pytest

from tfim_ed import ground_state

from entropy import entanglement_entropy, schmidt_values, entropy_profile


def test_product_and_bell():
    psi = np.zeros(16)
    psi[0] = 1.0
    assert entanglement_entropy(psi, 4, 2) == 0.0
    bell = np.array([1, 0, 0, 1]) / np.sqrt(2)
    assert abs(entanglement_entropy(bell, 2, 1) - np.log(2)) < 1e-12


def test_random_state_reversal_and_rho():
    L = 8
    rng = np.random.default_rng(7)
    psi = rng.standard_normal(2**L) + 1j * rng.standard_normal(2**L)
    psi /= np.linalg.norm(psi)
    rev = psi.reshape([2] * L).transpose(range(L - 1, -1, -1)).reshape(-1)
    for l in range(L + 1):
        assert abs(entanglement_entropy(rev, L, L - l) - entanglement_entropy(psi, L, l)) < 1e-10
    prof = entropy_profile(psi, L)

    def ent(ev):
        ev = ev[ev > 1e-16]
        return -np.sum(ev * np.log(ev))

    for l in range(1, L):
        M = psi.reshape(2**l, 2**(L - l))
        s1 = ent(np.linalg.eigvalsh(M @ M.conj().T))
        s2 = ent(np.linalg.eigvalsh(M.conj().T @ M))
        assert abs(s1 - s2) < 1e-10
        assert abs(s1 - prof[l]) < 1e-10
    sv = schmidt_values(psi, L, 4)
    assert np.all(np.diff(sv) <= 1e-15) and abs(sv.sum() - 1) < 1e-12


@pytest.mark.parametrize("boundary", ["open", "periodic"])
@pytest.mark.parametrize("h", [0.5, 1.0, 2.0])
def test_ground_state_reflection_symmetry(boundary, h):
    L = 10
    _, psi = ground_state(L, h, boundary)
    S = entropy_profile(psi, L)
    for l in range(L + 1):
        assert abs(S[l] - S[L - l]) < 1e-10


def test_bit_order():
    L = 3
    psi = np.zeros(2**L)
    psi[4] = 1.0
    # site 0 is the most significant bit, so state 4 = 100 has only site 0 equal to 1
    assert (4 >> (L - 1 - 0)) & 1 == 1
    assert (4 >> (L - 1 - 1)) & 1 == 0
