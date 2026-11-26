import numpy as np
import pytest

from tfim_ed import build_hamiltonian, ground_state, sector_energies, sector_states
from entropy import entanglement_entropy

X = np.array([[0, 1], [1, 0]], dtype=float)
Z = np.array([[1, 0], [0, -1]], dtype=float)
I2 = np.eye(2)


def op(L, i, P):
    mats = [I2] * L
    mats[i] = P
    out = np.array([[1.0]])
    for m in mats:
        out = np.kron(out, m)
    return out


def full_h(L, h, boundary):
    H = np.zeros((2**L, 2**L))
    bonds = [(i, i + 1) for i in range(L - 1)]
    if boundary == "periodic":
        bonds.append((L - 1, 0))
    for i, j in bonds:
        H -= op(L, i, X) @ op(L, j, X)
    for i in range(L):
        H -= h * op(L, i, Z)
    return H


@pytest.mark.parametrize("L", [2, 3, 4])
@pytest.mark.parametrize("boundary", ["open", "periodic"])
def test_sectors_match_kron(L, boundary):
    h = 0.83
    Hfull = full_h(L, h, boundary)
    Hrec = np.zeros_like(Hfull)
    for parity in (0, 1):
        H, states = build_hamiltonian(L, h, boundary, parity)
        Hrec[np.ix_(states, states)] = H.toarray()
    assert np.max(np.abs(Hrec - Hfull)) < 1e-12


def test_symmetric():
    for boundary in ("open", "periodic"):
        H, _ = build_hamiltonian(7, 0.9, boundary, 0)
        assert abs(H - H.T).max() < 1e-14


@pytest.mark.parametrize("L", [4, 6, 8, 10, 12])
@pytest.mark.parametrize("h", [0.5, 1.0, 2.0])
def test_periodic_energy_formula(L, h):
    k = (2 * np.arange(L) + 1) * np.pi / L
    exact = -np.sum(np.sqrt(1 + h**2 - 2 * h * np.cos(k)))
    E_even, _ = sector_energies(L, h, "periodic")
    assert abs(E_even - exact) < 1e-9
    if h == 1.0:
        assert abs(E_even + 2 / np.sin(np.pi / (2 * L))) < 1e-9


@pytest.mark.parametrize("L", [4, 5, 6, 7, 8])
@pytest.mark.parametrize("h", [0.3, 0.7, 1.0, 1.5])
@pytest.mark.parametrize("boundary", ["open", "periodic"])
def test_even_sector_lowest(L, h, boundary):
    E_even, E_odd = sector_energies(L, h, boundary)
    assert E_even <= E_odd + 1e-10


@pytest.mark.parametrize("L", [6, 8])
@pytest.mark.parametrize("boundary", ["open", "periodic"])
def test_h0_cat_state(L, boundary):
    _, psi = ground_state(L, 0.0, boundary)
    for l in range(1, L):
        assert abs(entanglement_entropy(psi, L, l) - np.log(2)) < 1e-9


def test_large_field_product():
    _, psi = ground_state(8, 1000.0, "open")
    assert entanglement_entropy(psi, 8, 4) < 1e-4


def test_bad_boundary():
    with pytest.raises(ValueError):
        build_hamiltonian(4, 1.0, "closed", 0)


def test_sector_states():
    s = sector_states(3, 0)
    assert list(s) == [0, 3, 5, 6]
