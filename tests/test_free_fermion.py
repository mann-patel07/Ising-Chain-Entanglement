import numpy as np
import pytest

from tfim_ed import ground_state
from entropy import entropy_profile
from free_fermion import majorana_matrix, ground_state_covariance, block_entropy


@pytest.mark.parametrize("L", [4, 6, 8, 10, 12])
@pytest.mark.parametrize("h", [1.0, 1.5, 2.0])
def test_matches_ed(L, h):
    E_ed, psi = ground_state(L, h, "open")
    Gamma, E_ff = ground_state_covariance(majorana_matrix(L, h))
    assert abs(E_ed - E_ff) < 1e-9
    S_ed = entropy_profile(psi, L)
    for l in range(L + 1):
        assert abs(S_ed[l] - block_entropy(Gamma, l)) < 1e-8


def test_covariance_properties():
    Gamma, _ = ground_state_covariance(majorana_matrix(50, 1.0))
    assert np.max(np.abs(Gamma + Gamma.T)) < 1e-12
    assert np.max(np.abs(Gamma @ Gamma + np.eye(100))) < 1e-8


def test_single_site():
    _, E0 = ground_state_covariance(majorana_matrix(1, 0.7))
    assert abs(E0 + 0.7) < 1e-12


def test_near_zero_mode():
    with pytest.raises(ValueError):
        ground_state_covariance(majorana_matrix(100, 0.5))
