import numpy as np
import pytest

from tfim_ed import ground_state
from entropy import entropy_profile
import json
import csv

from free_fermion import majorana_matrix, majorana_matrix_periodic, ground_state_covariance, block_entropy
from central_charge import c_from_two_sizes


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


@pytest.mark.parametrize("L,h", [(1024, 0.5), (256, 0.9)])
def test_near_zero_mode_raises_at_large_L(L, h):
    with pytest.raises(ValueError):
        ground_state_covariance(majorana_matrix(L, h))


def test_critical_chain_does_not_raise():
    Gamma, E0 = ground_state_covariance(majorana_matrix(512, 1.0))
    assert np.isfinite(E0)


@pytest.mark.parametrize("L", [4, 6, 8, 10, 12])
@pytest.mark.parametrize("h", [0.5, 1.0, 2.0])
def test_periodic_matches_ed(L, h):
    E_ed, psi = ground_state(L, h, "periodic")
    Gamma, E_ff = ground_state_covariance(majorana_matrix_periodic(L, h))
    assert abs(E_ed - E_ff) < 1e-9
    S_ed = entropy_profile(psi, L)
    for l in range(L + 1):
        assert abs(S_ed[l] - block_entropy(Gamma, l)) < 1e-8


@pytest.mark.parametrize("h", [0.5, 1.0, 2.0])
def test_periodic_energy_matches_exact_formula(h):
    L = 64
    _, E0 = ground_state_covariance(majorana_matrix_periodic(L, h))
    k = (2 * np.arange(L) + 1) * np.pi / L
    exact = -np.sum(np.sqrt(1.0 + h * h - 2.0 * h * np.cos(k)))
    assert abs(E0 - exact) < 1e-9


def test_periodic_guard_does_not_fire_in_the_ordered_phase():
    Gamma, E0 = ground_state_covariance(majorana_matrix_periodic(1024, 0.5))
    assert np.isfinite(E0)


def test_two_size_estimator_on_the_committed_open_chain_data():
    with open("results/ff_half_chain.csv") as f:
        rows = {int(r["L"]): float(r["S_half"]) for r in csv.DictReader(f) if float(r["h"]) == 1.0}
    with open("results/central_charge.json") as f:
        c_final = json.load(f)["c_final"]
    assert abs(c_from_two_sizes(rows[1024], 1024, rows[2048], 2048, "open") - c_final) < 1e-12
