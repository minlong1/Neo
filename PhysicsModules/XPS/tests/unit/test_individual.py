"""Round-trip tests for Individual get/set (enabled by Phase 2)."""

import pathlib
import random
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from PhysicsModules.XPS.xps_neo.xps_individual import Individual  # noqa: E402
from PhysicsModules.XPS.xps_neo.xps import XPS_GA  # noqa: E402


def one_peak_pars_range():
    """pars_range in the exact shape XPS_GA.initialize_range builds."""
    return {
        "Binding Energy": [[-0.2, 0.2, 0.01]],
        "BE": [421.2],
        "BE_limited": [False],
        "BE_correlated": [1],
        "BE_correlated_mult": [1],
        "Gaussian": [[-0.3, 0.3, 0.001]],
        "Sigma": [1.0],
        "sigma_limited": [False],
        "sigma_correlated": [1],
        "sigma_correlated_mult": [1],
        "Lorentzian": [[-0.2, 0.2, 0.001]],
        "Gamma": [0.5],
        "gamma_limited": [False],
        "gamma_correlated": [1],
        "gamma_correlated_mult": [1],
        "Amplitude": [[-200.0, 200.0, 0.05]],
        "Amp": [1000.0],
        "amp_limited": [False],
        "amp_correlated": [1],
        "amp_correlated_mult": [1],
        "Asymmetry": [[-0.00001, 0.00001, 0.00001]],
        "Asym": [1.0],
        "asym_limited": [True],
        "asym_correlated": [1],
        "asym_correlated_mult": [1],
        "Asymmetry Doniach-Sunjic": [0.0, 1.0, 0.01],
        "k_range": [0.01, 0.1, 0.001],
        "Background": [0.0, 5000.0, 0.05],
        "Slope": [0.0, 200.0, 0.01],
        "npeaks": 1,
        "baseline": [0.0, 100.0, 0.1],
        "br_range": [[-0.05, 0.05, 0.01]],
        "BR": [0.5],
        "br_limited": [False],
        "br_correlated": [1],
        "br_correlated_mult": [1],
        "is_singlet": [True],
        "is_coster_kronig": [False],
        "Lorentzian Coster-Kronig": [0.1, 1.5, 0.001],
        "Gamma Coster-Kronig": [0.5],
        "Asymmetry Range Coster-Kronig": [0.0, 10.0, 0.01],
        "Asymmetry Coster-Kronig": [1.0],
        "spinOrbitSplitting": [[-0.02, 0.02, 0.01]],
        "SOS": [0.0],
        "sos_limited": [False],
        "sos_correlated": [1],
        "sos_correlated_mult": [1],
        "photoline": "4d",
    }


def make_individual(seed):
    random.seed(seed)
    np.random.seed(seed)
    return Individual(["Baseline"], ["Voigt"], False, one_peak_pars_range())


def test_get_params_roundtrip_through_setpars():
    """get_params -> split_into_x -> setPars reproduces the peak params
    (the exact path the GA's DE crossover uses)."""
    a = make_individual(seed=1)
    b = make_individual(seed=2)
    params_a = a.get_params()
    assert params_a != b.get_params(), "different seeds should differ"

    XPS_GA.setPars(b, XPS_GA.split_into_x(list(params_a)))

    peaks_a = XPS_GA.split_into_x(list(a.get_params()))[:-1]  # bg group
    peaks_b = XPS_GA.split_into_x(list(b.get_params()))[:-1]  # is not set
    assert peaks_b == peaks_a


def test_len_matches_param_count():
    ind = make_individual(seed=3)
    assert len(ind) == len(ind.get_params())


def test_individual_construction_is_seed_deterministic():
    p1 = make_individual(seed=7).get_params()
    p2 = make_individual(seed=7).get_params()
    assert p1 == p2
