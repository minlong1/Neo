"""Unit tests for INI parsing (enabled by the Phase 2 refactor)."""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from conftest import load_case_config  # noqa: E402


def test_single_peak_case_values():
    cfg = load_case_config("voigt_baseline_svsc")
    assert cfg["npaths"] == 1
    assert cfg["BE_guess"] == [421.2]
    assert cfg["amp_guess"] == [1000.0]
    assert cfg["peak_type"] == ["Voigt"]
    assert cfg["background_type"] == ["Baseline", "SVSC"]
    assert cfg["is_singlet"] == [True]
    assert cfg["is_coster_kronig"] == [False]
    assert cfg["asym_limited"] == [True]
    assert cfg["gen_alt"] == 3
    assert cfg["size_population"] == 200
    assert cfg["number_of_generation"] == 5
    assert cfg["steady_state"] is False  # optional key, default applies
    assert cfg["k_range"] == [0.0, 0.1, 0.001]


def test_two_peak_case_per_peak_lists():
    cfg = load_case_config("two_peak_limited_correlated")
    assert cfg["npaths"] == 2
    assert cfg["BE_guess"] == [421.2, 434.2]
    assert cfg["BE_limited"] == [False, True]
    assert cfg["amp_guess"] == [1000.0, 850.0]
    # sigma of peak 2 is correlated to peak 1
    assert cfg["sigma_correlated"] == ["1", " 1"]
    # every per-peak list has one entry per peak
    for key in (
        "BE_guess",
        "sigma_guess",
        "gamma_guess",
        "amp_guess",
        "asym_guess",
        "sos_guess",
        "br_guess",
        "BE_limited",
        "sigma_limited",
        "gamma_limited",
        "amp_limited",
        "is_singlet",
        "is_coster_kronig",
    ):
        assert len(cfg[key]) == 2, key


def test_config_lists_are_shared_mutable_objects():
    """addPeak/removePeak mutate the config lists in place; the parser must
    hand out real lists, not copies-per-access (design.md Phase 2)."""
    cfg = load_case_config("voigt_baseline_svsc")
    be = cfg["BE_guess"]
    be.append(999.0)
    assert cfg["BE_guess"] is be
    assert cfg["BE_guess"][-1] == 999.0


def test_doublet_flags_parse_to_bools():
    cfg = load_case_config("double_lorentzian_doublet")
    assert cfg["is_singlet"] == [False]
    assert cfg["sos_guess"] == [13.0]
    assert cfg["br_guess"] == [0.57]
