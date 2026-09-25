"""Contract tests for config/run_config.json and the path convention."""
import pytest

from arsenic_hra import paths


def test_config_matches_frozen_decisions():
    cfg = paths.load_run_config()
    assert cfg["arsenic"]["selected_districts"] == [
        "Dhaka", "Chittagong", "Rajshahi", "Khulna", "Barisal",
        "Sylhet", "Rangpur", "Mymensingh", "Comilla", "Bogra",
    ]
    assert cfg["arsenic"]["candidate_distributions"] == ["normal", "lognormal", "gamma"]
    assert cfg["bodyweight"]["candidate_distributions"] == ["normal", "lognormal", "gamma", "triangular"]
    assert cfg["arsenic"]["censored_substitution"] == "none"
    assert cfg["arsenic"]["ugL_to_mgL_divisor"] == 1000.0
    sim = cfg["simulation"]
    assert sim["primary_N"] == 10000
    assert sim["convergence_N"] == [1000, 5000, 10000, 20000]
    assert isinstance(sim["master_seed"], int)
    assert sim["input_dependence"] == "independent"


def test_raw_inputs_exist_under_data_raw():
    for relative in paths.load_run_config()["raw_inputs"].values():
        assert paths.raw_path(relative).is_file(), relative


def test_raw_path_rejects_outside_data_raw():
    with pytest.raises(ValueError):
        paths.raw_path("results/arsenic/x.csv")


def test_output_dir_refuses_raw():
    with pytest.raises(ValueError):
        paths.ensure_output_dir(paths.DATA_RAW / "processed")
