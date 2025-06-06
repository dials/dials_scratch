"""Test components of the systematic_absences program."""
import os
import pytest
import procrunner
from dxtbx.serialize import load
from dials.array_family import flex
from dials_scratch.jbe.sys_abs.laue_groups_info import (
    score_space_groups,
    laue_groups,
    score_screw_axes,
)


def test_systematic_absences_script(dials_data, run_in_tmpdir):
    """Test the command line script with real data. Proteinase K in P41"""
    location = dials_data("vmxi_proteinase_k_sweeps")

    command = [
        "dials_scratch.systematic_absences",
        "output.experiments=symmetrized.expt",
    ]
    command.append(location.join("experiments_0.json").strpath)
    command.append(location.join("reflections_0.pickle").strpath)

    result = procrunner.run(command)
    assert result["exitcode"] == 0
    assert result["stderr"] == ""
    assert os.path.exists("absences.html")
    assert os.path.exists("symmetrized.expt")
    exps = load.experiment_list("symmetrized.expt", check_format=False)
    assert str(exps[0].crystal.get_space_group().info()) == "P 41"

    # Now try with a d_min
    command += ["d_min=4.0"]
    result = procrunner.run(command)
    assert result["exitcode"] == 0
    assert result["stderr"] == ""
    exps = load.experiment_list("symmetrized.expt", check_format=False)
    assert str(exps[0].crystal.get_space_group().info()) == "P 41"


def test_score_screw_axes_equivalent_axes():
    """Test the score_screw_axes function (when equivalent axes present)."""
    laue_group_info = laue_groups["P m -3"]

    # unique axis 21a, but should pick up 21b, 21c
    reflections = flex.reflection_table()
    reflections["miller_index"] = flex.miller_index(
        [(0, 1, 0), (0, 2, 0), (0, 0, 1), (0, 0, 2)]
    )
    reflections["intensity"] = flex.double([0.05, 100.0, 0.02, 100.0])
    reflections["variance"] = flex.double([1.0, 1.0, 1.0, 1.0])

    axes, scores = score_screw_axes(laue_group_info, reflections)
    assert len(scores) == 1
    assert len(axes) == 1
    assert axes[0].name == "21a"
    assert scores[0] > 0.99


def test_score_space_group():
    """Test scoring of space groups by combining axis scores."""

    # simple case, P 1 21 2, single score
    laue_group = laue_groups["P 1 2/m 1"]
    axis_scores = [0.98]  # 21 > P41212
    space_groups, scores = score_space_groups(axis_scores, laue_group)
    for sg, score in zip(space_groups, scores):
        if sg == "P 21":
            assert score == pytest.approx(0.98)
        elif sg == "P 2":
            assert score == pytest.approx(0.02)

    # More complex case - only use 42 score in space groups where not 41
    laue_group = laue_groups["P 4/m m m"]
    axis_scores = [0.95, 1.0, 0.95]  # 41, 21 and 42 score high > P41212
    space_groups, scores = score_space_groups(axis_scores, laue_group)
    for sg, score in zip(space_groups, scores):
        if sg == "P 41 21 2":
            assert score == pytest.approx(0.95)
        elif sg == "P 42 21 2":
            assert score == pytest.approx(0.05 * 0.95)
        elif sg == "P 4 21 2":
            assert score == pytest.approx(0.05 * 0.05)
        else:
            assert score == pytest.approx(0.0)

    # More complex case - only use 62, 63 score in space groups where not 61
    laue_group = laue_groups["P 6/m"]
    axis_scores = [0.95, 0.9, 0.85]  # 61, 62, 63
    space_groups, scores = score_space_groups(axis_scores, laue_group)
    for sg, score in zip(space_groups, scores):
        if sg == "P 61":
            assert score == pytest.approx(0.95)
        elif sg == "P 62":
            assert score == pytest.approx(0.05 * 0.9)
        elif sg == "P 63":
            assert score == pytest.approx(0.05 * 0.85)
        elif sg == "P 6":
            assert score == pytest.approx(0.05 * 0.1 * 0.15)
