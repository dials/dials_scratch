from dials.array_family import flex

# use case - a list of single reflection tables and experiments are given, want to
# produce a list of single reflection tables and experiments with unique identifiers
# - also want to have unique dataset id for each reflection table so that I can
# easily create a combined reflection table by extending existing ones.


def test_extend_reflection_table():
    """Simple case, two reflection tablss with different ids and key."""

    rt = flex.reflection_table()
    rt["id"] = flex.int([0])
    rt.experiment_identifiers()[0] = "a"
    # rt.experiment_identifiers()[1] = 'b' < passes if these two lines are included.
    # rt.experiment_identifiers()[2] = 'c' < Should this be necessary
    assert rt.are_experiment_identifiers_consistent()

    rt2 = flex.reflection_table()
    rt2["id"] = flex.int([1, 2])
    rt2.experiment_identifiers()[1] = "b"
    rt2.experiment_identifiers()[2] = "c"
    assert rt2.are_experiment_identifiers_consistent()

    rt.extend(
        rt2
    )  # currently fails here - "Experiment identifiers do not match" - bug?
    assert rt.are_experiment_identifiers_consistent()

    assert list(rt.experiment_identifiers().keys()) == [0, 1, 2]
    assert list(rt.experiment_identifiers().values()) == ["a", "b", "c"]


# For the test below, this is more a question of what behaviour we want.
# The behaviour below is not necessary if one preprepares the datasets to have
# unique identifiers an values for ['id'] in the reflection tables. This is
# what I have decided to do for now, in
# dials.algorithms.scaling.scaling_utilities.assign_unique_identifiers


def test_conflicting_extend_reflection_table():
    """Simple case, two experiments, but where they both have conflicting ids.
    Should this be up to the user to ensure that the ids are unique before, are
    should it be assumed that if the 'id' is the same then this should be the same
    experiment?"""

    rt = flex.reflection_table()
    rt["id"] = flex.int([0])
    rt.experiment_identifiers()[0] = "a"
    assert rt.are_experiment_identifiers_consistent()

    rt2 = flex.reflection_table()
    rt2["id"] = flex.int([0])
    rt2.experiment_identifiers()[0] = "b"
    assert rt2.are_experiment_identifiers_consistent()

    rt.extend(rt2)
    # would require internally updating ids to match new experiment identifiers

    assert list(rt["id"]) == [0, 1]  # currently gives [0, 0]
    assert list(rt.experiment_identifiers().keys()) == [0, 1]  # currently gives [0, 0]
    assert list(rt.experiment_identifiers().values()) == [
        "a",
        "b",
    ]  # currently gives ['b']
