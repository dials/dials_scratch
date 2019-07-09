"""Test scoring of screw axes."""
from dials.array_family import flex
from dials_scratch.jbe.sys_abs.screw_axes import ScrewAxis41c, ScrewAxis42c,\
    ScrewAxis61c, ScrewAxis62c, ScrewAxis63c

def test_fourfold_screw():
    """Test the discernment between a 41 and 42 with weak pattern."""
    intensity_pattern_41 = flex.double([0.05, 0.03, 0.07, 10.0, 0.02, 0.04, 0.06, 12.0])
    weak_intensity_pattern_42 = flex.double([0.05, 2.0, 0.07, 10.0, 0.02, 1.4, 0.06, 12.0])

    reflection_table = flex.reflection_table()
    reflection_table['miller_index'] = flex.miller_index([
        (0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 0, 4),
        (0, 0, 5), (0, 0, 6), (0, 0, 7), (0, 0, 8),
    ])
    reflection_table['variance'] = flex.double(8, 1.0)

    # Test 41 and 42 screw
    reflection_table['intensity'] = intensity_pattern_41
    score_41 = ScrewAxis41c().score_axis(reflection_table)
    score_42 = ScrewAxis42c().score_axis(reflection_table)

    assert score_41 > 0.95
    assert score_42 < 0.05

    reflection_table['intensity'] = weak_intensity_pattern_42

    score_41 = ScrewAxis41c().score_axis(reflection_table)
    score_42 = ScrewAxis42c().score_axis(reflection_table)

    assert score_42 > score_41 # score_41 still fairly high as looks almost like 41
    assert score_42 > 0.8 #weak but should still be found

def test_sixfold_screw():
    """Test the discernment between a 61 and 62 with weak pattern."""
    intensity_pattern_61 = flex.double([0.05, 0.03, 0.07, 0.04, 0.02, 10.0,
        0.06, 0.12, 0.15, 0.02, 0.07, 12.0])
    weak_intensity_pattern_62 = flex.double([0.05, 0.03, 2.0, 0.04, 0.02, 10.0,
        0.06, 0.12, 1.5, 0.02, 0.07, 12.0])
    weak_intensity_pattern_63 = flex.double([0.05, 2.0, 0.07, 1.5, 0.02, 10.0,
        0.06, 1.2, 0.15, 1.7, 0.07, 12.0])
    intensity_pattern_62 = flex.double([0.05, 0.03, 5.0, 0.04, 0.02, 10.0,
        0.06, 0.12, 7.0, 0.02, 0.07, 12.0])

    reflection_table = flex.reflection_table()
    reflection_table['miller_index'] = flex.miller_index([
        (0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 0, 4),
        (0, 0, 5), (0, 0, 6), (0, 0, 7), (0, 0, 8),
        (0, 0, 9), (0, 0, 10), (0, 0, 11), (0, 0, 12),
    ])
    reflection_table['variance'] = flex.double(12, 1.0)

    # Test 61, 62 and 63 screw
    reflection_table['intensity'] = intensity_pattern_61
    score_61 = ScrewAxis61c().score_axis(reflection_table)
    score_62 = ScrewAxis62c().score_axis(reflection_table)
    score_63 = ScrewAxis63c().score_axis(reflection_table)

    assert score_61 > 0.92
    assert score_62 < 0.05
    assert score_63 < 0.05

    reflection_table['intensity'] = weak_intensity_pattern_62

    score_61 = ScrewAxis61c().score_axis(reflection_table)
    score_62 = ScrewAxis62c().score_axis(reflection_table)
    score_63 = ScrewAxis63c().score_axis(reflection_table)

    assert score_62 > score_61
    assert score_62 > score_63
    assert score_62 > 0.75

    reflection_table['intensity'] = intensity_pattern_62

    score_61 = ScrewAxis61c().score_axis(reflection_table)
    score_62 = ScrewAxis62c().score_axis(reflection_table)
    score_63 = ScrewAxis63c().score_axis(reflection_table)

    assert score_62 > score_61
    assert score_62 > score_63
    assert score_62 > 0.90

    reflection_table['intensity'] = weak_intensity_pattern_63

    score_61 = ScrewAxis61c().score_axis(reflection_table)
    score_62 = ScrewAxis62c().score_axis(reflection_table)
    score_63 = ScrewAxis63c().score_axis(reflection_table)

    assert score_63 > score_61
    assert score_63 > score_62
    assert score_63 > 0.9
