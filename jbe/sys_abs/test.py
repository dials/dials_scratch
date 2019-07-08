from dials.array_family import flex
from dials_scratch.jbe.sys_abs.laue_groups_info import laue_groups, score_screw_axes, score_space_groups
from dials.util.filter_reflections import ScaleIntensityReducer
from libtbx.table_utils import simple_table
import time

# load some data
r1 = flex.reflection_table.from_pickle(
    "/Users/whi10850/Documents/test_data/th_8_2/th_8_2_analysis/080719/scaled.pickle")
#r2 = flex.reflection_table.from_pickle(
#    "/Users/whi10850/Documents/example_data/x4-wide/xia2_dials_refactor/DEFAULT/scale/19_scaled_reflections.pickle")
r = ScaleIntensityReducer.reduce_on_intensities(r1)
print("Number of reflections in dataset: %s" % r.size())
r['intensity'] = r['intensity.scale.value']
r['variance'] = r['intensity.scale.variance']

# Get the laue class from the space group
lg = 'P 4/m m m'

start = time.time()

screw_axis_scores = score_screw_axes(laue_groups[lg])
score_space_groups(screw_axis_scores, laue_groups[lg])

end = time.time()
print(end-start)
