import sys
from dxtbx.model.experiment_list import ExperimentListFactory

experiments = ExperimentListFactory.from_filenames(sys.argv[1:])
imageset = experiments[0].imageset
start, end = imageset.get_array_range()

for j in range(start, end):
    image = imageset.get_raw_data(j)[0]
    huge = 0x7FFFFFFF
    print(j, (image == huge).count(True))
