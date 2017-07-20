import os
from iotbx import mtz

m = mtz.object(file_name=os.path.join(os.environ["CEXAM"], "data", "insulin_unmerged.mtz"))
h = m.extract_miller_indices()

# get safe access to the crystal owned by m
x = mtz.crystal(m, 0)

# get safe access to the dataset owned by x
d = mtz.dataset(x, 0)

# add a column
c = d.add_column("my_column", "R")

# copy the intensity column
intensities = m.get_column("I").extract_values()

# put double the intensity in my_column
c.set_values(intensities * 2)

# write out the file
m.write("cctbx_output_insulin_unmerged.mtz")
