from __future__ import division, absolute_import, print_function
from dxtbx.format import setup_hdf5_plugin_path

setup_hdf5_plugin_path()

import h5py
import numpy
import sys

from scitbx.array_family import flex

f = h5py.File(sys.argv[1])
dataset = f[sys.argv[2]]

depth, height, width = dataset.shape

for j in range(depth):
    frame = flex.int(numpy.int32(dataset[j, :, :])).as_1d()
    print(j, flex.sum(frame.select(frame < 65535)))
