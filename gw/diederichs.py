#!/usr/bin/env libtbx.python

import iotbx.merging_statistics
import sys
import numpy

hklins = sys.argv[1:]

merged_data = []

from scitbx import matrix

cc_matrix = matrix.sqr([0.0 for j in range(len(hklins) ** 2)]).as_numpy_array()

for j, hklin in enumerate(hklins):
    iobs = iotbx.merging_statistics.select_data(hklin, data_labels=None)
    iobs = iobs.customized_copy(anomalous_flag=True, info=iobs.info())
    merged = iobs.merge_equivalents().array()

    for k, m in enumerate(merged_data):
        cc = merged.correlation(m).coefficient()
        cc_matrix[k, j] = cc
        cc_matrix[j, k] = cc

    merged_data.append(merged)
    ms = iotbx.merging_statistics.dataset_statistics(i_obs=iobs)
    overall = ms.overall
    cc_star = overall.cc_star
    cc_matrix[j, j] = cc_star * cc_star

from numpy import linalg

from matplotlib import pyplot

pyplot.imshow(cc_matrix)
pyplot.show()

eigenvalues, eigenvectors = linalg.eig(cc_matrix)
