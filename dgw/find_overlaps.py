#!/bin/env dials.python
"""Find overlapped reflections from a scaled data set, where the first
two lattices come from a single set of images and contain overlaps"""

from dxtbx.model.experiment_list import ExperimentList
import numpy as np

from annlib_ext import AnnAdaptor
from scitbx.array_family import flex

from dials.array_family import flex
from dials.algorithms.scaling.outlier_rejection import reject_outliers
from dials.algorithms.scaling.Ih_table import IhTable

from matplotlib import pyplot as plt


def moving_average(x, w):
    return np.convolve(x, np.ones(w), "valid") / w


# Open the scaled data and set required columns
expt = ExperimentList.from_file("scaled.expt")
refl = flex.reflection_table.from_file("scaled.refl")
refl["intensity"] = refl["intensity.scale.value"]
refl["variance"] = refl["intensity.scale.variance"]

# Remove bad reflections
refl = refl.select(refl["variance"] > 0)
refl = refl.select(refl["inverse_scale_factor"] > 0)

# Determine merged intensities
space_group = expt[0].crystal.get_space_group()
ih = IhTable([refl], space_group, nblocks=1)
ih = ih.blocked_data_list[0].as_reflection_table()
ih.sort("loc_indices")
refl["Ih"] = ih["Ih_values"]

# Remove reflections with weak merged intensities
refl = refl.select(refl["Ih"] > 0.05)

# Determine maximum z-score for which each reflection is still an outlier
refl["z_score"] = flex.double(len(refl), 0)
for z in np.arange(0.01, 7, 0.01):
    print(z)
    refl = reject_outliers(refl, expt[0], zmax=z)
    sel = refl.get_flags(refl.flags.outlier_in_scaling)
    refl["z_score"].set_selected(sel, z)

# Select just the reflections from the first two crystal models
r1 = refl.select(refl["id"] == 0)
r2 = refl.select(refl["id"] == 1)

# key = "xyzobs.px.value"
key = "xyzcal.px"

# Create the KD Tree and find the nearest neighbours
ann = AnnAdaptor(r2[key].as_double(), dim=3, k=1)
ann.query(r1[key].as_double())

# Select only close neighbours that might be overlaps. Average sigma_b
# is about 0.0095°. Detector distance is 644.4 mm. So, sigma_b covers
# about a single 0.1 mm pixel. Consider overlaps to be close when
# centroids are separated by 3 sigma_b, i.e. 3 pixels
r1["distances"] = flex.sqrt(ann.distances)
sel = r1["distances"] < 6
r1 = r1.select(sel)
r1.sort("distances")

# Plot z_score vs distances
plt.scatter(r1["distances"], r1["z_score"], marker=".")
plt.axhline(y=6.0, color="r", linestyle="-")
plt.xlabel("Distance (pixel space)")
plt.ylabel("Z-score")
plt.title("Outlier rejection of overlaps")
plt.savefig("z-score.pdf")

plt.clf()
scaled_intensity = r1["intensity"] / r1["inverse_scale_factor"]
ratio = scaled_intensity / r1["Ih"]
ma = moving_average(ratio, 60)
plt.scatter(r1["distances"], ratio, marker=".")
plt.plot(r1["distances"][30 : (30 + len(ma))], ma, color="black")
plt.axhline(y=1.0, color="r", linestyle="-")
plt.xlabel("Distance (pixel space)")
plt.ylabel("I / Ih")
plt.ylim(-2, 6)
plt.title("Intensity inflation of overlaps")
plt.savefig("inflation.pdf")
