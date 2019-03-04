#!/usr/bin/env python

from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import cPickle as pickle
import enum
import dials.array_family.flex as flex
import sys

if "-h" in sys.argv or "--help" in sys.argv or len(sys.argv[1:]) != 1:
    print("Usage: plot_centroid_per_shoebox.py <indexed.pickle>")
    sys.exit(0)

with open(sys.argv[1]) as f:
    reflections = pickle.load(f)

used_in_refine = reflections.get_flags(reflections.flags.used_in_refinement)
refined = reflections.select(used_in_refine)
# not_refined = reflections.select(~used_in_refine)

plt.figure()
plt.hist([x.zsize() for x in reflections["shoebox"]], bins=8, range=(1, 9))
plt.xlabel("Shoebox width")
plt.ylabel("N")

plt.figure()
plt.suptitle("Centroid differences for reflections that span N images")
for n in range(1, 7):
    plt.subplot(2, 3, n)
    select = flex.bool([x.zsize() == n for x in refined["shoebox"]])
    refs = refined.select(select)

    _, _, zc = [x.as_numpy_array() for x in refs["xyzcal.px"].parts()]
    _, _, zo = [x.as_numpy_array() for x in refs["xyzobs.px.value"].parts()]
    zd = zo - zc
    plt.hist2d(zc % 1, zd, range=[[0, 1], [-1, 1]], bins=[30, 30])
    plt.xlim(0, 1)
    plt.ylim(-0.5, 0.5)
    plt.title(n)
    plt.xlabel("$Z_c \% 1$")
    plt.ylabel("$Z_o-Z_c$")
plt.tight_layout()
plt.subplots_adjust(top=0.89)

plt.show()
