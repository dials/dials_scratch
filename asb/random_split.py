from __future__ import division
import sys, os, random
from libtbx import easy_pickle
from dials.array_family import flex

"""
Jiffy script to split a reflection table into two random halves.
Usage: libtbx.python random_split.py reflections.pickle
"""

filename = sys.argv[1]
name = os.path.basename(filename)
base, ext = os.path.splitext(name)
filea = base + "_a" + ext
fileb = base + "_b" + ext

data = easy_pickle.load(filename)

data_a = data.select(flex.random_permutation(len(data)))[:len(data)//2]
data_b = data.select(flex.random_permutation(len(data)))[len(data)//2:]
data_a = data_a.select(flex.sort_permutation(data_a['id']))
data_b = data_b.select(flex.sort_permutation(data_b['id']))

assert len(data_a) + len(data_b) == len(data)

easy_pickle.dump(filea, data_a)
easy_pickle.dump(fileb, data_b)

