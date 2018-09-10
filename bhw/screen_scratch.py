#!/usr/bin/env python
# coding=utf-8

from dials.util.phil import ReflectionTableConverters
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex
from matplotlib import pyplot as plt

rtable = ReflectionTableConverters().from_string('indexed.pickle').data
elist = ExperimentListFactory.from_json_file('experiments.json')
unit_cell = elist.crystals()[0].get_unit_cell()
indexed = rtable.select(rtable.get_flags(rtable.flags.indexed))
d = unit_cell.d(indexed['miller_index'])

plt.plot(1/flex.pow2(d), indexed['intensity.sum.value'], '.')
plt.semilogy()
plt.show()