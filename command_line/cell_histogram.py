from __future__ import division, print_function

def load_experiment(experiment_file):
  from dxtbx.model.experiment_list import ExperimentListFactory
  return ExperimentListFactory.from_json_file(experiment_file,
                                              check_format=False) 

def cell_hist():
  from dials.array_family import flex
  import sys
  
  a = flex.double()
  b = flex.double()
  c = flex.double()
  al = flex.double()
  be = flex.double()
  ga = flex.double()

  for arg in sys.argv[1:]:
    expt = load_experiment(arg)
    for xtal in expt.crystals():
      cell = xtal.get_unit_cell().parameters()
      a.append(cell[0])
      b.append(cell[1])
      c.append(cell[2])
      al.append(cell[3])
      be.append(cell[4])
      ga.append(cell[5])

  a_h = flex.histogram(a, data_min=0, data_max=100, n_slots=1000)
  b_h = flex.histogram(b, data_min=0, data_max=100, n_slots=1000)
  c_h = flex.histogram(c, data_min=0, data_max=100, n_slots=1000)
  al_h = flex.histogram(al, data_min=0, data_max=100, n_slots=1000)
  be_h = flex.histogram(be, data_min=0, data_max=100, n_slots=1000)
  ga_h = flex.histogram(ga, data_min=0, data_max=100, n_slots=1000)

  for v in zip(a_h.slot_centers(), a_h.slots(), b_h.slots(), c_h.slots(),
               al_h.slots(), be_h.slots(), ga_h.slots()):
    print('%5.2f %5d %5d %5d %5d %5d %5d' % v)


if __name__ == '__main__':
  cell_hist()
