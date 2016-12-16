

if __name__ == '__main__':

  from dials.array_family import flex
  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  import sys

  reflections = flex.reflection_table.from_pickle(sys.argv[1])
  experiments = ExperimentListFactory.from_json_file(sys.argv[2],
                                                     check_format=False)

  beam = experiments[0].beam
  detector = experiments[0].detector
  goniometer = experiments[0].goniometer
  scan = experiments[0].scan

  selection = reflections.get_flags(reflections.flags.used_in_refinement)
  reflections = reflections.select(selection)
  reflections.compute_zeta(experiments[0])

  phi = reflections['xyzcal.mm'].parts()[2]
  sbox = reflections['shoebox']
  zeta = reflections['zeta']

  print "Num Refl: ", len(reflections)

  a_list = []
  b_list = []
  n_list = []
  i_list = []
  for p, s, z in zip(phi, sbox, zeta):
    z0 = s.bbox[4]

    for k in range(s.data.all()[0]):
      phi0 = scan.get_angle_from_array_index(z0+k, deg=False)
      phi1 = scan.get_angle_from_array_index(z0+k+1, deg=False)

      a = (p - phi0) * z#abs(z)
      b = (p - phi1) * z#abs(z)
      if a > b:
        a, b = b, a

      sum_frames = 0
      for j in range(s.data.all()[1]):
        for i in range(s.data.all()[2]):
          if s.mask[k,j,i] != 0:
            sum_frames += s.data[k,j,i]
            #n_list.append(s.data[k,j,i])
      n_list.append(sum_frames)
      a_list.append(a)
      b_list.append(b)

    i_list.append(s.data.all()[0])

a_list = flex.double(a_list)
b_list = flex.double(b_list)
n_list = flex.double(n_list)
i_list = flex.size_t(i_list)


from math import floor, ceil, pi
a_list = a_list * 180.0 / pi
b_list = b_list * 180.0 / pi


import pickle

pickle.dump((a_list, b_list, n_list, i_list), open("extract.pickle", "w"))

# from matplotlib import pylab
# pylab.plot(X, Y)
# pylab.show()
