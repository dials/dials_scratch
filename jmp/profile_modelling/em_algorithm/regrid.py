def compute_centroid(a, b, p):
  from dials.array_family import flex
  from math import sqrt
  x = (a + b) / 2
  xc = flex.sum(p * x) / flex.sum(p)
  xv = flex.sum(p * (x - xc)**2) / (flex.sum(p)-1)
  return xc, sqrt(xv)


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
  i0 = 0
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

    i_list.append(s.data.all()[0] + i0)
    i0 = i0 + s.data.all()[0]

a_list = flex.double(a_list)
b_list = flex.double(b_list)
n_list = flex.double(n_list)
i_list = flex.size_t(i_list)


from math import pi, floor, ceil


from dials.algorithms.statistics import BinnedGMMSingle1DFixedMean, BinnedGMMSingle1D

from time import time



min_a = min(a_list)
max_b = max(b_list)
num_bins = 100
s1 = (num_bins-1) / (max_b - min_a)
s0 = -s1 * min_a

x1 = [min_a + i * (max_b-min_a) / num_bins for i in range(num_bins)]
x2 = [min_a + (i+1) * (max_b-min_a) / num_bins for i in range(num_bins)]
x1 = [xx * 180 / pi for xx in x1]
x2 = [xx * 180 / pi for xx in x2]
y = [0] * num_bins
for a, b, n in zip(a_list, b_list, n_list):
  i0 = s0 + s1 * a
  i1 = s0 + s1 * b
  assert i0 < i1
  for i in range(int(floor(i0)), int(ceil(i1))):
    j0 = max(i0, i)
    j1 = min(i1, i+1)
    if j0 < j1:
      fraction = (j1 - j0) / (i1 - i0)
      y[i] += n*fraction

mean, sigma =  compute_centroid(flex.double(x1), flex.double(x2), flex.double(y))


import pickle
pickle.dump((x1,x2,y), open(sys.argv[3], "w"))

# def gaussian(x, mu, sigma):
#   from math import pi, sqrt, exp
#   return [exp(-(xx - mu)**2 / (2*sigma**2)) for xx in x]

# from matplotlib import pylab
# pylab.plot(x1, y)
# pylab.legend()
# pylab.show()
