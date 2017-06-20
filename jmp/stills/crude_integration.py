
def integrate(experiment):
  from dials.algorithms.spot_prediction import PixelToMillerIndex
  from dials.array_family import flex
  from math import floor, sqrt
  from collections import defaultdict

  detector = experiment.detector
  assert len(detector) == 1
  panel = detector[0]

  xsize, ysize = panel.get_image_size()

  transform = PixelToMillerIndex(
    experiment.beam,
    experiment.detector,
    experiment.crystal)

  data = experiment.imageset.get_raw_data(0)[0]

  mask = flex.bool(flex.grid(ysize, xsize), False)
  reflections = defaultdict(list)

  print "Doing pixel labelling"
  for j in range(ysize):
    for i in range(xsize):
      h = transform.h(0, i, j)
      h0 = tuple(map(lambda x: int(floor(x+0.5)), h))

      d = sqrt(sum(map(lambda x,y: (x-y)**2, h, h0)))
      # if not hasattr(reflections[h0], "xd"):
      #   reflections[h0].xd = d
      #   reflections[h0].xc = i
      #   reflections[h0].yc = j
      # elif reflections[h0].xd > d:
      #   reflections[h0].xd = d
      #   reflections[h0].xc = i
      #   reflections[h0].yc = j
      
      if d < 0.3:
        mask[j,i] = True
      reflections[h0].append((j,i))

#   from matplotlib import pylab
#   pylab.imshow(mask.as_numpy_array(), interpolation='none')
#   pylab.show()

  print "Integrating reflections"
  miller_index = flex.miller_index()
  intensity = flex.double()
  variance = flex.double()
  bbox = flex.int6()
  xyz = flex.vec3_double()
  for h, r in reflections.iteritems():
    
    # xc = r.xc
    # yc = r.yc

    b_sum = 0
    f_sum = 0
    b_cnt = 0
    f_cnt = 0
    for i in range(len(r)):
      y, x = r[i]
      m = mask[y,x]
      if data[y,x] >= 0:
        if m:
          f_sum += data[y,x]
          f_cnt += 1
        else:
          b_sum += data[y,x]
          b_cnt += 1
    Y, X = zip(*r)
    x0, x1, y0, y1 = min(X), max(X), min(Y), max(Y)
    if f_cnt > 0 and b_cnt > 0:
      B = b_sum / b_cnt
      I = f_sum - B * f_cnt
      V = f_sum + B * (1 + f_cnt / b_cnt)
      miller_index.append(h)
      intensity.append(I)
      variance.append(V)
      bbox.append((x0, x1, y0, y1, 0, 1))
      # xyz.append((xc, yc, 0))
  
  print "Integrated %d reflections" % len(reflections)
  print flex.min(intensity), flex.max(intensity), flex.mean(intensity)
  reflections = flex.reflection_table()
  reflections["miller_index"] = miller_index
  reflections["intensity.sum.value"] = intensity
  reflections["intensity.sum.variance"] = variance
  reflections["bbox"] = bbox
  reflections["panel"] = flex.size_t(len(reflections), 0)
  reflections["id"] = flex.size_t(len(reflections), 0)
  # reflections["xyzcal.px"] = xyz
  # reflections["xyzobs.px"] = xyz
  reflections.set_flags(flex.size_t(range(len(reflections))),
                        reflections.flags.integrated_sum)
  return reflections




if __name__ == '__main__':

  import sys
  from dxtbx.model.experiment_list import ExperimentListFactory

  experiments_filename = sys.argv[1]
  reflections_filename = sys.argv[2]

  assert experiments_filename != ""
  assert reflections_filename != ""

  experiments = ExperimentListFactory.from_json_file(experiments_filename)
  assert len(experiments) == 1
  experiment = experiments[0]

  reflections = integrate(experiment)

  reflections.as_pickle(reflections_filename)

