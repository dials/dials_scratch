def main():
  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  import libtbx.load_env

  usage = "%s [options] image_*.cbf" % (
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    read_datablocks=True,
    read_datablocks_from_images=True)

  params, options = parser.parse_args()
  datablocks = flatten_datablocks(params.input.datablock)

  if len(datablocks) == 0:
    parser.print_help()
    exit()

  assert(len(datablocks) == 1)

  datablock = datablocks[0]
  imagesets = datablock.extract_imagesets()

  assert(len(imagesets) == 1)

  imageset = imagesets[0]

  images = imageset.indices()

  previous = signal(imageset, images[0])

  # to do:
  # given grid dimensions, compute common-map across up / left / down / right
  # (is abelian)

  for indx in images[1:]:
    next = signal(imageset, indx)
    print indx, (next & previous).count(True)
    previous = next

def signal(imageset, indx):
  from dials.array_family import flex
  from libtbx.phil import parse

  detectors = imageset.get_detector()
  assert(len(detectors) == 1)
  detector = detectors[0]
  trusted = detector.get_trusted_range()

  data = imageset.get_raw_data(indx)
  assert(len(data) == 1)
  data = data[0]
  negative = (data < 0)
  hot = (data > int(round(trusted[1])))
  bad = negative | hot

  from dials.algorithms.spot_finding.factory import SpotFinderFactory
  from dials.algorithms.spot_finding.factory import phil_scope

  data = data.as_double()

  from dxtbx import datablock

  spot_params = phil_scope.fetch(source=parse("")).extract()
  threshold_function = SpotFinderFactory.configure_threshold(
    spot_params, datablock.DataBlock([imageset]))
  peak_pixels = threshold_function.compute_threshold(data, ~bad)
  return peak_pixels

if __name__ == '__main__':
  main()

