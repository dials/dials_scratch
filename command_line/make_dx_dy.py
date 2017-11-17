from __future__ import division, print_function
from iotbx import phil

scope = phil.parse('''
  dx = 0.0
    .type = float
  dy = 0.0
    .type = float
''')

def main():
  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  import libtbx.load_env

  usage = "%s [options] image_*.cbf" % (
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=scope,
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

  image = imageset[images[0]]

  from dials.array_family import flex
  import cPickle as pickle

  dx = []
  dy = []

  for block in image:
    dx.append(flex.double(flex.grid(block.focus()), params.dx))
    dy.append(flex.double(flex.grid(block.focus()), params.dy))

  dx = tuple(dx)
  dy = tuple(dy)

  with open('dx.pickle', 'w') as f:
    pickle.dump(dx, f)

  with open('dy.pickle', 'w') as f:
    pickle.dump(dy, f)

if __name__ == '__main__':
  main()
