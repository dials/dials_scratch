from __future__ import division

import libtbx.phil

phil_scope= libtbx.phil.parse("""
data = *raw corrected
  .type = choice
""")



def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks

  parser = OptionParser(
    read_datablocks=True,
    read_datablocks_from_images=True,
    phil=phil_scope,
    check_format=True)

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)
  assert len(datablocks) == 1
  imagesets = datablocks[0].extract_imagesets()

  img_count = 0
  import time
  t0 = time.time()
  for imgset in imagesets:
    for i in range(len(imgset)):
      if params.data == 'raw':
        imgset.get_raw_data(i)
      else:
        imgset.get_corrected_data(i)
      img_count += 1
      print "Read %i images" %img_count
  t1 = time.time()
  t = t1 - t0
  print "Read %i images in %.2fs (%.1f images/s)" %(
    img_count, t, img_count/t)

  return



if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  import sys
  run(sys.argv[1:])
