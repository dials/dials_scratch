from __future__ import division

import iotbx.phil

help_message = '''
dials.copy_metrology datablock.json reference=reference_datablock.json
'''

phil_scope = iotbx.phil.parse("""
max_delta_distance = 1
  .type = float(value_min=0)
input {
  reference = None
    .type = path
}
output {
  datablock = metrology_corrected_datablock.json
    .type = path
}
""", process_includes=True)


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  import libtbx.load_env

  usage = "%s [options] datablock.json reference=reference_datablock.json" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)

  if len(datablocks) == 0:
    parser.print_help()
    exit()

  # Load reference geometry
  reference_detector = None
  if params.input.reference is not None:
    from dxtbx.serialize import load
    try:
      reference_experiments = load.experiment_list(
        params.input.reference, check_format=False)
      assert len(reference_experiments.detectors()) == 1
      reference_detector = reference_experiments.detectors()[0]
    except Exception, e:
      reference_datablocks = load.datablock(params.input.reference)
      assert len(reference_datablocks) == 1
      imageset = reference_datablocks[0].extract_imagesets()[0]
      reference_detector = imageset.get_detector()

  assert len(datablocks) == 1
  imageset = datablocks[0].extract_imagesets()[0]
  detector = imageset.get_detector()

  h = detector.hierarchy()
  href = reference_detector.hierarchy()

  assert len(h) == len(href)

  assert abs(h.get_distance() - href.get_distance()) < params.max_delta_distance

  for panel, panel_ref in zip(h.children(), href.children()):
    panel.set_local_frame(
      panel_ref.get_fast_axis(),
      panel_ref.get_slow_axis(),
      panel_ref.get_origin()
    )

  print 'Writing metrology-corrected datablock to %s' %params.output.datablock
  from dxtbx.serialize import dump
  dump.datablock(datablocks, params.output.datablock)

  return

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
