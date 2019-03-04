from __future__ import division
from __future__ import print_function

import os
import iotbx.phil
from libtbx.phil import command_line
from cctbx import sgtbx
from cctbx.array_family import flex
from iotbx.reflection_file_reader import any_reflection_file


phil_scope = iotbx.phil.parse('''\
change_of_basis_op = a,b,c
  .type = str

output {
  suffix = "_reindexed"
    .type = str
}
''', process_includes=True)

def run(args):
  import libtbx
  from libtbx import easy_pickle

  cmd_line = command_line.argument_interpreter(master_params=phil_scope)
  working_phil, files = cmd_line.process_and_fetch(
    args=args, custom_processor="collect_remaining")
  working_phil.show()
  params = working_phil.extract()

  datasets = []

  for file_name in files:

    reader = any_reflection_file(file_name)
    assert reader.file_type() == 'ccp4_mtz'
    mtz_object = reader.file_content()

    cb_op = sgtbx.change_of_basis_op(params.change_of_basis_op)
    basename = os.path.basename(file_name)
    out_name = os.path.splitext(basename)[0] + params.output.suffix + '.mtz'
    print("reindexing %s (%s)" %(file_name, cb_op.as_xyz()))
    mtz_object.change_basis_in_place(cb_op)
    mtz_object.write(out_name)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
