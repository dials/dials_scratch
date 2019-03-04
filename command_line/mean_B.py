from __future__ import print_function
import iotbx.cif

def get_B(cif_file):
  cif_text = open(cif_file, 'r').read()

  assert 'data_' in cif_text

  for record in cif_text.split('\n'):
    if record.startswith('data_'):
      name = record.strip().replace('data_', '')
      break

  # strip out Uij == 0 for H atoms, then rebuild cif_text
  non_zero = []
  for record in cif_text.split('\n'):
    if record.startswith('H') and '0 0 0 0 0 0' in record:
      continue
    non_zero.append(record)
  cif_text = '\n'.join(non_zero)

  if name is '':
    cif_text = cif_text.replace('data_', 'data_35dnba')
    name = '35dnba'

  # read structure, remove hydrogen, make isotropic for analysis
  cs = iotbx.cif.reader(input_string=cif_text).build_crystal_structures()[name]
  csh = cs.element_selection('H')

  hs = cs.select(csh)
  hs.convert_to_isotropic()

  cs = cs.select(~csh)
  cs.convert_to_isotropic()

  return cs.b_iso_min_max_mean()[2], hs.b_iso_min_max_mean()[2]

if __name__ == '__main__':
  import sys

  bs = { }
  for filename in sys.argv[1:]:
    bs[filename] = get_B(filename)

  for filename in sorted(bs):
    print(filename, '%.4f %.4f' % bs[filename])
