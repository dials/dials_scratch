from __future__ import division

from table import Table

class ReflectionTable(Table):

  def spatial_index(self, column):
    from dials.algorithms.spatial_indexing import make_spatial_index
    return make_spatial_index(self._columns[column])

def ReflectionTableMTZ(mtz_file):
  '''Create a ReflectionTable from an MTZ file.'''

  from iotbx import mtz

  table = ReflectionTable()

  m = mtz.object(mtz_file)
  mi = m.extract_miller_indices()

  table = ReflectionTable()
  table["miller_index"] = mi

  for c in m.columns():
    if c.type() == 'H':
      continue
    if c.type() in 'BYI':
      table[c.label()] = m.extract_integers(c.label()).data
    else:
      table[c.label()] = c.extract_values()

  crystals = m.crystals()

  table.metadata['unit_cell'] = crystals[0].unit_cell()
  table.metadata['max_min_resolution'] = m.max_min_resolution()

  # FIXME add more metadata?

  return table

if __name__ == '__main__':
  from cctbx.array_family import flex

  # Create a reflection table with miller indices and intensities
  table = ReflectionTable([
    ('miller_index',
     flex.miller_index([(0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 0, 2), (0, 0, 1)])
    ),
    ('intensity',
     flex.double([1, 2, 3, 4, 5])
    ),
    ('centroid_px',
     flex.vec3_double([(0, 0, 0), (-1, -1, 0), (-1, 1, 0), (1, 1, 0), (1, -1, 0)])
    )
  ])

  print '\nGet the list of intensities for each miller index'
  indexer = table.index_map('miller_index')
  for key, indices in indexer.iteritems():
    print "Miller Index: {0}, Intensities: {1}".format(
      key, [table[i].intensity for i in indices])

  print '\nGet a list of reflections recorded in a volume'
  indexer = table.spatial_index('centroid_px')

  vol = (-1, 1, -1, 1, -1, 1)
  print "Volume: ", vol
  for i in indexer.query_range(vol):
    print "  Miller Index: {0}, Centroid: {1}".format(
      table[i].miller_index, table[i].centroid_px)

  vol = (0, 1, -1, 1, -1, 1)
  print "Volume: ", vol
  for i in indexer.query_range(vol):
    print "  Miller Index: {0}, Centroid: {1}".format(
      table[i].miller_index, table[i].centroid_px)

  print '\nSet some global metadata'
  table.metadata['some_data'] = "Some global meta data"
  print table.metadata