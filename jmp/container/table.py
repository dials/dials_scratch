from __future__ import division


class Table(object):

  def __init__(self, columns=None, metadata=None):
    from collections import OrderedDict
    if columns == None:
      self._columns = OrderedDict()
    else:
      self._columns = OrderedDict(columns)
    if not self.is_consistent():
      raise RuntimeError('inconsistent column sizes')
    if metadata == None:
      self._metadata = dict()
    else:
      self._metadata = metadata

  def __getitem__(self, index):
    if isinstance(index, int):
      d = dict((k, c[index]) for k, c in self._columns.iteritems())
      return self.make_record(**d)
    elif isinstance(index, slice):
      return Table(
        [(k, c[index]) for k, c in self._columns.iteritems()],
        self._metadata)
    elif isinstance(index, str):
      return type(self._columns[index])(self._columns[index])
    else:
      raise TypeError('expected int, str or slice, got %s' % type(index))

  def __setitem__(self, index, item):
    if isinstance(index, int):
      if isinstance(item, tuple):
        for i, key in enumerate(self._columns):
          self._columns[key][index] = item[i]
      elif isinstance(item, self.record_type()):
        for key in self._columns:
          self._columns[key][index] = getattr(item, key)
      else:
        raise TypeError('expected tuple or record, got %s' % type(item))
    elif isinstance(index, str):
      if self.ncols() != 0 and len(item) != self.nrows():
        raise RuntimeError('inconsistent column sizes')
      self._columns[index] = type(item)(item)
    else:
      raise TypeError('expected int, str or slice, got %s' % type(index))

  def __len__(self):
    return self.nrows()

  def __repr__(self):
    return '%s: %s' % (self.__class__, str(self._columns))

  def __iter__(self):
    for row in self.iterrows():
      yield row

  def __contains__(self, column):
    return column in self._columns

  def iterrows(self):
    for i in range(len(self)):
      yield self.__getitem__(i)

  def iterheaders(self):
    for header in self._columns:
      yield header

  def itercolumns(self):
    for header, column in self._columns.iteritems():
      yield (header, column)

  def append(self, item):
    if isinstance(item, tuple):
      for i, key in enumerate(self._columns):
        self._columns[key].append(item[i])
    elif isinstance(item, self.record_type()):
      for key in self._columns:
        self._columns[key].append(getattr(item, key))
    else:
      raise TypeError('expected tuple or record, got %s' % type(item))

  def remove(self, index):
    for key in self._columns:
      self._columns.pop(index)

  def remove_column(self, name):
    self._columns.remove(name)

  def extend(self, table):
    for key in self._columns:
      self._columns[key].extend(table._columns[key])

  def make_record(self, **kwargs):
    return self.record_type()(**kwargs)

  def record_type(self):
    from collections import namedtuple
    return namedtuple("Record", self._columns.iterkeys())

  def extract(self, items, ltype=None, func=None):
    columns = [self._columns[name] for name in items]
    if func is not None:
      result = func(columns)
      assert(self.is_consistent())
      return result
    if ltype is None:
      ltype = list
    return ltype(zip(*columns))

  def is_consistent(self):
    length = [len(c) for c in self._columns.itervalues()]
    if len(length) == 0:
      return True
    return length.count(length[0]) == len(length)

  def nrows(self):
    if len(self._columns) == 0:
      return 0
    else:
      return len(self._columns[self._columns.keys()[0]])

  def ncols(self):
    return len(self._columns)

  def index_map(self, column):
    from collections import defaultdict
    index = defaultdict(list)
    for i, item in enumerate(self._columns[column]):
      index[item].append(i)
    return index

  @property
  def metadata(self):
    return self._metadata


if __name__ == '__main__':
  from scitbx.array_family import flex

  # Create a table with a list of columns
  table = Table([
    ("column_1", flex.int()),
    ("column_2", flex.int())])

  # Append rows to the columns either as tuples or named tuples
  table.append((1, 10))
  table.append((2, 20))
  table.append(table.make_record(column_1=3, column_2=30))

  # Extend one table with another
  table.extend(table)

  # Get the number of rows and columns
  print "Access some data"
  print "N Rows: ", len(table)
  print "N Rows: ", table.nrows()
  print "N Cols: ", table.ncols()

  # Access a row of the table
  print "Row 2: ", table[2]

  # Access a columns of the table
  print "Column 2: ", table['column_2']

  # Add additional columns (must be the correct length)
  table['column_3'] = flex.double(6)
  table['column_4'] = flex.double(6)

  # Set the values of particular rows by tuple or named tuple
  table[3] = (100, 200, 300.0, 400.0)
  table[2] = table.make_record(column_1=1, column_2=2, column_3=3, column_4=4)

  # Iterate through stuff
  print '\nIterate through the rows of the table'
  for row in table:
      print row

  print '\nIterate through the column headings of the table'
  for header in table.iterheaders():
      print header

  print '\nIterate through the columns of the table'
  for header, data in table.itercolumns():
      print header, data

  print '\nSlice the table (copies data)'
  table2 = table[2:4]
  print "N Rows: ", table2.nrows()

  print '\nExtract columns into arbitrary lists'
  table['column_5'] = flex.double(6)
  l = table.extract(('column_3', 'column_4', 'column_5'), ltype=flex.vec3_double)
  print "Extracted: ", type(l), len(l)
  l = table.extract(('column_3', 'column_4', 'column_5'), func=lambda x: list(zip(*x)))
  print "Extracted: ", type(l), len(l)

  print '\nTest column is available'
  print "Is Column 4 available: ", 'column_4' in table
  print "Is Column 8 available: ", 'column_8' in table

  print '\nTry bad column size'
  try:
    table['columns_6'] = flex.int(10)
  except Exception, e:
    print "Exception raised: ", e

  print '\nGet the indices for unique column_1 values'
  indexer = table.index_map('column_1')
  for key, indices in indexer.iteritems():
    print key, [table[i] for i in indices]

  print '\nSet some metadata'
  table.metadata['some_data'] = "hello world"
  print table.metadata

  # Create table with nothing in
  table = Table()

  # Set items
  table['a'] = flex.int(1000000)
  table['b'] = flex.int(1000000)
  table['c'] = flex.int(1000000)