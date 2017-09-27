import h5py
import numpy

def extract(filename, dataset_list):

  f = h5py.File(filename, 'r')

  data = { }

  for d in dataset_list:
    data[d] = f[d]

  datasets = [data[d] for d in dataset_list]

  print ' '.join(dataset_list)

  for j, record in enumerate(zip(*datasets)):
    print '%d %s' % (j, ' '.join(map(str, record)))

def main(filename):

  f = h5py.File(filename, 'r')

  datasets = []

  def visitor(name, obj):
    if isinstance(obj, h5py.Dataset):
      datasets.append(name)
      print 'dataset:', name
      for k in obj.attrs:
        print k, obj.attrs[k]
    elif isinstance(obj, h5py.Group):
      print 'group:  ', name

    print

  f.visititems(visitor)

  for d in datasets:
    print d
    print f[d]
    print

if __name__ == '__main__':
  import sys
  if len(sys.argv) == 2:
    main(sys.argv[1])
  else:
    extract(sys.argv[1], sys.argv[2:])
