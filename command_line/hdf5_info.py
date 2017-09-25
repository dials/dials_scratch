import h5py
import numpy

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
  main(sys.argv[1])
