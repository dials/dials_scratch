import libtbx.load_env

def run(args):
  assert len(args) == 1
  json_file = args[0]
  import json

  with open(json_file, 'rb') as f:
    results = json.load(f)

  n_indexed = []
  n_spots = []

  for r in results:
    n_spots.append(r['n_spots_total'])
    if 'n_indexed' in r:
      n_indexed.append(r['n_indexed'])
    else:
      n_indexed.append(0)

  from matplotlib import pyplot
  pyplot.scatter(n_spots, n_indexed)
  pyplot.plot([0, max(n_spots)], [0, max(n_spots)], c='red')
  pyplot.show()

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
