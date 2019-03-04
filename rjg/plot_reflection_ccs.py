from __future__ import division
from __future__ import print_function
import json
import os

def run(args):
  from matplotlib import pyplot
  fig = pyplot.figure()

  input_files = [arg for arg in args if os.path.isfile(arg)]
  labels = [arg for arg in args if not os.path.isfile(arg)]
  assert len(input_files) == len(labels)
  for f, l in zip(input_files, labels):
    print(l, f)
    d = json.load(open(f, 'rb'))
    reflection_ccs = d['reference']['reflection_cc_vs_resolution']
    data = reflection_ccs['data'][0]
    pyplot.plot(data['x'], data['y'], label=l)

  layout = reflection_ccs['layout']
  pyplot.xticks(layout['xaxis']['tickvals'], layout['xaxis']['ticktext'])
  pyplot.xlabel(layout['xaxis']['title'])
  pyplot.ylabel(layout['yaxis']['title'])
  pyplot.legend(loc='upper right')
  #fig.savefig('reflection_cc_vs_resolution.png')

  pyplot.show()

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])

