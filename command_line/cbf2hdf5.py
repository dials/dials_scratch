from __future__ import division

if __name__ == '__main__':

  import sys
  from dxtbx.model.experiment_list import ExperimentListFactory
  
  experiments = ExperimentListFactory.from_json_file(sys.argv[1])
  assert len(experiments) == 1
  assert len(experiments[0].imageset) == 1

  from dials.util.nexus import get_entry  
  from dials.util.nexus import nx_mx
  entry = get_entry("data.h5", "w")
  experiment_names = nx_mx.dump(entry, experiments)

  from dials.array_family import flex
  data = experiments[0].imageset.get_raw_data(0)[0]
  height, width = data.all()
  data.reshape(flex.grid(1, height, width))

  nx_data = entry[experiment_names[0]].create_group("data")
  nx_data.attrs['NX_class'] = 'NXdata'
  nx_data['data'] = data.as_numpy_array()
  
  entry[experiment_names[0]]['names'] = ['test_cbf'] + ["" for i in range(199)]
