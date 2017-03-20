from __future__ import division

def estimate_model_parameters(experiments, reflections):
  '''
  Estimate the model parameters

  '''
  from scitbx import simplex
  from copy import deepcopy
  from dials.array_family import flex
  from random import uniform

  # The number of parameters
  num_parameters = 6

  # Setup the starting simplex
  starting_simplex = []
  for i in range(num_parameters+1):
    starting_simplex.append(
      flex.log(flex.double([uniform(0.0001,0.01) for j in range(num_parameters)]))
    )

  class Evaluator(object):

    def __init__(self):
      from dials_scratch.jmp.profile_modelling import MLTarget3D
      self.func = MLTarget3D(experiments[0], reflections)
      self.count = 1

    def target(self, log_parameters):
      from dials.array_family import flex
      parameters = flex.exp(log_parameters)

      logL = self.func.log_likelihood(parameters)
 
      self.count += 1

      print self.count, list(parameters), logL

      # Return negative log likelihood 
      return -logL

  # Setup the simplex optimizer
  optimizer = simplex.simplex_opt(
    num_parameters,
    matrix    = starting_simplex,
    evaluator = Evaluator(),
    tolerance = 1e-7)

  # Get the solution
  parameters = flex.exp(optimizer.get_solution())

  # Return the current model
  return parameters


def generate_profile_model(experiments, reflections):
  '''
  Generate the reflection profile model

  '''

  def select_used_in_refinement(reflections):
    '''
    Select reflections to use

    '''
    selection = reflections.get_flags(reflections.flags.used_in_refinement)
    print 'Selecting %d/%d strong reflections' % (
      len(reflections), selection.count(True))
    print ""
    return reflections.select(selection)

  def sort_by_intensity(reflections):
    '''
    Sort the reflections by intensity

    '''
    reflections.sort("intensity.sum.value", reverse=True)
    return reflections

  def select_subset(reflections, num):
    '''
    Select most intense reflections

    '''
    print 'Selecting %d/%d reflections' % (num, len(reflections))
    return reflections[0:num]

  # def display(experiments, reflections, model, num):
  #   '''
  #   Display some shoeboxes

  #   '''
  #   from dials_scratch.jmp.viewer import show_image_stack_multi_view
  #   from random import sample
  #   from dials.array_family import flex

  #   # Sample from reflections
  #   reflections = reflections.select(flex.size_t(sample(range(len(reflections)), num)))

  #   # Simulate the reflection profiles from the current model
  #   simulated = simulate(experiments, reflections, model)
   
  #   # Display stuff
  #   for model_sbox, data_sbox in zip(simulated, reflections['shoebox']):
  #     model = model_sbox.data
  #     data = data_sbox.data
  #     show_image_stack_multi_view(model.as_numpy_array(), vmax=flex.max(model))
  #     show_image_stack_multi_view(data.as_numpy_array(), vmax=flex.max(data))


  # Select the strong reflections
  reflections = select_used_in_refinement(reflections)
  reflections = sort_by_intensity(reflections)

  # Select a subset of strong reflections
  subset = select_subset(reflections, 100)
  
  # Display a few
  # if False:
  #   display(experiments, subset, initial_model, 5)

  # Estimate the model parameters
  final_model = estimate_model_parameters(
    experiments, 
    subset)

  # Display a few
  # if True:
  #   display(experiments, subset, final_model, 5)

  # Return the final model
  return final_model


if __name__ == '__main__':

  def read_experiments(filename):
    '''
    Read the experiments file

    '''
    from dxtbx.model.experiment_list import ExperimentListFactory
    print "Reading experiments from %s" % filename
    print ""
    return ExperimentListFactory.from_json_file(filename)

  def read_reflections(filename):
    '''
    Read the reflections file

    '''
    from dials.array_family import flex
    print "Reading reflections from %s" % filename
    print ""
    return flex.reflection_table.from_pickle(filename)

  # Hard code the filenames
  experiments_filename = '/home/upc86896/Data/bag_training/processed_profile/profile_model/experiments.json'
  reflections_filename = '/home/upc86896/Data/bag_training/processed_profile/profile_model/reflections.pickle'

  # Read the experiments and reflections
  experiments = read_experiments(experiments_filename)
  reflections = read_reflections(reflections_filename)

  # Generate the profile model
  final_model = generate_profile_model(experiments, reflections)
    
  print "Generated final model:"
  print ""
  print list(final_model)
  print ""

