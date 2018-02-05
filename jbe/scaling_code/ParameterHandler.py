import logging
from dials.array_family import flex

logger = logging.getLogger('dials')


class ParameterlistFactory(object):
  '''
  Factory to create parameter lists to pass to a minimiser.
  The methods should return a nested list
  i.e [['g_scale','g_decay','g_absorption']]
  or [['g_scale','g_decay'],['g_absorption']]
  '''
  @classmethod
  def full_active_list(cls, experiments):
    '''create a list with all params to include'''
    
    param_name = []
    for param in experiments.scaling_model.configdict['corrections']:
      param_name.append('g_'+str(param))
    if not param_name:
      assert 0, 'no parameters have been chosen for scaling, aborting process'
    return [param_name]

  @classmethod
  def consecutive_list(cls, experiments):
    '''create a nested list with all params to include'''
    param_name = []
    corrections = experiments.scaling_model.configdict['corrections']
    if experiments.scaling_model.id_ == 'aimless':
      if 'scale' in corrections and 'decay' in corrections:
        param_name.append(['g_scale', 'g_decay'])
      elif 'scale' in corrections:
        param_name.append(['g_scale'])
      elif 'decay' in corrections:
        param_name.append(['g_decay'])
      if 'absorption' in experiments.scaling_model.configdict['corrections']:
        param_name.append(['g_absorption'])
    else:
      for param in corrections:
        param_name.append(['g_'+str(param)])
    return param_name


class active_parameter_manager(object):
  ''' object to manage the current active parameters during minimisation.
  Separated out to provide a consistent interface between the scaler and
  minimiser. Takes in a scaler, needed to access SF objects through
  g_parameterisation, and a param_name list indicating the active parameters.'''
  def __init__(self, scaler, param_name):
    self.constant_g_values = None
    self.x = flex.double([])
    self.active_parameterisation = []
    self.n_params_list = [] #no of params in each SF
    self.n_cumul_params_list = [0]
    self.active_derivatives = None
    for p_type, scalefactor in scaler.g_parameterisation.iteritems():
      if p_type in param_name:
        self.x.extend(scalefactor.parameters)
        self.active_parameterisation.append(p_type)
        self.n_params_list.append(scalefactor.n_params)
        self.n_cumul_params_list.append(len(self.x))
      else:
        if self.constant_g_values is None:
          self.constant_g_values = scalefactor.inverse_scales
        else:
          self.constant_g_values *= scalefactor.inverse_scales
    self.n_active_params = len(self.x)
    logger.info(('Set up parameter handler for following corrections: {0}\n'
      ).format(''.join(i.lstrip('g_')+' ' for i in self.active_parameterisation)))


class multi_active_parameter_manager(object):
  ''' object to manage the current active parameters during minimisation
  for multiple datasets that are simultaneously being scaled.'''
  def __init__(self, multiscaler, param_name):
    self.apm_list = []
    for scaler in multiscaler.single_scalers:
      self.apm_list.append(active_parameter_manager(scaler, param_name))
    self.active_parameterisation = []
    for apm in self.apm_list:
      self.active_parameterisation.extend(apm.active_parameterisation)
    self.x = flex.double([])
    self.n_params_in_each_apm = []
    self.n_cumul_params_list = [0]
    self.active_derivatives = None
    for apm in self.apm_list:
      self.x.extend(apm.x)
      self.n_params_in_each_apm.append(len(apm.x))
      self.n_cumul_params_list.append(len(self.x))
    self.n_active_params = len(self.x)
    logger.info(('Set up joint parameter handler for following corrections: {0}\n'
      ).format(''.join(i.lstrip('g_')+' ' for i in self.active_parameterisation)))
