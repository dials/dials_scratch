import logging
from dials.array_family import flex

logger = logging.getLogger('dials')

class ActiveParameterFactory(object):
  '''class to create and return appropriate apm factory.'''

  @classmethod
  def create(cls, scaler):
    '''create the factory based on phil options'''
    if scaler.params.scaling_options.concurrent_scaling:
      return ConcurrentParameterFactory(scaler)
    else:
      return ConsecutiveParameterFactory(scaler)


class ConcurrentParameterFactory(object):
  '''factory to set up apms for concurrent scaling of all active parameters'''

  def __init__(self, scaler):
    self.scaler = scaler
    self.apm = None
    self.param_lists = []
    self.create_active_list()
    self.n_cycles = 1

  def create_active_list(self):
    '''return a list indicating the names of active parameters'''
    #replace these methods with some kind of entry point mechanism to allow
    #scalers derived from new base scalers in future?
    from dials_scratch.jbe.scaling_code import Scaler
    if isinstance(self.scaler, Scaler.SingleScalerBase):
      param_name = []
      for param in self.scaler.corrections:
        param_name.append(str(param))
      if not param_name:
        assert 0, 'no parameters have been chosen for scaling, aborting process'
      self.param_lists = param_name
      self.apm = active_parameter_manager(self.scaler, self.param_lists)

    elif isinstance(self.scaler, Scaler.MultiScalerBase):
      if self.scaler.id_ == 'target':
        scalers = self.scaler.unscaled_scalers
      elif self.scaler.id_ == 'multi':
        scalers = self.scaler.single_scalers
      else:
        assert 0, 'unrecognised scaler id_ for scaler derived from MultiScalerBase'

      for scaler in scalers:
        param_name = []
        for param in scaler.corrections:
          param_name.append(str(param))
        if not param_name:
          assert 0, 'no parameters have been chosen for scaling, aborting process'
        self.param_lists.append(param_name)
      if self.scaler.id_ == 'target':
        self.apm = target_active_parameter_manager(self.scaler, self.param_lists)
      elif self.scaler.id_ == 'multi':
        self.apm = multi_active_parameter_manager(self.scaler, self.param_lists)

    else:
      assert 0, '''scaler not derived from Scaler.SingleScalerBase or
        Scaler.MultiScalerBase. An additional option must be defined in ParameterHandler'''

  def make_next_apm(self):
    'method to call to return the apm'
    return self.apm

class ConsecutiveParameterFactory(object):
  '''factory to set up apms for consecutive scaling of all active parameters'''

  def __init__(self, scaler):
    self.scaler = scaler
    self.param_lists = []
    self.n_cycles = None
    self.create_consecutive_list()
    print(self.param_lists)
    print(self.n_cycles)

  def create_consecutive_list(self):
    '''return a list indicating the names of active parameters'''
    #replace these methods with some kind of entry point mechanism to allow
    #scalers derived from new base scalers in future?
    from dials_scratch.jbe.scaling_code import Scaler
    if isinstance(self.scaler, Scaler.SingleScalerBase):
      corrections = self.scaler.experiments.scaling_model.configdict['corrections']
      for cycle in self.scaler.consecutive_scaling:
        corrlist = []
        for corr in cycle:
          if corr in corrections:
            corrlist.append(corr)
        self.param_lists.append(corrlist)
      self.n_cycles = sum([1 for i in self.param_lists if i])

    elif isinstance(self.scaler, Scaler.MultiScalerBase):
      if self.scaler.id_ == 'target':
        scalers = self.scaler.unscaled_scalers
      elif self.scaler.id_ == 'multi':
        scalers = self.scaler.single_scalers
      else:
        assert 0, 'unrecognised scaler id_ for scaler derived from MultiScalerBase'

      for scaler in scalers:
        ind_param_list = []
        corrections = scaler.experiments.scaling_model.configdict['corrections']
        for cycle in scaler.consecutive_scaling:
          corrlist = []
          for corr in cycle:
            if corr in corrections:
              corrlist.append(corr)
          ind_param_list.append(corrlist)
        self.param_lists.append(ind_param_list)
      #now need to calculate the max number of cycles needed across all scalers
      is_cycle_active = []
      for p_list in self.param_lists:
        for i, cycle in enumerate(p_list):
          if cycle:
            is_cycle_active.append(i)
      self.n_cycles = len(set(is_cycle_active))
      #now make sure all lists are same length
      max_len = max([len(i) for i in self.param_lists])
      for p_list in self.param_lists:
        for _ in range(max_len - len(p_list)):
          p_list.append([])
    else:
      assert 0, '''scaler not derived from Scaler.SingleScalerBase or
        Scaler.MultiScalerBase. An additional option must be defined in ParameterHandler'''

  def make_next_apm(self):
    '''generate a valid apm for minimisation (contains some active parameters,
    but not necessarily for all datasets)'''
    from dials_scratch.jbe.scaling_code import Scaler
    if isinstance(self.scaler, Scaler.SingleScalerBase):
      apm = active_parameter_manager(self.scaler, self.param_lists[0])
      self.param_lists = self.param_lists[1:]
    elif isinstance(self.scaler, Scaler.MultiScalerBase):
      params = []
      for i in range(0, len(self.param_lists)):
        params.append(self.param_lists[i][0])
        self.param_lists[i] = self.param_lists[i][1:] #remove first element
      if self.scaler.id_ == 'multi':
        apm = multi_active_parameter_manager(self.scaler, params)
      elif self.scaler.id_ == 'target':
        apm = target_active_parameter_manager(self.scaler, params)
      else:
        assert 0, 'unrecognised scaler id_ for scaler derived from MultiScalerBase'
    else:
      assert 0, '''scaler not derived from SingleScalerBase or MultiScalerBase.
        An additional option must be defined in ParameterHandler'''
    if not apm.active_parameterisation: #no active parameters, so iterate
      apm = self.make_next_apm()
    return apm

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
    for p_type, scalefactor in scaler.experiments.scaling_model.components.iteritems():
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
      ).format(''.join(str(i)+' ' for i in self.active_parameterisation)))


class multi_active_parameter_manager(object):
  ''' object to manage the current active parameters during minimisation
  for multiple datasets that are simultaneously being scaled.'''
  def __init__(self, multiscaler, param_list):
    self.apm_list = []
    for i, scaler in enumerate(multiscaler.single_scalers):
      self.apm_list.append(active_parameter_manager(scaler, param_list[i]))
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
      ).format(''.join(str(i)+' ' for i in self.active_parameterisation)))

class target_active_parameter_manager(object):
  ''' object to manage the current active parameters during minimisation
  for multiple datasets that are simultaneously being scaled.'''
  def __init__(self, targetscaler, param_list):
    self.apm_list = []
    for i, scaler in enumerate(targetscaler.unscaled_scalers):
      self.apm_list.append(active_parameter_manager(scaler, param_list[i]))
    self.active_parameterisation = []
    for apm in self.apm_list:
      self.active_parameterisation.extend(apm.active_parameterisation)
    self.x = flex.double([])
    self.n_params_in_each_apm = []
    self.n_cumul_params_list = [0]
    for apm in self.apm_list:
      self.x.extend(apm.x)
      self.n_params_in_each_apm.append(len(apm.x))
      self.n_cumul_params_list.append(len(self.x))
    self.n_active_params = len(self.x)
    logger.info(('Set up joint parameter handler for following corrections: {0}\n'
      ).format(''.join(str(i)+' ' for i in self.active_parameterisation)))
