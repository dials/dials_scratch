import iotbx.phil

phil_scope = iotbx.phil.parse('''

  parameterisation {
    scale_term = True
      .type = bool
      .help = "Option to turn off decay correction (only for KB scaling)"
    rotation_interval = 15.0
      .type = float
      .help = "User specified rotation (phi) interval in degrees for phi binning
              for the scale term"
    decay_term = True
      .type = bool
      .help = "Option to turn off decay correction"
    B_factor_interval = 20.0
      .type = float
      .help = "User specified rotation (phi) interval in degrees for phi binning
              for the decay term"
    absorption_term = True
      .type = bool
      .help = "Option to turn off absorption correction"
    lmax = 4
      .type = int
      .help = "Number of spherical harmonics to include for absorption correction,
              recommended to be no more than 6."
  }
  reflection_selection {
    E2min = 0.8
      .type = float
      .help = "Minimum normalised E^2 value to select reflections for scaling"
    E2max = 5.0
      .type = float
      .help = "Maximum normalised E^2 value to select reflections for scaling"
    Isigma_min = -5.0
      .type = float
      .help = "Option to use a I/sigma subset of reflections to determine scale factors"
    d_min = 0.0
      .type = float
      .help = "Option to use a d-value subset of reflections to determine scale factors"
  }
  scaling_options {
    force_space_group = None
      .type = str
      .help = "Option to specify space group for scaling"
    concurrent_scaling = True
      .type = bool
      .help = "Option to allow absorption correction after decay/scale, 
              if concurrent_scaling is set to False"
    optimise_error_model = True
      .type = bool
      .help = "Option to allow optimisation of weights for scaling. Performs
               and additional scale factor minimisation after adjusting weights."
    error_model_params = None
      .type = floats(size=2)
      .help = "Ability to force an error model adjustment, using the model 
              in aimless - factors are called SDFac, SDadd in aimless."
    reject_outliers = True
      .type = bool
      .help = "Option to turn on outlier rejection"
    verbosity = 1
      .type = int(value_min=0)
      .help = "The verbosity level"
    integration_method = 'prf'
      .type = str
      .help = "Option to choose from profile fitted intensities (prf)
              or summation integrated intensities (sum)"
    minimisation_parameterisation = 'standard'
      .type = str
      .help = "Choice of 'standard' (multiplicative) or 'log' g-value 
               minimisation parameterisation"
    target = None
      .type = str
      .help = "Choice to specify a target dataset for scaling"
  }

  ''')