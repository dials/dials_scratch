
class FractionOfObservedIntensity(object):
  '''Calculate the fraction of observed intensity for different sigma_m.'''

  def __init__(self, crystal, beam, detector, goniometer, scan, reflections):
    '''Initialise the algorithm. Calculate the list of tau and zetas.

    Params:
        reflections The list of reflections
        experiment The experiment object

    '''
    from dials.array_family import flex
    from math import sqrt

    # Get the oscillation width
    dphi2 = scan.get_oscillation(deg=False)[1] / 2.0

    # Calculate a list of angles and zeta's
    tau, zeta = self._calculate_tau_and_zeta(crystal, beam, detector,
                                             goniometer, scan, reflections)

    # Calculate zeta * (tau +- dphi / 2) / sqrt(2)
    self.e1 = (tau + dphi2) * flex.abs(zeta)
    self.e2 = (tau - dphi2) * flex.abs(zeta)

  def _calculate_tau_and_zeta(self, crystal, beam, detector, goniometer, scan, reflections):
    '''Calculate the list of tau and zeta needed for the calculation.

    Params:
        reflections The list of reflections
        experiment The experiment object.

    Returns:
        (list of tau, list of zeta)

    '''
    from scitbx.array_family import flex

    # Calculate the list of frames and z coords
    bbox = reflections['bbox']
    phi = reflections['xyzcal.mm'].parts()[2]

    # Calculate the zeta list
    zeta = reflections['zeta']

    # Calculate the list of tau values
    tau = []
    zeta2 = []
    scan = scan
    for b, p, z in zip(bbox, phi, zeta):
      for f in range(b[4], b[5]):
        phi0 = scan.get_angle_from_array_index(int(f), deg=False)
        phi1 = scan.get_angle_from_array_index(int(f)+1, deg=False)
        tau.append((phi1 + phi0) / 2.0 - p)
        zeta2.append(z)

    # Return the list of tau and zeta
    return flex.double(tau), flex.double(zeta2)

  def __call__(self, sigma_m):
    '''Calculate the fraction of observed intensity for each observation.

    Params:
        sigma_m The mosaicity

    Returns:
        A list of log intensity fractions

    '''
    from math import sqrt
    from scitbx.array_family import flex
    import scitbx.math

    # Tiny value
    TINY = 1e-10
    assert(sigma_m > TINY)

    # Calculate the two components to the fraction
    a = scitbx.math.erf(self.e1 / (sqrt(2) * sigma_m))
    b = scitbx.math.erf(self.e2 / (sqrt(2) * sigma_m))

    # Calculate the fraction of observed reflection intensity
    R = (a - b) / 2.0
    

    # Set any points <= 0 to 1e-10 (otherwise will get a floating
    # point error in log calculation below).
    assert(R.all_ge(0))
    mask = R < TINY
    #print "NUM: ", mask.count(True)
    assert(mask.count(True) < len(mask))
    R.set_selected(mask, TINY)

    # Return the logarithm of r
    return flex.log(R)


 
