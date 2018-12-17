from __future__ import division
from libtbx import phil
from collections import Counter, defaultdict
from dials.array_family import flex
from math import floor, sqrt, ceil
import cPickle as pickle

phil_scope = phil.parse('''

  tolerance = 1
    .type = float
    .help = "Tolerance to test for convergence (difference in log(L))"
  
  max_iter = 10000
    .type = int
    .help = "Maximum number of iterations"

  grid_size = 50,20
    .type = ints(size=2)
    .help = "The grid of partiality priors"

''')

def compute_Ep(rho0, var_rho0, mu, I, V):
  var_rho1 = 1.0 / (1/var_rho0 + mu**2/V)
  rho1 = (rho0/var_rho0 + I*mu/V)*var_rho1
  return rho1, var_rho1+rho1**2


class Scaler(object):

  def __init__(self, 
               tolerance=1e-3, 
               max_iter=10000,
               grid_size=(20,20)):

    # Set the parameters
    self.tolerance = tolerance
    self.max_iter = max_iter
    self.grid_size = grid_size

    # Initialise the reflections
    self.reflections = None

    # Initialise the merged intensities
    self.h_unique = None
    self.I_mean = None
    self.V_mean = None

    # Initialise the scales
    self.image_scale = None
    self.epsilon_mean = None
    self.epsilon_variance = None

  def scale(self, experiments, reflections):
  
    # Preprocess the data
    self.reflections = self._preprocess(experiments, reflections)

    # Solve the equations
    self._solve(
      self.reflections['intensity.sum.value'],
      self.reflections['intensity.sum.variance'],
      self.reflections['image'],
      self.reflections['reflection_id'],
      self.reflections['epsilon_id'])

  def _preprocess(self, experiments, reflections):

    # Filter reflections by:
    # - integrated
    # - non-zero variance
    # - valid id
    selection1 = reflections.get_flags(reflections.flags.integrated_sum)
    selection2 = reflections['intensity.sum.variance'] > 0
    selection3 = reflections['id'] >= 0
    selection = selection1 & selection2 & selection3
    reflections = reflections.select(selection)

    # Convert id to image 
    reflections['image'] = flex.size_t(list(reflections['id']))

    def compute_epsilon(reflections):
      s2 = reflections['s2']
      s1 = reflections['s1']
      epsilon = (s2.norms() - s1.norms())
      reflections['epsilon'] = epsilon
  
    # Compute the resolution
    reflections.compute_d(experiments)
  
    # Compute distance to Ewald sphere
    compute_epsilon(reflections)

    # Compute the miller indices in the asu
    reflections.compute_miller_indices_in_asu(experiments)

    # Compute the unique miller indices
    h_unique = flex.miller_index(sorted(list(set(reflections['miller_index_asu']))))
    h_lookup = dict((h,i) for i,h in enumerate(h_unique)) 
    h_index = flex.size_t()
    for h in reflections['miller_index_asu']:
      h_index.append(h_lookup[h])

    # Save the reflection data
    reflections['reflection_id'] = h_index
    
    # Create a lookup for reflection index
    lookup_index = [flex.size_t() for i in range(len(h_unique))]
    for i, h in enumerate(reflections['reflection_id']):
      lookup_index[h].append(i)

    # Create a lookup for image number
    num_images = len(set(reflections['id']))
    assert num_images == flex.max(reflections['id'])+1
    lookup_image = [flex.size_t() for i in range(num_images)]
    for i, j in enumerate(reflections['id']):
      lookup_image[j].append(i)
  
    # Compute the epsilon bins
    M = self.grid_size[0]
    abs_epsilon = flex.abs(reflections['epsilon'])
    sorted_index = flex.size_t(sorted(range(len(abs_epsilon)), key=lambda x: abs_epsilon[x]))
    eps_index = flex.size_t(len(reflections))
    lookup_epsilon = []
    step = int(ceil(len(sorted_index) / M))
    for j in range(M):
      m1 = j * step
      m2 = (j+1)*step
      selection = sorted_index[m1:m2]
      lookup_epsilon.append(selection)
      eps_index.set_selected(selection, j)

    # Set the epsilon index
    reflections['epsilon_id'] = eps_index

    # Save the arrays
    self.h_unique = h_unique
    self.lookup_index = lookup_index
    self.lookup_image = lookup_image
    self.lookup_epsilon = lookup_epsilon
  
    # Compute multiplicity and # obs on image
    multiplicity = flex.int(len(x) for x in lookup_index)
    obs_on_image = flex.int(len(x) for x in lookup_image)
    obs_on_epbin = flex.int(len(x) for x in lookup_epsilon)
    print len(reflections), flex.sum(obs_on_epbin)
    assert flex.sum(multiplicity) == len(reflections)
    assert flex.sum(obs_on_image) == len(reflections)
    assert flex.sum(obs_on_epbin) == len(reflections)

    # Print some info
    print "# Total reflections:  %d" % len(reflections)
    print "# Unique reflections: %d" % len(self.h_unique)
    print "# Images: %d" % len(self.lookup_image)
    print "Min multiplicity: %d" % flex.min(multiplicity)
    print "Max multiplicity: %d" % flex.max(multiplicity)
    print "Min observations on image: %d" % flex.min(obs_on_image)
    print "Max observations on image: %d" % flex.max(obs_on_image)
    print "Min observations in epsilon bin: %d" % flex.min(obs_on_epbin)
    print "Max observations in epsilon bin: %d" % flex.max(obs_on_epbin)

    # Return reflections
    return reflections


  def _solve(self, I, V, image, ref_index, eps_index):

    # Get number of classes
    num_images = len(self.lookup_image)
    num_epbins = len(self.lookup_epsilon)

    # Initialise the scale factors
    g_old = flex.double(num_images, 1.0)
    m_old = flex.double(num_epbins, 1.0)
    s_old = flex.double(num_epbins, 0.05**2)

    # Initialise the intensities
    c_old = flex.double()
    for selection in self.lookup_index:
      I_i = I.select(selection)
      c_old.append(flex.sum(I_i)/len(I_i))
    
    # Set the starting logL to something small
    logL_old = -1e50
    
    # Iterate until max iter has been reached, the change in log likelihood is
    # below the tolerance indicating convergence or a keyboard interrupt is
    # received.
    for iteration in range(self.max_iter):
      try:

        # Get the scale factors for each reflection
        c = c_old.select(ref_index)
        m = m_old.select(eps_index)
        s = s_old.select(eps_index)
        g = g_old.select(image)

        # Compute the expected value of rho and rho**2
        E_p = flex.double()
        E_p2 = flex.double()
        for i in range(len(c)):
          m1, m2 = compute_Ep(m[i], s[i], g[i]*c[i], I[i], V[i])
          E_p.append(m1)
          E_p2.append(m2)

        # Compute the new values for the epsilon scale factors
        m_new = flex.double()
        s_new = flex.double()
        for i, selection in enumerate(self.lookup_epsilon):
          E_p_i = E_p.select(selection)
          E_p2_i = E_p2.select(selection)
          n = len(E_p_i)
          mm = flex.sum(E_p_i) / n
          ss = flex.sum((E_p2_i-2*m_old[i]*E_p_i+m_old[i]**2)) / n
          m_new.append(mm)
          s_new.append(ss)
        m_new /= m_new[0]
        s_new = s_old

        # Compute the new values for the intensities
        c_new = flex.double()
        for selection in self.lookup_index:
          E_p_i = E_p.select(selection)
          E_p2_i = E_p2.select(selection)
          yy = I.select(selection)
          ss = V.select(selection)
          gg = g.select(selection)
          cc = flex.sum(yy*gg*E_p_i/ss)/flex.sum(gg**2*E_p2_i/ss)
          c_new.append(cc)

        # Compute the new values for the image scale factors 
        g_new = flex.double()
        for selection in self.lookup_image:
          E_p_i = E_p.select(selection)
          E_p2_i = E_p2.select(selection)
          yy = I.select(selection)
          ss = V.select(selection)
          cc = c.select(selection)
          gg = flex.sum(yy*cc*E_p_i/ss)/flex.sum(cc**2*E_p2_i/ss)
          g_new.append(gg)
        g_new /= g_new[0]

        # Compute the log likelihood
        logL_new = \
          flex.sum(-0.5*(I-c*g*E_p)**2/V) + \
          flex.sum(-0.5*(E_p-m)**2/s) - \
          flex.sum(flex.log(s))

        # assert logL > logL_old
        eps = abs(logL_new-logL_old)
        print \
          iteration, \
          eps, \
          flex.sum(c_new)/len(c_new), \
          flex.min(g_new), \
          flex.max(g_new), \
          " ".join("%.2f" % mm for mm in m_new), \
          " ".join("%.2f" % sqrt(ss) for ss in s_new)
        if eps < self.tolerance:
          break
        g_old = g_new
        m_old = m_new
        s_old = s_new
        c_old = c_new
        logL_old = logL_new
      except KeyboardInterrupt, e:
        break
    
    # Set the results
    self.I_mean = c_new
    self.image_scales = g_new
    self.epsilon_mean = m_new
    self.epsilon_variance = s_new

    # Set some reflection data
    self.reflections['E_p'] = E_p
    self.reflections['E_p2'] = E_p2
    self.reflections['image_scale'] = g
    self.reflections['epsilon_mean'] = m
    self.reflections['epsilon_variance'] = s
    self.reflections['Imean'] = c


def scale_and_merge(experiments, reflections, params):

  # Configure the scaler
  scaler = Scaler(
    tolerance = params.tolerance,
    max_iter  = params.max_iter,
    grid_size = params.grid_size)

  # Run the scaler
  scaler.scale(experiments, reflections)

  # Get the reflections
  reflections = scaler.reflections

  # Get the results
  h_unique = scaler.h_unique
  I_mean = scaler.I_mean
  V_mean = scaler.V_mean
 
  # Get some scales
  image_scale = scaler.image_scale
  epsilon_mean = scaler.epsilon_mean
  epsilon_variance = scaler.epsilon_variance

  result = {
    'h_unique' : h_unique,
    'I_mean' : I_mean,
    'image_scale' : image_scale,
    'epsilon_mean' : epsilon_mean,
    'epsilon_variance' : epsilon_variance
  }

  reflections.as_pickle("output_reflections.pickle")

  pickle.dump(result, open('scales.pickle', 'wb'))


