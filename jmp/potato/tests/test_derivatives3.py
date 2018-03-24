from __future__ import division
from scitbx import matrix
from math import log, pi
from random import uniform, randint
from dials_scratch.jmp.potato.model import compute_change_of_basis_operation
from dials_scratch.jmp.potato.parameterisation import SimpleMosaicityParameterisation
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.refinement.parameterisation.crystal_parameters import CrystalUnitCellParameterisation
from dials.algorithms.refinement.parameterisation.crystal_parameters import CrystalOrientationParameterisation
from dials_scratch.jmp.potato.profile_refiner import ReflectionData
from dials.array_family import flex

def first_derivative(func, x, h):
  return (-func(x+2*h)+8*func(x+h)-8*func(x-h)+func(x-2*h)) / (12*h)

def second_derivative(func, x, y=None, h=None):
  if y is None:
    A = func(x+2*h)
    B = func(x+h)
    C = func(x)
    D = func(x-h)
    E = func(x-2*h)
    return (-(1/12)*(A+E) + (4/3)*(B+D) -(5/2)*C) / h**2
  else:
    A = func(x-h,y-h)
    B = func(x-h,y)
    C = func(x,y-h)
    D = func(x,y)
    E = func(x,y+h)
    F = func(x+h,y)
    G = func(x+h,y+h)
    return (A-B-C+2*D-E-F+G)/(2*h**2)


def generate_data(experiments, reflections):

  from random import seed

  seed(0)

  index = randint(0, len(reflections))

  h = reflections[index]['miller_index']

  s0 = matrix.col(experiments[0].beam.get_s0())

  U_param = CrystalOrientationParameterisation(experiments[0].crystal)
  B_param = CrystalUnitCellParameterisation(experiments[0].crystal)

  U = matrix.sqr(experiments[0].crystal.get_U())
  B = matrix.sqr(experiments[0].crystal.get_B())
  r = U*B*matrix.col(h)
  s2 = s0 + r
  mobs = (s2 +matrix.col((
    uniform(0, 1e-3),
    uniform(0, 1e-3),
    uniform(0, 1e-3)))).normalize()*s0.length()


  b1, b2, b3, b4, b5, b6 = (
    uniform(1e-3,3e-3),
    uniform(0.0,1e-3),
    uniform(1e-3,3e-3),
    uniform(0.0,1e-3),
    uniform(0.0,1e-3),
    uniform(1e-3,3e-3))

  params = (b1, b2, b3, b4, b5, b6)

  S_param = SimpleMosaicityParameterisation(params)

  ctot = randint(100,1000)

  T = matrix.sqr((
    uniform(1e-3,2e-3), 0,
    uniform(1e-6,2e-6), uniform(1e-3,2e-3)))
  Sobs = T*T.transpose()

  params = [S_param, U_param, B_param]


  return params, s0, h, ctot, mobs, Sobs


def test_first_derivatives(experiment, models, s0, h, ctot, mobs, Sobs):

  def compute_L(parameters):
    S_model = models[0]
    U_model = models[1]
    B_model = models[2]

    NU = U_model.num_total()
    NB = B_model.num_total()
    S_model.set_parameters(parameters[0:6])
    U_model.set_param_vals(parameters[6:6+NU])
    B_model.set_param_vals(parameters[6+NU:])

    U = matrix.sqr(experiment.crystal.get_U())
    B = matrix.sqr(experiment.crystal.get_B())

    q = U*B*h
    s2 = s0 + q

    R = compute_change_of_basis_operation(s0, mobs)

    sigma = S_model.sigma()
    sigmap= R*sigma*R.transpose()
    sigma11 = matrix.sqr((
      sigmap[0], sigmap[1],
      sigmap[3], sigmap[4]))
    sigma12 = matrix.col((sigmap[2], sigmap[5]))
    sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
    sigma22 = sigmap[8]

    mu = R*s2
    mu1 = matrix.col((mu[0], mu[1]))
    mu2 = mu[2]

    z = s0.length()
    mubar = mu1 + sigma12*(1/sigma22)*(z - mu2)
    sigma_bar = sigma11 - sigma12*(1/sigma22)*sigma21

    d = z-mu2
    c_d = mubar - matrix.col((0,0))
    A = log(sigma22)
    B = (1/sigma22)*d**2
    C = log(sigma_bar.determinant())*ctot
    D = (sigma_bar.inverse() * ctot*Sobs).trace()
    E = (sigma_bar.inverse() * ctot*c_d*c_d.transpose()).trace()
    return -0.5 * (A + B + C + D + E)


  step = 1e-8

  parameters = []
  parameters.extend(models[0].parameters())
  parameters.extend(models[1].get_param_vals())
  parameters.extend(models[2].get_param_vals())

  dL_num = []
  for i in range(len(parameters)):
    def f(x):
      p = [pp for pp in parameters]
      p[i] = x
      return compute_L(p)
    dL_num.append(first_derivative(f, parameters[i], step))

  def compute_dL_dS(i):

    S_model = models[0]
    U_model = models[1]
    B_model = models[2]
    NU = U_model.num_total()
    NB = B_model.num_total()

    U = matrix.sqr(experiment.crystal.get_U())
    B = matrix.sqr(experiment.crystal.get_B())

    q = U*B*h
    s2 = s0 + q

    reflection_model = ReflectionData(
      S_model,
      s0,
      s2,
      ctot,
      mobs,
      Sobs)
    return reflection_model.first_derivatives()[i]

  def compute_dL_dU(i, U_or_B):

    S_model = models[0]
    U_model = models[1]
    B_model = models[2]
    NU = U_model.num_total()
    NB = B_model.num_total()

    U = matrix.sqr(experiment.crystal.get_U())
    B = matrix.sqr(experiment.crystal.get_B())

    q = U*B*h
    s2 = s0 + q

    R = compute_change_of_basis_operation(s0, mobs)
    sigma = S_model.sigma()
    sigmap= R*sigma*R.transpose()
    S11 = matrix.sqr((
      sigmap[0], sigmap[1],
      sigmap[3], sigmap[4]))
    S12 = matrix.col((sigmap[2], sigmap[5]))
    S21 = matrix.col((sigmap[6], sigmap[7])).transpose()
    S22 = sigmap[8]

    mu = R*s2
    mu1 = matrix.col((mu[0], mu[1]))
    mu2 = mu[2]

    S22_inv = 1/S22

    epsilon = s0.length() - mu2
    mubar = mu1 + S12*(1/S22)*epsilon
    Sbar = S11 - S12*(1/S22)*S21

    Sbar_inv = Sbar.inverse()

    if U_or_B == "U":
      dU = U_model.get_ds_dp()[i]
      ds2 = dU*B*h
    else:
      dB = B_model.get_ds_dp()[i]
      ds2 = U*dB*h
    dmu = R*ds2
    dmu1 = matrix.col((dmu[0], dmu[1]))
    dmu2 = dmu[2]
    dep = -dmu2

    dmbar = dmu1 + S12*(1/S22)*dep
    dSbar = matrix.sqr((0, 0, 0, 0))
    dS22 = 0


    c_d = matrix.col((0,0)) - mubar
    I = matrix.sqr((
      1, 0,
      0, 1))

    U = S22_inv*dS22*(1 - S22_inv*epsilon**2)+2*S22_inv*epsilon*dep
    V = (Sbar_inv*dSbar*ctot*(I - Sbar_inv*(Sobs+c_d*c_d.transpose()))).trace()
    W = (-2*ctot*Sbar_inv*c_d*dmbar.transpose()).trace()

    return -0.5*(U+V+W)


  dL_cal = []
  # for i in range(6):
  #   dL_cal.append(compute_dL_dS(i))
  for i in range(models[1].num_total()):
    dL_cal.append(compute_dL_dU(i, U_or_B="U"))
  for i in range(models[2].num_total()):
    dL_cal.append(compute_dL_dU(i, U_or_B="B"))

  print dL_num
  print dL_cal

def read_experiments():
  experiments = ExperimentListFactory.from_json_file("experiments.json")
  return experiments


def test():

  experiments = read_experiments()
  experiments[0].scan.set_oscillation((0,0.01), deg=True)
  reflections = flex.reflection_table.from_predictions_multi(experiments)

  for i in range(1):

    models, s0, h, ctot, mobs, Sobs = generate_data(experiments, reflections)

    test_first_derivatives(experiments[0], models, s0, h, ctot, mobs, Sobs)

if __name__ == '__main__':
  test()
