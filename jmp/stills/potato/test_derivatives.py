from __future__ import division
from scitbx import matrix
from math import log, pi
from random import uniform, randint

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


def compute_sigma(b1, b2, b3, b4, b5, b6):
  M = matrix.sqr((
    b1, 0, 0,
    b2, b3, 0,
    b4, b5, b6))
  sigma = M*M.transpose()
  return sigma

def compute_dSdb1(b1, b2, b3, b4, b5, b6):
  return matrix.sqr((
    2*b1,b2,b4,
    b2,0,0,
    b4,0,0))

def compute_dSdb2(b1, b2, b3, b4, b5, b6):
  return matrix.sqr((
    0,b1,0,
    b1,2*b2,b4,
    0,b4,0))

def compute_dSdb3(b1, b2, b3, b4, b5, b6):
  return matrix.sqr((
    0,0,0,
    0,2*b3,b5,
    0,b5,0))

def compute_dSdb4(b1, b2, b3, b4, b5, b6):
  return matrix.sqr((
    0,0,b1,
    0,0,b2,
    b1,b2,2*b4))

def compute_dSdb5(b1, b2, b3, b4, b5, b6):
  return matrix.sqr((
    0,0,0,
    0,0,b3,
    0,b3,2*b5))

def compute_dSdb6(b1, b2, b3, b4, b5, b6):
  return matrix.sqr((
    0,0,0,
    0,0,0,
    0,0,2*b6))


def compute_d2Sdbij(i, j, b1, b2, b3, b4, b5, b6):
  data = {
    1 : [
      matrix.sqr((
        2,0,0,
        0,0,0,
        0,0,0)),
      matrix.sqr((
        0,1,0,
        1,0,0,
        0,0,0)),
      None,
      matrix.sqr((
        0, 0, 1,
        0, 0, 0,
        1, 0, 0)),
      None,
      None
    ],
    2 : [
      matrix.sqr((
        0, 1, 0,
        1, 0, 0,
        0, 0, 0)),
      matrix.sqr((
        0, 0, 0,
        0, 2, 0,
        0, 0, 0)),
      None,
      matrix.sqr((
        0, 0, 0,
        0, 0, 1,
        0, 1, 0)),
      None,
      None
    ],
    3 : [
      None,
      None,
      matrix.sqr((
        0, 0, 0,
        0, 2, 0,
        0, 0, 0)),
      None,
      matrix.sqr((
        0, 0, 0,
        0, 0, 1,
        0, 1, 0)),
      None
    ],
    4 : [
      matrix.sqr((
        0, 0, 1,
        0, 0, 0,
        1, 0, 0)),
      matrix.sqr((
        0, 0, 0,
        0, 0, 1,
        0, 1, 0)),
      None,
      matrix.sqr((
        0, 0, 0,
        0, 0, 0,
        0, 0, 2)),
      None,
      None
    ],
    5 : [
      None,
      None,
      matrix.sqr((
        0, 0, 0,
        0, 0, 1,
        0, 1, 0)),
      None,
      matrix.sqr((
        0, 0, 0,
        0, 0, 0,
        0, 0, 2)),
      None
    ],
    6 : [
      None,
      None,
      None,
      None,
      None,
      matrix.sqr((
        0, 0, 0,
        0, 0, 0,
        0, 0, 2))
    ],
  }

  m = data[i+1][j]

  if m is None:
    m = matrix.sqr((0, 0, 0, 0, 0, 0, 0, 0, 0))

  return m

def generate_data():

  from random import seed

  seed(0)

  axis = matrix.col((
    uniform(-1,1),
    uniform(-1,1),
    uniform(-1,1)))
  angle=uniform(-pi/2,pi/2)

  R = axis.axis_and_angle_as_r3_rotation_matrix(angle=angle,deg=False)
  mu2 = uniform(0.8,1.2)
  r = 1

  T = matrix.sqr((
    uniform(0.5,2.5), 0,
    uniform(0,0.1), uniform(0.5,2.5)))
  S = T*T.transpose()

  #b1, b2, b3, b4, b5, b6 = 1, 0.1, 2, 0.2, 0.3, 3
  b1, b2, b3, b4, b5, b6 = (
    uniform(0.5,3.5),
    uniform(0.0,0.5),
    uniform(0.5,3.5),
    uniform(0.0,0.5),
    uniform(0.0,0.5),
    uniform(0.5,3.5))

  ctot = randint(100,1000)

  return (b1, b2, b3, b4, b5, b6), R, mu2, r, S, ctot

def test_dSdb_22():

  (b1, b2, b3, b4, b5, b6), R, mu2, r, Sobs, ctot = generate_data()

  def compute_sigma22(b1, b2, b3, b4, b5, b6):
    sigma = compute_sigma(b1, b2, b3, b4, b5, b6)
    sigmap= R*sigma*R.transpose()
    sigma22 = sigmap[8]
    return sigma22

  def f1(x):
    return compute_sigma22(x, b2, b3, b4, b5, b6)

  def f2(x):
    return compute_sigma22(b1, x, b3, b4, b5, b6)

  def f3(x):
    return compute_sigma22(b1, b2, x, b4, b5, b6)

  def f4(x):
    return compute_sigma22(b1, b2, b3, x, b5, b6)

  def f5(x):
    return compute_sigma22(b1, b2, b3, b4, x, b6)

  def f6(x):
    return compute_sigma22(b1, b2, b3, b4, b5, x)

  h = 0.001

  dSdb1_22_num = first_derivative(f1, b1, h)
  dSdb2_22_num = first_derivative(f2, b2, h)
  dSdb3_22_num = first_derivative(f3, b3, h)
  dSdb4_22_num = first_derivative(f4, b4, h)
  dSdb5_22_num = first_derivative(f5, b5, h)
  dSdb6_22_num = first_derivative(f6, b6, h)

  dSdb1 = compute_dSdb1(b1, b2, b3, b4, b5, b6)
  dSdb2 = compute_dSdb2(b1, b2, b3, b4, b5, b6)
  dSdb3 = compute_dSdb3(b1, b2, b3, b4, b5, b6)
  dSdb4 = compute_dSdb4(b1, b2, b3, b4, b5, b6)
  dSdb5 = compute_dSdb5(b1, b2, b3, b4, b5, b6)
  dSdb6 = compute_dSdb6(b1, b2, b3, b4, b5, b6)

  dSdb1_22_cal = (R*dSdb1*R.transpose())[8]
  dSdb2_22_cal = (R*dSdb2*R.transpose())[8]
  dSdb3_22_cal = (R*dSdb3*R.transpose())[8]
  dSdb4_22_cal = (R*dSdb4*R.transpose())[8]
  dSdb5_22_cal = (R*dSdb5*R.transpose())[8]
  dSdb6_22_cal = (R*dSdb6*R.transpose())[8]

  assert abs(dSdb1_22_num - dSdb1_22_cal) < 1e-7
  assert abs(dSdb2_22_num - dSdb2_22_cal) < 1e-7
  assert abs(dSdb3_22_num - dSdb3_22_cal) < 1e-7
  assert abs(dSdb4_22_num - dSdb4_22_cal) < 1e-7
  assert abs(dSdb5_22_num - dSdb5_22_cal) < 1e-7
  assert abs(dSdb6_22_num - dSdb6_22_cal) < 1e-7

  #print 'OK'

def test_dS_bar_db():

  (b1, b2, b3, b4, b5, b6), R, mu2, r, Sobs, ctot = generate_data()

  def compute_sigma_bar(b1, b2, b3, b4, b5, b6):
    sigma = compute_sigma(b1, b2, b3, b4, b5, b6)
    sigmap= R*sigma*R.transpose()
    sigma11 = matrix.sqr((
      sigmap[0], sigmap[1],
      sigmap[3], sigmap[4]))
    sigma12 = matrix.col((sigmap[2], sigmap[5]))
    sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
    sigma22 = sigmap[8]
    return sigma11 - sigma12*(1/sigma22)*sigma21

  def f1(x):
    return compute_sigma_bar(x, b2, b3, b4, b5, b6)

  def f2(x):
    return compute_sigma_bar(b1, x, b3, b4, b5, b6)

  def f3(x):
    return compute_sigma_bar(b1, b2, x, b4, b5, b6)

  def f4(x):
    return compute_sigma_bar(b1, b2, b3, x, b5, b6)

  def f5(x):
    return compute_sigma_bar(b1, b2, b3, b4, x, b6)

  def f6(x):
    return compute_sigma_bar(b1, b2, b3, b4, b5, x)

  h = 0.001

  dS_bar_db1_num = first_derivative(f1, b1, h)
  dS_bar_db2_num = first_derivative(f2, b2, h)
  dS_bar_db3_num = first_derivative(f3, b3, h)
  dS_bar_db4_num = first_derivative(f4, b4, h)
  dS_bar_db5_num = first_derivative(f5, b5, h)
  dS_bar_db6_num = first_derivative(f6, b6, h)

  def compute_dS_bar_db(S, dS):
    RSR = R*S*R.transpose()
    RdSR = R*dS*R.transpose()

    S11 = matrix.sqr((
      RSR[0], RSR[1],
      RSR[3], RSR[4]))
    S12 = matrix.col((RSR[2], RSR[5]))
    S21 = matrix.col((RSR[6], RSR[7])).transpose()
    S22 = RSR[8]

    dS11 = matrix.sqr((
      RdSR[0], RdSR[1],
      RdSR[3], RdSR[4]))
    dS12 = matrix.col((RdSR[2], RdSR[5]))
    dS21 = matrix.col((RdSR[6], RdSR[7])).transpose()
    dS22 = RdSR[8]

    A = dS11
    B = S12*(1/S22)*dS22*(1/S22)*S21
    C = S12*(1/S22)*dS21
    D = dS12*(1/S22)*S21
    return A + B - (C+D)

  dSdb1 = compute_dSdb1(b1, b2, b3, b4, b5, b6)
  dSdb2 = compute_dSdb2(b1, b2, b3, b4, b5, b6)
  dSdb3 = compute_dSdb3(b1, b2, b3, b4, b5, b6)
  dSdb4 = compute_dSdb4(b1, b2, b3, b4, b5, b6)
  dSdb5 = compute_dSdb5(b1, b2, b3, b4, b5, b6)
  dSdb6 = compute_dSdb6(b1, b2, b3, b4, b5, b6)

  S = compute_sigma(b1, b2, b3, b4, b5, b6)

  dS_bar_db1_cal = compute_dS_bar_db(S, dSdb1)
  dS_bar_db2_cal = compute_dS_bar_db(S, dSdb2)
  dS_bar_db3_cal = compute_dS_bar_db(S, dSdb3)
  dS_bar_db4_cal = compute_dS_bar_db(S, dSdb4)
  dS_bar_db5_cal = compute_dS_bar_db(S, dSdb5)
  dS_bar_db6_cal = compute_dS_bar_db(S, dSdb6)

  assert all(abs(a-b) < 1e-7 for a, b in zip(dS_bar_db1_num, dS_bar_db1_cal))
  assert all(abs(a-b) < 1e-7 for a, b in zip(dS_bar_db2_num, dS_bar_db2_cal))
  assert all(abs(a-b) < 1e-7 for a, b in zip(dS_bar_db3_num, dS_bar_db3_cal))
  assert all(abs(a-b) < 1e-7 for a, b in zip(dS_bar_db4_num, dS_bar_db4_cal))
  assert all(abs(a-b) < 1e-7 for a, b in zip(dS_bar_db5_num, dS_bar_db5_cal))
  assert all(abs(a-b) < 1e-7 for a, b in zip(dS_bar_db6_num, dS_bar_db6_cal))

  #print 'OK'

def test_dLdb():

  (b1, b2, b3, b4, b5, b6), R, mu2, r, Sobs, ctot = generate_data()

  def compute_L(b1, b2, b3, b4, b5, b6):
    sigma = compute_sigma(b1, b2, b3, b4, b5, b6)
    sigmap= R*sigma*R.transpose()
    sigma11 = matrix.sqr((
      sigmap[0], sigmap[1],
      sigmap[3], sigmap[4]))
    sigma12 = matrix.col((sigmap[2], sigmap[5]))
    sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
    sigma22 = sigmap[8]

    sigma_bar = sigma11 - sigma12*(1/sigma22)*sigma21

    d = r-mu2
    A = log(sigma22)
    B = (1/sigma22)*d**2
    C = log(sigma_bar.determinant())*ctot
    D = (sigma_bar.inverse() * Sobs).trace()
    return -0.5 * (A + B + C + D)


  def f1(x):
    return compute_L(x, b2, b3, b4, b5, b6)

  def f2(x):
    return compute_L(b1, x, b3, b4, b5, b6)

  def f3(x):
    return compute_L(b1, b2, x, b4, b5, b6)

  def f4(x):
    return compute_L(b1, b2, b3, x, b5, b6)

  def f5(x):
    return compute_L(b1, b2, b3, b4, x, b6)

  def f6(x):
    return compute_L(b1, b2, b3, b4, b5, x)

  h = 0.001

  dLdb1_num = first_derivative(f1, b1, h)
  dLdb2_num = first_derivative(f2, b2, h)
  dLdb3_num = first_derivative(f3, b3, h)
  dLdb4_num = first_derivative(f4, b4, h)
  dLdb5_num = first_derivative(f5, b5, h)
  dLdb6_num = first_derivative(f6, b6, h)

  def compute_dLdb(S, dS):
    RSR = R*S*R.transpose()
    RdSR = R*dS*R.transpose()

    S11 = matrix.sqr((
      RSR[0], RSR[1],
      RSR[3], RSR[4]))
    S12 = matrix.col((RSR[2], RSR[5]))
    S21 = matrix.col((RSR[6], RSR[7])).transpose()
    S22 = RSR[8]

    Sbar = S11 - S12*(1/S22)*S21

    dS11 = matrix.sqr((
      RdSR[0], RdSR[1],
      RdSR[3], RdSR[4]))
    dS12 = matrix.col((RdSR[2], RdSR[5]))
    dS21 = matrix.col((RdSR[6], RdSR[7])).transpose()
    dS22 = RdSR[8]

    A = dS11
    B = dS12*(1/S22)*S21
    C = S12*(1/S22)*(dS22*(1/S22)*S21 - dS21)
    dSbar = A - B + C

    d = r-mu2
    I = matrix.sqr((
      1, 0,
      0, 1))
    U = (1/S22)*dS22*(1-(1/S22)*d**2)
    V = (Sbar.inverse()*dSbar*(ctot*I-Sbar.inverse()*Sobs)).trace()
    return -0.5 * (U+V)


  dSdb1 = compute_dSdb1(b1, b2, b3, b4, b5, b6)
  dSdb2 = compute_dSdb2(b1, b2, b3, b4, b5, b6)
  dSdb3 = compute_dSdb3(b1, b2, b3, b4, b5, b6)
  dSdb4 = compute_dSdb4(b1, b2, b3, b4, b5, b6)
  dSdb5 = compute_dSdb5(b1, b2, b3, b4, b5, b6)
  dSdb6 = compute_dSdb6(b1, b2, b3, b4, b5, b6)

  S = compute_sigma(b1, b2, b3, b4, b5, b6)

  dLdb1_cal = compute_dLdb(S, dSdb1)
  dLdb2_cal = compute_dLdb(S, dSdb2)
  dLdb3_cal = compute_dLdb(S, dSdb3)
  dLdb4_cal = compute_dLdb(S, dSdb4)
  dLdb5_cal = compute_dLdb(S, dSdb5)
  dLdb6_cal = compute_dLdb(S, dSdb6)

  assert abs(dLdb1_num - dLdb1_cal) < 1e-7
  assert abs(dLdb2_num - dLdb2_cal) < 1e-7
  assert abs(dLdb3_num - dLdb3_cal) < 1e-7
  assert abs(dLdb4_num - dLdb4_cal) < 1e-7
  assert abs(dLdb5_num - dLdb5_cal) < 1e-7
  assert abs(dLdb6_num - dLdb6_cal) < 1e-7

  #print 'OK'

def test_d2S_dbij(i, j):

  (b1, b2, b3, b4, b5, b6), R, mu2, r, Sobs, ctot = generate_data()

  def compute_S(b1, b2, b3, b4, b5, b6):
    sigma = compute_sigma(b1, b2, b3, b4, b5, b6)
    sigmap= R*sigma*R.transpose()
    sigma11 = matrix.sqr((
      sigmap[0], sigmap[1],
      sigmap[3], sigmap[4]))
    sigma12 = matrix.col((sigmap[2], sigmap[5]))
    sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
    sigma22 = sigmap[8]
    return sigma22

  def f1(x):
    params = [b1, b2, b3, b4, b5, b6]
    params[i] = x
    return compute_S(*params)

  def f2(x, y):
    params = [b1, b2, b3, b4, b5, b6]
    params[i] = x
    params[j] = y
    return compute_S(*params)

  h = 0.001

  x = (b1, b2, b3, b4, b5, b6)[i]
  y = (b1, b2, b3, b4, b5, b6)[j]

  if i == j:
    d2Sdb_22_num = second_derivative(f1, x=x, h=h)
  else:
    d2Sdb_22_num = second_derivative(f2, x=x, y=y, h=h)


  d2Sdb = compute_d2Sdbij(i, j, b1, b2, b3, b4, b5, b6)

  d2Sdb_22_cal = (R*d2Sdb*R.transpose())[8]

  assert abs(d2Sdb_22_num - d2Sdb_22_cal) < 1e-7

  #print 'OK'

def test_d2S_bar_dbij(i, j):

  (b1, b2, b3, b4, b5, b6), R, mu2, r, Sobs, ctot = generate_data()

  def compute_Sbar(b1, b2, b3, b4, b5, b6):
    sigma = compute_sigma(b1, b2, b3, b4, b5, b6)
    sigmap= R*sigma*R.transpose()
    sigma11 = matrix.sqr((
      sigmap[0], sigmap[1],
      sigmap[3], sigmap[4]))
    sigma12 = matrix.col((sigmap[2], sigmap[5]))
    sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
    sigma22 = sigmap[8]
    return sigma11 - sigma12*(1/sigma22)*sigma21

  def f1(x):
    params = [b1, b2, b3, b4, b5, b6]
    params[i] = x
    return compute_Sbar(*params)

  def f2(x, y):
    params = [b1, b2, b3, b4, b5, b6]
    params[i] = x
    params[j] = y
    return compute_Sbar(*params)

  h = 0.001

  x = (b1, b2, b3, b4, b5, b6)[i]
  y = (b1, b2, b3, b4, b5, b6)[j]

  if i == j:
    d2S_bar_db_num = second_derivative(f1, x=x, h=h)
  else:
    d2S_bar_db_num = second_derivative(f2, x=x, y=y, h=h)

  def compute_d2S_bar_dbij(S, dSi, dSj, d2S):

    RSR = R*S*R.transpose()
    RdSiR = R*dSi*R.transpose()
    RdSjR = R*dSj*R.transpose()
    Rd2SR = R*d2S*R.transpose()

    S11 = matrix.sqr((
      RSR[0], RSR[1],
      RSR[3], RSR[4]))
    S12 = matrix.col((RSR[2], RSR[5]))
    S21 = matrix.col((RSR[6], RSR[7])).transpose()
    S22 = RSR[8]

    dSi11 = matrix.sqr((
      RdSiR[0], RdSiR[1],
      RdSiR[3], RdSiR[4]))
    dSi12 = matrix.col((RdSiR[2], RdSiR[5]))
    dSi21 = matrix.col((RdSiR[6], RdSiR[7])).transpose()
    dSi22 = RdSiR[8]

    dSj11 = matrix.sqr((
      RdSjR[0], RdSjR[1],
      RdSjR[3], RdSjR[4]))
    dSj12 = matrix.col((RdSjR[2], RdSjR[5]))
    dSj21 = matrix.col((RdSjR[6], RdSjR[7])).transpose()
    dSj22 = RdSjR[8]

    d2S11 = matrix.sqr((
      Rd2SR[0], Rd2SR[1],
      Rd2SR[3], Rd2SR[4]))
    d2S12 = matrix.col((Rd2SR[2], Rd2SR[5]))
    d2S21 = matrix.col((Rd2SR[6], Rd2SR[7])).transpose()
    d2S22 = Rd2SR[8]

    A = d2S11
    B = dSj12*(1/S22)*dSi22*(1/S22)*S21
    C = S12*(1/S22)*dSj22*(1/S22)*dSi22*(1/S22)*S21

    D = S12*(1/S22)*d2S22*(1/S22)*S21
    E = S12*(1/S22)*dSi22*(1/S22)*dSj22*(1/S22)*S21
    F = S12*(1/S22)*dSi22*(1/S22)*dSj21

    G = dSj12*(1/S22)*dSi21
    H = S12*(1/S22)*dSj22*(1/S22)*dSi21
    I = S12*(1/S22)*d2S21

    J = d2S12*(1/S22)*S21
    K = dSi12*(1/S22)*(dSj22)*(1/S22)*S21
    L = dSi12*(1/S22)*dSj21

    return A+B-C+D-E+F-G+H-I-J+K-L

  S = compute_sigma(b1, b2, b3, b4, b5, b6)

  dS_func = {
    0 : compute_dSdb1(b1, b2, b3, b4, b5, b6),
    1 : compute_dSdb2(b1, b2, b3, b4, b5, b6),
    2 : compute_dSdb3(b1, b2, b3, b4, b5, b6),
    3 : compute_dSdb4(b1, b2, b3, b4, b5, b6),
    4 : compute_dSdb5(b1, b2, b3, b4, b5, b6),
    5 : compute_dSdb6(b1, b2, b3, b4, b5, b6)
  }

  dSi = dS_func[i]
  dSj = dS_func[j]

  d2Sdb = compute_d2Sdbij(i, j, b1, b2, b3, b4, b5, b6)

  d2S_bar_db_cal = compute_d2S_bar_dbij(S, dSi, dSj, d2Sdb)

  assert all(abs(a-b) < 1e-5 for a, b in zip(d2S_bar_db_num, d2S_bar_db_cal))

  #print 'OK'


def test_d2L_dbij(i, j):

  (b1, b2, b3, b4, b5, b6), R, mu2, r, Sobs, ctot = generate_data()

  def compute_L(b1, b2, b3, b4, b5, b6):
    sigma = compute_sigma(b1, b2, b3, b4, b5, b6)
    sigmap= R*sigma*R.transpose()
    sigma11 = matrix.sqr((
      sigmap[0], sigmap[1],
      sigmap[3], sigmap[4]))
    sigma12 = matrix.col((sigmap[2], sigmap[5]))
    sigma21 = matrix.col((sigmap[6], sigmap[7])).transpose()
    sigma22 = sigmap[8]

    sigma_bar = sigma11 - sigma12*(1/sigma22)*sigma21

    d = r-mu2
    A = log(sigma22)
    B = (1/sigma22)*d**2
    C = log(sigma_bar.determinant())*ctot
    D = (sigma_bar.inverse() * Sobs).trace()
    return -0.5 * (A + B + C + D)

  def f1(x):
    params = [b1, b2, b3, b4, b5, b6]
    params[i] = x
    return compute_L(*params)

  def f2(x, y):
    params = [b1, b2, b3, b4, b5, b6]
    params[i] = x
    params[j] = y
    return compute_L(*params)

  h = 0.001

  x = (b1, b2, b3, b4, b5, b6)[i]
  y = (b1, b2, b3, b4, b5, b6)[j]

  if i == j:
    d2Ldb_num = second_derivative(f1, x=x, h=h)
  else:
    d2Ldb_num = second_derivative(f2, x=x, y=y, h=h)

  def compute_dS_bar_db(S, dS):
    RSR = R*S*R.transpose()
    RdSR = R*dS*R.transpose()

    S11 = matrix.sqr((
      RSR[0], RSR[1],
      RSR[3], RSR[4]))
    S12 = matrix.col((RSR[2], RSR[5]))
    S21 = matrix.col((RSR[6], RSR[7])).transpose()
    S22 = RSR[8]

    dS11 = matrix.sqr((
      RdSR[0], RdSR[1],
      RdSR[3], RdSR[4]))
    dS12 = matrix.col((RdSR[2], RdSR[5]))
    dS21 = matrix.col((RdSR[6], RdSR[7])).transpose()
    dS22 = RdSR[8]

    A = dS11
    B = S12*(1/S22)*dS22*(1/S22)*S21
    C = S12*(1/S22)*dS21
    D = dS12*(1/S22)*S21
    return A + B - (C+D)

  def compute_d2S_bar_dbij(S, dSi, dSj, d2S):

    RSR = R*S*R.transpose()
    RdSiR = R*dSi*R.transpose()
    RdSjR = R*dSj*R.transpose()
    Rd2SR = R*d2S*R.transpose()

    S11 = matrix.sqr((
      RSR[0], RSR[1],
      RSR[3], RSR[4]))
    S12 = matrix.col((RSR[2], RSR[5]))
    S21 = matrix.col((RSR[6], RSR[7])).transpose()
    S22 = RSR[8]

    dSi11 = matrix.sqr((
      RdSiR[0], RdSiR[1],
      RdSiR[3], RdSiR[4]))
    dSi12 = matrix.col((RdSiR[2], RdSiR[5]))
    dSi21 = matrix.col((RdSiR[6], RdSiR[7])).transpose()
    dSi22 = RdSiR[8]

    dSj11 = matrix.sqr((
      RdSjR[0], RdSjR[1],
      RdSjR[3], RdSjR[4]))
    dSj12 = matrix.col((RdSjR[2], RdSjR[5]))
    dSj21 = matrix.col((RdSjR[6], RdSjR[7])).transpose()
    dSj22 = RdSjR[8]

    d2S11 = matrix.sqr((
      Rd2SR[0], Rd2SR[1],
      Rd2SR[3], Rd2SR[4]))
    d2S12 = matrix.col((Rd2SR[2], Rd2SR[5]))
    d2S21 = matrix.col((Rd2SR[6], Rd2SR[7])).transpose()
    d2S22 = Rd2SR[8]

    A = d2S11
    B = dSj12*(1/S22)*dSi22*(1/S22)*S21
    C = S12*(1/S22)*dSj22*(1/S22)*dSi22*(1/S22)*S21

    D = S12*(1/S22)*d2S22*(1/S22)*S21
    E = S12*(1/S22)*dSi22*(1/S22)*dSj22*(1/S22)*S21
    F = S12*(1/S22)*dSi22*(1/S22)*dSj21

    G = dSj12*(1/S22)*dSi21
    H = S12*(1/S22)*dSj22*(1/S22)*dSi21
    I = S12*(1/S22)*d2S21

    J = d2S12*(1/S22)*S21
    K = dSi12*(1/S22)*(dSj22)*(1/S22)*S21
    L = dSi12*(1/S22)*dSj21

    return A+B-C+D-E+F-G+H-I-J+K-L

  def compute_dLdb(S, dSi, dSj, d2S):

    RSR = R*S*R.transpose()
    RdSiR = R*dSi*R.transpose()
    RdSjR = R*dSj*R.transpose()
    Rd2SR = R*d2S*R.transpose()

    S11 = matrix.sqr((
      RSR[0], RSR[1],
      RSR[3], RSR[4]))
    S12 = matrix.col((RSR[2], RSR[5]))
    S21 = matrix.col((RSR[6], RSR[7])).transpose()
    S22 = RSR[8]

    dSi11 = matrix.sqr((
      RdSiR[0], RdSiR[1],
      RdSiR[3], RdSiR[4]))
    dSi12 = matrix.col((RdSiR[2], RdSiR[5]))
    dSi21 = matrix.col((RdSiR[6], RdSiR[7])).transpose()
    dSi22 = RdSiR[8]

    dSj11 = matrix.sqr((
      RdSjR[0], RdSjR[1],
      RdSjR[3], RdSjR[4]))
    dSj12 = matrix.col((RdSjR[2], RdSjR[5]))
    dSj21 = matrix.col((RdSjR[6], RdSjR[7])).transpose()
    dSj22 = RdSjR[8]

    d2S11 = matrix.sqr((
      Rd2SR[0], Rd2SR[1],
      Rd2SR[3], Rd2SR[4]))
    d2S12 = matrix.col((Rd2SR[2], Rd2SR[5]))
    d2S21 = matrix.col((Rd2SR[6], Rd2SR[7])).transpose()
    d2S22 = Rd2SR[8]

    d = r-mu2
    Sbar = S11 - S12*(1/S22)*S21
    dSbari = compute_dS_bar_db(S, dSi)
    dSbarj = compute_dS_bar_db(S, dSj)
    d2Sbar = compute_d2S_bar_dbij(S, dSi, dSj, d2S)

    I = matrix.sqr((
      1, 0,
      0, 1))

    A1 = (1/S22)*d2S22*(1 - (1/S22)*d**2)
    A2 = (1/S22)*dSj22*(1/S22)*dSi22*(1 - 2*(1/S22)*d**2)
    B1 = Sbar.inverse() * d2Sbar*(ctot*I - Sbar.inverse()*Sobs)
    B2 = Sbar.inverse() * dSbarj * Sbar.inverse() * dSbari * (ctot*I - 2*Sbar.inverse()*Sobs)
    A = A1-A2
    B = (B1-B2).trace()
    return -0.5 * (A+B)

  S = compute_sigma(b1, b2, b3, b4, b5, b6)

  dS_func = {
    0 : compute_dSdb1(b1, b2, b3, b4, b5, b6),
    1 : compute_dSdb2(b1, b2, b3, b4, b5, b6),
    2 : compute_dSdb3(b1, b2, b3, b4, b5, b6),
    3 : compute_dSdb4(b1, b2, b3, b4, b5, b6),
    4 : compute_dSdb5(b1, b2, b3, b4, b5, b6),
    5 : compute_dSdb6(b1, b2, b3, b4, b5, b6)
  }

  dSi = dS_func[i]
  dSj = dS_func[j]

  d2Sdb = compute_d2Sdbij(i, j, b1, b2, b3, b4, b5, b6)

  d2Ldb_cal = compute_dLdb(S, dSi, dSj, d2Sdb)

  assert abs(d2Ldb_num - d2Ldb_cal) < 1e-5

  #print 'OK'

def test_first_derivatives():
  test_dSdb_22()
  test_dS_bar_db()
  test_dLdb()


def test_second_derivatives():
  for j in range(6):
    for i in range(6):
      test_d2S_dbij(i,j)
      test_d2S_bar_dbij(i,j)
      test_d2L_dbij(i,j)

if __name__ == '__main__':

  for i in range(10):
    print i
    test_first_derivatives()
    test_second_derivatives()
