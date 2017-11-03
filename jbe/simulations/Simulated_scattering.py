from dials.array_family import flex


dimensions = (2,2,5)

def generate_crystal(dimensions):
  crystal_points = flex.vec3_double([])
  for i in range(dimensions[0]*10):
    for j in range(dimensions[1]*10):
      for k in range(dimensions[2]*10):
        crystal_points.append(((i+0.5)/10.0,(j+0.5)/10.0,(k+0.5)/10.0))
  return crystal_points

def crystal_planes(dimensions):
  #define each plane by (point in plane, plane normal)
  #define here a cuboidal crystal
  p1 = (dimensions, (1, 0, 0))
  p2 = ((0, 0, 0), (-1, 0, 0))
  p3 = (dimensions, (0, 1, 0))
  p4 = ((0, 0, 0), (0, -1, 0))
  p5 = (dimensions, (0, 0, 1))
  p6 = ((0, 0, 0), (0, 0, -1))
  return (p1, p2, p3, p4, p5, p6)

def calc_dist_to_surface(k, p0, planes):
  #calculate distance along direction of +k to edge of crystal defined by planes
  s_list = []
  for plane in planes:
    n_u = ((plane[1][0] * k[0]) + (plane[1][1] * k[1]) + (plane[1][2] * k[2]))
    if n_u != 0.0:
      s = -((plane[1][0] * (p0[0] - plane[0][0])) 
           +(plane[1][1] * (p0[1] - plane[0][1]))
           +(plane[1][2] * (p0[2] - plane[0][2])))/n_u
      if s > 0:
        s_list.append(s)
  return min(s_list)

def calc_path_through_crystal(k0, k1, p0, planes):
  k0_pr = (-1.0 * k0[0], -1.0 * k0[1], -1.0 * k0[2])
  d1 = calc_dist_to_surface(k0_pr, p0, planes)
  d2 = calc_dist_to_surface(k1, p0, planes)
  return d1 + d2

def normalise_vec(v):
  s = ((v[0]**2) + (v[1]**2) + (v[2]**2))**0.5
  return (v[0]/s, v[1]/s, v[2]/s)

def calculate_intensity_of_reflection(k0, k1, crystal_points, cryst_planes):
  from math import exp
  mu = 1.0
  I = 0.0
  k0 = normalise_vec(k0)
  k1 = normalise_vec(k1)
  for point in crystal_points:
    d = calc_path_through_crystal(k0, k1, point, cryst_planes)
    I += exp(- mu * d)
  return I

def rotation_matrix(phi):
  from math import cos, sin, pi
  return ((cos(pi*phi/180.0),-1.0*sin(pi*phi/180.0),0.0),
          (sin(pi*phi/180.0),cos(pi*phi/180.0),0.0),
          (0.0,0.0,1.0))

def multiply_matrices(M1, M2):
  a11 = (M1[0][0] * M2[0][0]) + (M1[0][1] * M2[1][0]) + (M1[0][2] * M2[2][0])
  a12 = (M1[0][0] * M2[0][1]) + (M1[0][1] * M2[1][1]) + (M1[0][2] * M2[2][1])
  a13 = (M1[0][0] * M2[0][2]) + (M1[0][1] * M2[1][2]) + (M1[0][2] * M2[2][2])
  a21 = (M1[1][0] * M2[0][0]) + (M1[1][1] * M2[1][0]) + (M1[1][2] * M2[2][0])
  a22 = (M1[1][0] * M2[0][1]) + (M1[1][1] * M2[1][1]) + (M1[1][2] * M2[2][1])
  a23 = (M1[1][0] * M2[0][2]) + (M1[1][1] * M2[1][2]) + (M1[1][2] * M2[2][2])
  a31 = (M1[2][0] * M2[0][0]) + (M1[2][1] * M2[1][0]) + (M1[2][2] * M2[2][0])
  a32 = (M1[2][0] * M2[0][1]) + (M1[2][1] * M2[1][1]) + (M1[2][2] * M2[2][1])
  a33 = (M1[2][0] * M2[0][2]) + (M1[2][1] * M2[1][2]) + (M1[2][2] * M2[2][2])
  return ((a11, a12, a13), (a21, a22, a23), (a31, a32, a33))

def rotate_vector(M, V):
  a11 = (M[0][0] * V[0]) + (M[0][1] * V[1]) + (M[0][2] * V[2])
  a21 = (M[1][0] * V[0]) + (M[1][1] * V[1]) + (M[1][2] * V[2])
  a31 = (M[2][0] * V[0]) + (M[2][1] * V[1]) + (M[2][2] * V[2])
  return (a11, a21, a31)

def transform_to_lab_frame(UB, R, millerset):
  RUB = multiply_matrices(R, UB)
  rotated_indices = flex.vec3_double([])
  for index in millerset.indices():
    rotated_indices.append(rotate_vector(RUB, index))
  return rotated_indices

def calculate_scattering_vectors(k0, rlps):
  scattering_vectors = flex.vec3_double([])
  for rlp in rlps:
    k1_0 = k0[0]+rlp[0]
    k1_1 = k0[1]+rlp[1]
    k1_2 = k0[2]+rlp[2]
    scattering_vectors.append((k1_0, k1_1, k1_2))
  return scattering_vectors

def calculate_spot_positions(scattering_vectors, ms):
  x_det = 10.0 #distance from crystal along beam direction
  y_det = 4.0 #distance detector extends in y_dir (size = 2*y_det)
  z_det = 4.0 #distance detector extends in z_dir (size = 2*z_det)
  spots_y = []
  spots_z = []
  indices = ms.indices()
  miller_indices = []
  spot_vectors = []
  for idx, vector in enumerate(scattering_vectors):
    if vector[0] > 0.0:
      vector_size = ((vector[0]**2) + (vector[1]**2) + (vector[2]**2))**0.5
      if 0.995 < vector_size < 1.005:
        y_at_det = x_det * vector[1]/vector[0]
        z_at_det = x_det * vector[2]/vector[0]
        if abs(y_at_det) < y_det and abs(z_at_det) < z_det:
          spots_y.append(y_at_det)
          spots_z.append(z_at_det)
          spot_vectors.append(vector)
          miller_indices.append(indices[idx])
  return spots_y, spots_z, miller_indices, spot_vectors


from cctbx import miller
from cctbx import crystal

def simulate_dataset():
  a = 4.0
  b = 4.0
  c = 6.0
  dimensions = (2,2,5)
  crystal_points = generate_crystal(dimensions)
  cryst_planes = crystal_planes(dimensions)

  ms = miller.build_set(crystal_symmetry=crystal.symmetry(space_group_symbol="P1",
    unit_cell=(a, b, c, 90, 90, 90)), anomalous_flag=False, d_min=0.5)
  B = ((1.0/a, 0.0, 0.0), (0.0, 1.0/b, 0.0), (0.0, 0.0, 1.0/c))
  U = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
  UB = multiply_matrices(U, B)
  #define direction of incoming beam
  k0 = (1.0, 0.0, 0.0)
  reflections_list = []
  miller_indices_list = []
  spot_vectors_list = []
  for p in range(0, 40, 1):
    phi = float(p/2)
    R = rotation_matrix(phi)
    reciprocal_lattice_points = transform_to_lab_frame(UB, R, ms)
    #calculate scattering vectors of all rlps
    scattering_vectors = calculate_scattering_vectors(k0, reciprocal_lattice_points)
    #determine which scattering vectors will cause a spot on the detector.
    xs, ys, miller_indices, spot_vectors = calculate_spot_positions(scattering_vectors, ms)
    phis = [phi]*len(xs)
    xyzs = zip(xs, ys, phis)
    reflections_list.extend(xyzs)
    miller_indices_list.extend(miller_indices)
    spot_vectors_list.extend(spot_vectors)
  print "calculated which spots would be detected during sweep - found %s instances" % (len(reflections_list))
  Ilist = []
  for idx, k1 in enumerate(spot_vectors_list):
    #rotate back to crystal frame first
    r = rotation_matrix(-1.0 * reflections_list[idx][2])
    k1prime = rotate_vector(r, k1)
    k0prime = rotate_vector(r, k0)
    #calculate the intensity based on paths through the crystal
    I = calculate_intensity_of_reflection(k0prime, k1prime, crystal_points, cryst_planes)
    Ilist.append(I)


  #apply LP correction
  from math import acos, sin
  for idx, I in enumerate(Ilist):
    k1 = spot_vectors_list[idx]
    twotheta = acos(((k1[0]*k0[0]) + (k1[1]*k0[1]) + (k1[2]*k0[2])) / normalise_vec(k1))
    LP = 1.0 / (sin(twotheta / 2.0) * sin(twotheta))
    Ilist[idx] = I * LP

  #return objects
  return reflections_list, miller_indices_list, spot_vectors_list, Ilist


refls, miller_idx, k1s, Is = simulate_dataset()

miller_indices = flex.miller_index(miller_idx)
intensities = flex.double(Is)
xyzpositions = flex.vec3_double(refls)

#now create a reflection table
reflections = flex.reflection_table()
reflections['miller_index'] = miller_indices
reflections['intensity'] = intensities
reflections['xyz'] = xyzpositions


print "(x,y,phi) pos, k1 vector, (h,k,l), Intensity"
for i in range(len(refls)):
  print (("(%.4f,%.4f,%.4f), (%.4f,%.4f,%.4f), %s, %s")
         % (refls[i][0], refls[i][1], refls[i][2],
            k1s[i][0], k1s[i][1], k1s[i][2],
            miller_idx[i], Is[i]))
