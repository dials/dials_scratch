
from scitbx import matrix

rlp = matrix.col((-0.549237619333,0.529065437127,0.33645893324))
s0 = matrix.col((-0.000141687134493,0.0,-1.02432777509))

#s0 = matrix.col((0, 0, 1))

#rlp = matrix.col((0, -1, 0))

s1 = s0 + rlp
e1 = s1.cross(s0).normalize()
e2 = -e1.cross(rlp.normalize())
e3 = rlp.normalize()

# print e1.dot(e2)
# print e1.dot(e3)
# print e2.dot(e3)

E = matrix.sqr((e1[0], e2[0], e3[0],
                e1[1], e2[1], e3[1],
                e1[2], e2[2], e3[2]))


L = matrix.sqr((
  0.1, 0, 0,
  0, 0.1, 0,
  0, 0, 0.1))

Sigma = L*L.transpose()

print E.transpose() * Sigma * E

print E.is_r3_rotation_matrix()
print "e1: ", tuple("%.2f" % e for e in e1)
print "e2: ", tuple("%.2f" % e for e in e2)
print "e3: ", tuple("%.2f" % e for e in e3)


E = E.transpose() * Sigma * E

print Sigma.inverse()




# print E.transpose() * L * E
# print E
# print tuple(Sigma)

# print E.transpose() * Sigma * E

# print E.inverse().is_r3_rotation_matrix()
# print ""

# F = E.transpose() * L * E
# print F
# print E
# print E * E.transpose()

