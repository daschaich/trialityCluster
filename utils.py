#!/usr/bin/python
import numpy as np
# Some basic utilities for triality cluster code
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Sanity check: Make sure sum of all occupation numbers equals nq
def check_nq(occupation, lattice):
  tot = 0
  tot = np.uint(tot)            # Set proper type
  for i, j, k in np.ndindex(occupation.shape):
#    print "occupation[%d][%d][%d] = %d" % (i, j, k, occupation[i][j][k])
    tot += occupation[i][j][k]
#  print tot
  if not tot == lattice['nq']:
    print "ERROR: Counted", tot, "rather than", nq, "quarks...",
    print "aborting"
    sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Randomly choose neighboring site
# Use MILC ordering conventions for directions
def get_neighbor(x, y, z, lattice):
  new_x = x
  new_y = y
  new_z = z
  ran_dir = lattice['prng'].randint(0, lattice['Ndir'])
  if ran_dir == 0:
    new_x += 1
  elif ran_dir == 1:
    new_y += 1
  elif ran_dir == 2:
    new_z += 1
  elif ran_dir == 3:
    new_z -= 1
  elif ran_dir == 4:
    new_y -= 1
  elif ran_dir == 5:
    new_x -= 1

  return new_x % lattice['nx'], new_y % lattice['ny'], new_z % lattice['nz']
# ------------------------------------------------------------------
