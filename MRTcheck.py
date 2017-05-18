#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from utils import follow_bond
# ------------------------------------------------------------------
# Check zero-density Potts model with Metropolis--Rosenbluth--Teller algorithm
# With weight exp[-S] = exp[gamma sum_<ij> \delta_{s_i, s_j}]
#   the acceptance probability is exp[gamma * (new# - cur#)]

# Parse arguments: 3d lattice volume,
# Potts coupling gamma, number of sweeps to do, RNG seed
# and directory for output data
if len(sys.argv) < 8:
  print "Usage:", str(sys.argv[0]), "<nx> <ny> <nz>"
  print "                     <gamma> <sweeps> <RNG seed> <output dir>"
  sys.exit(1)
nx = np.uint(sys.argv[1])
ny = np.uint(sys.argv[2])
nz = np.uint(sys.argv[3])
vol = nx * ny * nz
Ndim = 3                      # Number of dimension
Ndir = 2 * Ndim               # Number of directions (forward and backward)
Nstate = 3                    # Hard-code three-state Potts model
gamma = float(sys.argv[4])
Nsweep = np.uint(sys.argv[5])
seed = int(sys.argv[6])
outdir = sys.argv[7]
runtime = -time.time()

# TODO: Utilities for loading configuration...

# Create output directory if it doesn't exist already
if not os.path.isdir(outdir):
  print "Creating directory", outdir, "for output"
  os.makedirs(outdir)

# Save run parameters for posterity
PARAMS = open(outdir + '/params.csv', 'w')
print >> PARAMS, ' '.join(sys.argv)
PARAMS.close()

# Seed (Mersenne Twister) random number generator
# Use RandomState instead of (global) seed
# in case multiple independent streams may be needed in the future
prng = np.random.RandomState(seed)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Set up lattice
# First define arrays to store (x, y, z) indices of each site
x = np.empty(vol, dtype=np.uint)
y = np.empty(vol, dtype=np.uint)
z = np.empty(vol, dtype=np.uint)
for i, j, k in np.ndindex((nx, ny, nz)):
  index = np.uint(i + nx * (j + ny * k))
  x[index] = i
  y[index] = j
  z[index] = k

# Pack constant information into single variable for passing to subroutines
lattice = dict({'nx': nx, 'ny': ny, 'nz': nz, 'Ndim': Ndim, 'Ndir': Ndir,
                'vol': vol, 'prng': prng, 'x': x, 'y': y, 'z': z})

# Now for each site we need the following:
#   The state of the Potts 'spin'
# We start with randomly assigned values
config = np.empty(vol, dtype=np.uint)
for i in range(vol):
  config[i] = prng.randint(0, Nstate)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Open files for output
ACCEPT = open(outdir + '/accept.csv', 'w')
MAGNET = open(outdir + '/magnet.csv', 'w')
ACTION = open(outdir + '/action.csv', 'w')

# Print starting state
# Note S = -gamma sum_<ij> \delta_{s_i, s_j}]
tot_spin = -float(vol)      # Converts spins to (-1, 0, 1, ..., Nstate - 2)
tot_act = 0.0
for i in range(vol):
  tot_spin += config[i]
  for mu in range(Ndim):    # Only the forward neighbors
    neigh = config[follow_bond(i, mu, lattice)]
    if config[i] == neigh:
      tot_act -= gamma

# Print 'magnetization' and action,
# for each including both total and average over lattice volume
print >> MAGNET, tot_spin, tot_spin / float(vol)
print >> ACTION, tot_act, tot_act / float(vol)

# Loop over sweeps, printing some basic data after each one
for sweep in range(Nsweep):
  # Each sweep loops (randomly) over the lattice volume
  for i in range(vol):
    # Update: Try to change the state at the current site
    # The new state is allowed to be the current state
    ran = prng.randint(0, vol)
    cur = config[ran]
    new = prng.randint(0, Nstate)   # Proposed new state at site ran

    # Compute change in energy, if non-zero
    # With weight exp[-S] = exp[gamma sum_<ij> \delta_{s_i, s_j}]
    #   accept with probability exp[newE - curE]
    if new == cur:
      print >> ACCEPT, 1
    else:
      curE = 0.0
      newE = 0.0
      for mu in range(Ndir):
        neigh = config[follow_bond(ran, mu, lattice)]
        if cur == neigh:
          curE -= gamma
        if new == neigh:
          newE -= gamma

      if newE > curE:
        config[ran] = new
        print >> ACCEPT, 1
      elif prng.uniform(0, 1) < np.exp(newE - curE):
        config[ran] = new
        print >> ACCEPT, 1
      else:
        print >> ACCEPT, 0

  # Print some basic data after each sweep
  # (Can also run after each update if speed is not an issue)
  tot_spin = -float(vol)      # Converts spins to (-1, 0, 1, ..., Nstate - 2)
  tot_act = 0.0
  for i in range(vol):
    tot_spin += config[i]
    for mu in range(Ndim):    # Only the forward neighbors
      neigh = config[follow_bond(i, mu, lattice)]
      if config[i] == neigh:
        tot_act -= gamma

  # Print 'magnetization' and action,
  # for each including both total and average over lattice volume
  print >> MAGNET, tot_spin, tot_spin / float(vol)
  print >> ACTION, tot_act, tot_act / float(vol)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Clean up and close down
ACCEPT.close()
MAGNET.close()
ACTION.close()

# TODO: Utilities for saving configuration...

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
