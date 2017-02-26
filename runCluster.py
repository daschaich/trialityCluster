#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from utils import check_nq
# ------------------------------------------------------------------
# Run triality cluster simulation of Potts model for canonical heavy-dense QCD

# Parse arguments: 3d lattice volume, Potts coupling gamma,
# canonical sector in terms of number of (three-quark) baryons,
# number of sweeps to do, directory for output data
if len(sys.argv) < 9:
  print "Usage:", str(sys.argv[0]), "<nx> <ny> <nz> <#baryons>"
  print "                     <gamma> <sweeps> <RNG seed> <output dir>"
  sys.exit(1)
nx = int(sys.argv[1])
ny = int(sys.argv[2])
nz = int(sys.argv[3])
vol = nx * ny * nz
NB = int(sys.argv[4])       # Number of baryons
nq = 3 * NB                 # Number of quarks
gamma = float(sys.argv[5])
Nsweep = int(sys.argv[6])
seed = int(sys.argv[7])
outdir = sys.argv[8]

# Create output directory if it doesn't exist already
if not os.path.isdir(outdir):
  print "Creating directory", outdir, "for output"
  os.makedirs(outdir)

# Quick sanity check: Make sure all NB baryons can fit on the lattice
if NB > 2 * vol:
  print "ERROR: Cannot fit", NB, "baryons in",
  print nx, "x", ny, "x", nz, "lattice...",
  print "aborting"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Set up lattice
# Basic structure: Site as list of three coordinates
site = [('x', np.int), ('y', np.int), ('z', np.int)]

# For each site we need the following:
#   An occupation number (counting quarks, not baryons)
#   Six booleans to tell whether or not bonds are present
#   A pointer to the site at the root of its cluster
# If NB > volume, start with full lattice and remove NB baryons
# Otherwise start with empty lattice and add NB baryons
if NB > vol:
  occupation = np.full((nx, ny, nz), 6, dtype=np.uint)
else:
  occupation = np.zeros((nx, ny, nz), dtype=np.uint)
bond = np.zeros((nx, ny, nz, 6), dtype=bool)        # All False
root = np.empty((nx, ny, nz), dtype=site)

for i, j, k in np.ndindex((nx, ny, nz)):
  root[i][j][k] = (i, j, k)

# Seed random number generator
# Use RandomState instead of (global) seed
# in case multiple independent streams may be needed in the future
prng = np.random.RandomState(seed)

# Initialize quark configuration as described above
# Either add or remove baryons at randomly chosen sites
if NB > vol:                # Remove (2vol - NB) baryons from full lattice
  for i in range(2 * vol - NB):
    success = False         # Keep trying until success!
    while not success:
      ran_x = prng.randint(0, nx)
      ran_y = prng.randint(0, ny)
      ran_z = prng.randint(0, nz)
      if occupation[ran_x][ran_y][ran_z] > 2:
        occupation[ran_x][ran_y][ran_z] -= 3
        success = True
else:                       # Add NB baryons to empty lattice
  for i in range(NB):
    success = False
    while not success:
      ran_x = prng.randint(0, nx)
      ran_y = prng.randint(0, ny)
      ran_z = prng.randint(0, nz)
      if occupation[ran_x][ran_y][ran_z] < 4:
        occupation[ran_x][ran_y][ran_z] += 3
        success = True

# Check for successful layout
check_nq(occupation, nq)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Sweeps!
#for traj_done in range(Nsweep):
  # TODO: Update
  # TODO: Measure
# ------------------------------------------------------------------
