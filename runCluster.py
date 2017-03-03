#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from utils import *
# ------------------------------------------------------------------
# Run triality cluster simulation of Potts model for canonical heavy-dense QCD

# Parse arguments: 3d lattice volume, Potts coupling gamma,
# canonical sector in terms of number of (three-quark) baryons,
# number of sweeps to do, directory for output data
if len(sys.argv) < 9:
  print "Usage:", str(sys.argv[0]), "<nx> <ny> <nz> <#baryons>"
  print "                     <gamma> <sweeps> <RNG seed> <output dir>"
  sys.exit(1)
nx = np.uint(sys.argv[1])
ny = np.uint(sys.argv[2])
nz = np.uint(sys.argv[3])
vol = nx * ny * nz
Ndim = 3                      # Number of dimension
Ndir = 6                      # Number of directions (forward and backward)
NB = np.uint(sys.argv[4])     # Number of baryons
nq = 3 * NB                   # Number of quarks
gamma = float(sys.argv[5])
Nsweep = np.uint(sys.argv[6])
seed = int(sys.argv[7])
outdir = sys.argv[8]
runtime = -time.time()

# Create output directory if it doesn't exist already
if not os.path.isdir(outdir):
  print "Creating directory", outdir, "for output"
  os.makedirs(outdir)

# Save run parameters for posterity
PARAMS = open(outdir + '/params.csv', 'w')
print >> PARAMS, ' '.join(sys.argv)
PARAMS.close()

# Quick sanity check: Make sure all NB baryons can fit on the lattice
if NB > 2 * vol:
  print "ERROR: Cannot fit", NB, "baryons in",
  print nx, "x", ny, "x", nz, "lattice...",
  print "aborting"
  sys.exit(1)

# Seed random number generator
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
                'vol': vol, 'nq': nq, 'prng': prng, 'x': x, 'y': y, 'z': z})

# Now for each site we need the following:
#   An occupation number (counting quarks, not baryons)
#   Ndim booleans to tell whether or not bonds are present
#   A pointer to the site at the root of its cluster
# We start with vol single-site clusters
# This requires {0, 3, 6} quarks at each site
# If NB > vol, start with full lattice and remove (2vol - NB) baryons
# Otherwise start with empty lattice and add NB baryons
if NB > vol:
  occupation = np.full(vol, 6, dtype=np.uint)
else:
  occupation = np.zeros(vol, dtype=np.uint)
bond = np.zeros((vol, Ndim), dtype=bool)        # All False
root = np.arange(vol, dtype=np.uint)            # root[i] = i

# Some gross features of configuration: Average cluster size,
# size of largest cluster, total numbers of clusters and bonds
aveCluster = 1.0 / float(vol)
maxCluster = 1.0 / float(vol)
numBonds = np.uint(0)
numCluster = vol

# Initialize quark configuration as described above
# Either add or remove baryons at randomly chosen sites
if NB > vol:                # Remove (2vol - NB) baryons from full lattice
  for i in range(2 * vol - NB):
    success = False         # Keep trying until success!
    while not success:
      ran = prng.randint(0, vol)
      if occupation[ran] > 2:
        occupation[ran] -= 3
        success = True
else:                       # Add NB baryons to empty lattice
  for i in range(NB):
    success = False
    while not success:
      ran = prng.randint(0, vol)
      if occupation[ran] < 4:
        occupation[ran] += 3
        success = True

# Check that layout was successful
check_nq(occupation, nq)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Open files for output
ACCEPT = open(outdir + '/accept.csv', 'w')
MAXCLUSTER = open(outdir + '/maxcluster.csv', 'w')
AVECLUSTER = open(outdir + '/avecluster.csv', 'w')
NUMBONDS = open(outdir + '/numbonds.csv', 'w')

# Loop over sweeps, printing some basic data after each one
for sweep in range(Nsweep):
  # Each sweep loops (randomly) over the lattice volume
  for i in range(vol):
    # --------------------------------------------------------------
    # Update step 1: Try to move full baryon to neighboring site
    # Check that we have a baryon to move
    ran = prng.randint(0, vol)
    if occupation[ran] > 2:
      # Choose random neighbor and see if it can accept the baryon
      new = get_neighbor(ran, lattice)
      if occupation[new] < 4:
        occupation[ran] -= 3
        occupation[new] += 3
        print >> ACCEPT, 1,
      else:
        print >> ACCEPT, 0,
    else:
      print >> ACCEPT, 0,
    # --------------------------------------------------------------



    # --------------------------------------------------------------
    # Update step 2: Try to move quark within cluster
    # Check that we have a quark to move
    ran = prng.randint(0, vol)
    if occupation[ran] > 0:
      # Choose random neighbor and see if it can accept the quark
      new = get_neighbor(ran, lattice)
      if occupation[new] < 6:
        # See whether or not both sites are in the same cluster
        if get_root(root, ran) == get_root(root, new):
          occupation[ran] -= 1
          occupation[new] += 1
          print >> ACCEPT, 1,
        else:
          print >> ACCEPT, 0,
      else:
        print >> ACCEPT, 0,

    else:
      print >> ACCEPT, 0,
    # --------------------------------------------------------------



    # --------------------------------------------------------------
    # Update step 3: Try to change bond
    ran = prng.randint(0, vol)
    ran_dir = prng.randint(0, Ndim)

    # If bond is present, try to remove it
    if bond[ran][ran_dir]:
      #TODO
      pass

    # If bond is not present, try to add it
    else:
      #TODO
      pass

    print >> ACCEPT, 0

  # Print some basic data after each sweep
  # (Can also run after each update if speed is not an issue)
  # Sanity check: make sure our total occupation number remains correct
  check_nq(occupation, nq)

  # Make sure our count of clusters remains correct
  # count_clusters prints size of largest cluster
  count_clusters(root, numCluster, MAXCLUSTER)

  # Print average cluster size as fraction of total volume
  print >> AVECLUSTER, 1.0 / float(numCluster)

  # Make sure our count of bonds remains correct
  count_bonds(bond, numBonds)
  print >> NUMBONDS, numBonds
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Clean up and close down
ACCEPT.close()
MAXCLUSTER.close()
AVECLUSTER.close()
NUMBONDS.close()

# TODO: Utilities for saving bond configuration and occupation numbers...

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
