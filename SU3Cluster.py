#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from utils import *
# ------------------------------------------------------------------
# Run triality cluster simulation of SU(3) gauge theory

# Parse arguments: 4d lattice volume (will likely have nt=1 for now),
# canonical sector in terms of number of (three-quark) baryons,
# gauge coupling beta, number of sweeps to do, RNG seed
# and directory for output data
if len(sys.argv) < 10:
  print "Usage:", str(sys.argv[0]), "<nx> <ny> <nz> <nt> <#baryons>"
  print "                   <beta> <sweeps> <random_seed> <out_dir>"
  sys.exit(1)
nx = np.uint(sys.argv[1])
ny = np.uint(sys.argv[2])
nz = np.uint(sys.argv[3])
nt = np.uint(sys.argv[4])
vol = nx * ny * nz * nt
Ndim = 4                      # Number of dimension
Ndir = 2 * Ndim               # Number of directions (forward and backward)
NB = np.uint(sys.argv[5])     # Number of baryons
Nq = 3 * NB                   # Number of quarks
beta = float(sys.argv[6])
Nsweep = int(sys.argv[7])
seed = int(sys.argv[8])
outdir = sys.argv[9]
runtime = -time.time()

# In deconfined phase, a single cluster can fill the entire lattice
# So the recursive build_cluster routine in utils.py
# can easily exceed the default 10^3 recursion limit
# Need a little over vol in total
# since build_cluster calls follow_bond which then calls site_index
sys.setrecursionlimit(int(vol + 2))

# Compute and save these constant floats
TODO: TO BE UPDATED...
exp_mbe = np.exp(-beta)          # Also bond removal probability
add_prob = 1.0 - exp_mbe
split_prob = 3.0 * exp_mbe / (1.0 + 2.0 * exp_mbe)
merge_prob = add_prob / (1.0 + 2.0 * exp_mbe)
act_frac = -1.0 * beta / add_prob          # Factor for the action

# TODO: Utilities for loading bond configuration and occupation numbers...

# Create output directory if it doesn't exist already
if not os.path.isdir(outdir):
  print "Creating directory", outdir, "for output"
  os.makedirs(outdir)

# Save run parameters for posterity
PARAMS = open(outdir + '/params.txt', 'w')
print >> PARAMS, "python", ' '.join(sys.argv)

# Quick sanity check: Make sure all NB baryons can fit on the lattice
if NB > 2 * vol:
  print "ERROR: Cannot fit", NB, "baryons in",
  print nx, "x", ny, "x", nz, "lattice...",
  print "aborting"
  sys.exit(1)

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
t = np.empty(vol, dtype=np.uint)
for i, j, k, l in np.ndindex((nx, ny, nz, nt)):
  index = np.uint(i + nx * (j + ny * (k + nt * l)))
  x[index] = i
  y[index] = j
  z[index] = k
  t[index] = l

# Pack constant information into single variable for passing to subroutines
lattice = dict({'nx': nx, 'ny': ny, 'nz': nz, 'nt': nt,
                'Ndim': Ndim, 'Ndir': Ndir, 'vol': vol, 'prng': prng,
                'x': x, 'y': y, 'z': z, 't': t})

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
numBond = np.uint(0)
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
check_Nq(occupation, Nq)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Open files for output
ACCEPT = open(outdir + '/accept.csv', 'w')
print >> ACCEPT, "sweep,accept_mvB,accept_mvQ,accept_bond"
MAXCLUSTER = open(outdir + '/maxcluster.csv', 'w')
print >> MAXCLUSTER, "sweep,max_tot,max_rel"
AVECLUSTER = open(outdir + '/avecluster.csv', 'w')
print >> AVECLUSTER, "sweep,ave_tot,ave_rel"
NUMBONDS = open(outdir + '/numbonds.csv', 'w')
print >> NUMBONDS, "sweep,nb_tot,nb_rel"
ACTION = open(outdir + '/action.csv', 'w')
print >> ACTION, "sweep,action_tot,action_rel"

# Print starting state
count_clusters(root, numCluster, 0, MAXCLUSTER)

tot = float(vol) / float(numCluster)
rel = 1.0 / float(numCluster)
print >> AVECLUSTER, "0,%.8g,%.8g" % (tot, rel)

rel = float(numBond) / float(vol * Ndim)
print >> NUMBONDS, "0,%d,%.8g" % (numBond, rel)
if not beta == 0:
  tr = float(numBond) / add_prob
  print >> ACTION, "0,%.8g,%.8g" % (tr, tr / float(vol))
else:
  print >> ACTION, "0,0.0,0.0"

# Loop over sweeps, printing some basic data after each one
for sweep in range(1, Nsweep + 1):
  # Each sweep loops (randomly) over the lattice volume
  accept = [0.0, 0.0, 0.0]        # Initialize acceptance rate
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
        accept[0] += 1.0
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
          occupation[ran] -= np.uint(1)
          occupation[new] += np.uint(1)
          accept[1] += 1.0
    # --------------------------------------------------------------



    # --------------------------------------------------------------
    # Update step 3: Try to change bond
    ran = prng.randint(0, vol)
    ran_dir = prng.randint(0, Ndim)

    # Figure out the site on the other side of the bond
    neigh = follow_bond(ran, ran_dir, lattice)

    # If the bond is present, try to remove it
    if bond[ran][ran_dir]:
      bond[ran][ran_dir] = False      # Consequences to be checked...

      # Build cluster from ran, see if neigh is still in it
      # If no change in clusters, accept with probability exp_mbe
      ran_cluster = []
      connect = check_connect(bond, ran, ran_cluster, lattice, neigh)
      if connect > 0:         # No change in clusters
        if prng.uniform(0, 1) < exp_mbe:
          numBond -= np.uint(1)
          accept[2] += 1.0
        else:
          bond[ran][ran_dir] = True       # Reject!

      # If the cluster will be split we need to check the occupation numbers
      # connect<0 means that check_connect built the complete ran_cluster
      else:     # Cluster will be split
        ran_Nq = check_occupation(occupation, ran_cluster)
        if not np.mod(ran_Nq, 3) == 0:
          bond[ran][ran_dir] = True       # Reject!

        else:   # Accept with probability 3 * exp_mbe / (1 + 2 * exp_mbe)
                # (We already know that the other occupation number is fine)
          if prng.uniform(0, 1) < split_prob:
            accept[2] += 1.0
            numBond -= np.uint(1)
            numCluster += np.uint(1)

            # Build new cluster and reset roots
            neigh_cluster = []
            build_cluster(bond, neigh, neigh_cluster, lattice)
            for i in ran_cluster:
              root[i] = ran
            for i in neigh_cluster:
              root[i] = neigh

          else:   # The final reject!
            bond[ran][ran_dir] = True       # Reject!

    # If the bond is not present, try to add it
    else:
      ran_root = get_root(root, ran)
      neigh_root = get_root(root, neigh)
      # If both sites are already in the same cluster,
      # then add bond with probability (1 - exp_mbe)
      if ran_root == neigh_root:
        if prng.uniform(0, 1) < add_prob:
          bond[ran][ran_dir] = True
          numBond += np.uint(1)
          accept[2] += 1.0

      # Otherwise the addition decreases the number of clusters by one,
      # and so occurs with probability (1 - exp_mbe) / (1 + 2 * exp_mbe)
      else:
        if prng.uniform(0, 1) < merge_prob:
          bond[ran][ran_dir] = True
          numBond += np.uint(1)
          numCluster -= np.uint(1)
          root[neigh_root] = ran_root         # Merge clusters
          accept[2] += 1.0
    # --------------------------------------------------------------

  # Print some basic data after each sweep
  # (Can also run after each update if speed and output size aren't issues)
  # First print average acceptances for the sweep
  aB = accept[0] / float(vol)
  aQ = accept[1] / float(vol)
  aBond = accept[2] / float(vol)
  print >> ACCEPT, "%d,%.4g,%.4g,%.4g" % (sweep, aB, aQ, aBond)

  # Sanity check: make sure our total occupation number remains correct
  check_Nq(occupation, Nq)

  # Make sure our count of clusters remains correct
  # count_clusters prints size of largest cluster
  count_clusters(root, numCluster, sweep, MAXCLUSTER)

  # Print average cluster size, both absolute and as fraction of total volume
  tot = float(vol) / float(numCluster)
  rel = 1.0 / float(numCluster)
  print >> AVECLUSTER, "%d,%.8g,%.8g" % (sweep, tot, rel)

  # Make sure our count of bonds remains correct
  count_bonds(bond, numBond)

  # Print number of bonds, both absolute and as fraction of the total
  rel = float(numBond) / float(vol * Ndim)
  print >> NUMBONDS, "%d,%d,%.8g" % (sweep, numBond, rel)

  # Print action S = -beta * numBond / (1 - exp[-beta])
  # (Derived in Eq. 3.28 of Philippe Widmer's thesis)
  # Again, first total action then average divided by total volume
  # Note that numBond = 0 when beta = 0
  if not beta == 0:
    tr = float(numBond) * act_frac
    print >> ACTION, "%d,%.8g,%.8g" % (sweep, tr, tr / float(vol))
  else:
    print >> ACTION, "0,0.0,0.0"
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Clean up and close down
ACCEPT.close()
MAXCLUSTER.close()
AVECLUSTER.close()
NUMBONDS.close()
ACTION.close()

# TODO: Utilities for saving bond configuration and occupation numbers...

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
print >> PARAMS, "Runtime: %0.1f seconds" % runtime
PARAMS.close()
# ------------------------------------------------------------------
