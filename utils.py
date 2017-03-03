#!/usr/bin/python
import sys
import numpy as np
# Some basic utilities for triality cluster code
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Simple helper utility to convert from (x, y, z) to single unsigned int
# !!! Currently assuming that x<nx, etc.  Would be safer to check this
def site_index(x, y, z, lattice):
  return np.uint(x + lattice['nx'] * (y + lattice['ny'] * z))
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Randomly choose neighboring site
# Use MILC ordering conventions for directions
def get_neighbor(site, lattice):
  new_x = lattice['x'][site]
  new_y = lattice['y'][site]
  new_z = lattice['z'][site]
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

  # Keep new site within proper range of lattice volume
  new_x = np.uint(np.mod(new_x, lattice['nx']))
  new_y = np.uint(np.mod(new_y, lattice['ny']))
  new_z = np.uint(np.mod(new_z, lattice['nz']))
  return np.uint(site_index(new_x, new_y, new_z, lattice))
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Sanity check: Make sure sum of all occupation numbers equals nq
def check_nq(occupation, nq):
  tot = np.uint(0)        # Set proper type
  for i in range(len(occupation)):
#    print "occupation[%d] = %d" % (i, occupation[i])
    tot += occupation[i]
#  print tot
  if not tot == nq:
    print "ERROR: Counted", tot, "rather than", nq, "quarks... aborting"
    sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Simple helper utility to get root of given site
def get_root(root, site):
  # Follow pointers from each site to the root of its cluster
  ptr = root[site]
  while not root[ptr] == ptr:
    new = root[ptr]     # Let's be cautious about overwriting ptr...
    ptr = new
  return ptr
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Determine sizes of all clusters, printing size of largest
def count_clusters(root, numClusters, MAXCLUSTER):
  # Can have up to len(root) clusters -- maybe more than we need
  clusters = np.zeros(len(root), dtype=np.uint)
  for i in range(len(root)):
    ptr = get_root(root, i)

    # Increment size of corresponding cluster
    clusters[ptr] += 1

  # Check that all sites are accounted for
  if not np.sum(clusters) == len(root):
    print "ERROR: Counted %d sites in clusters rather than %d... aborting" \
          % (np.sum(clusters), len(root))
    sys.exit(1)

  # Count total number of clusters and check against numClusters
  tot = np.uint(0)
  for i in range(len(root)):
    if clusters[i] > 0:
      tot += 1

  if not tot == numClusters:
    print "ERROR: Counted", tot, "rather than", numClusters, "clusters...",
    print "aborting"
    sys.exit(1)

  # Print largest cluster size as fraction of total volume
  print >> MAXCLUSTER, np.amax(clusters) / float(len(root))
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Count total number of bonds
def count_bonds(bond, numBonds):
  tot = np.uint(0)        # Set proper type
  for i, mu in np.ndindex(bond.shape):
    if bond[i][mu]:
      tot += 1
  if not tot == numBonds:
    print "ERROR: Counted", tot, "rather than", numBonds, "bonds...",
    print "aborting"
    sys.exit(1)
# ------------------------------------------------------------------
