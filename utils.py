#!/usr/bin/python
import numpy as np
# Some basic utilities for triality cluster code
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Sanity check: Make sure sum of all occupation numbers equals nq
def check_nq(occupation, nq):
  tot = 0
  tot = np.uint(tot)            # Set proper type
  for i, j, k in np.ndindex(occupation.shape):
#    print "occupation[%d][%d][%d] = %d" % (i, j, k, occupation[i][j][k])
    tot += occupation[i][j][k]
#  print tot
  if not tot == nq:
    print "ERROR: Counted", tot, "rather than", nq, "quarks...",
    print "aborting"
    sys.exit(1)
# ------------------------------------------------------------------
