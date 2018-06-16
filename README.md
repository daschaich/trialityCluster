# trialityCluster

Canonical cluster algorithms for the 3d three-state Potts model and heavy-dense SU(3) QCD

The triality cluster algorithmic approach is presented by:\
Andrei Alexandru, Georg Bergner, David Schaich and Urs Wenger, "Solution of the sign problem in the Potts model at fixed fermion number", [*Physical Review* **D97** (2018) 114503](http://doi.org/10.1103/PhysRevD.97.114503) [[arXiv:1712.07585](https://arxiv.org/abs/1712.07585)]

Here we're working in python to (hopefully) simplify development and experimentation.
Initial implementations are straightforward and simplistic, and may need significant refinements and optimizations to enable meaningful computations

Numerical Python ([NumPy](https://github.com/numpy/numpy)) is required by all the applications described below.\
All applications hard-code the number of dimensions as `Ndim`.  Changing the dimensionality should only require modifying the definitions of `Ndim` and the lattice volume `vol`, along with the input summarized below.  This has not yet been tested.

## Triality cluster algorithm for 3d three-state Potts model

`PottsCluster.py` is the main file for canonical three-dimensional three-state Potts model computations using the triality cluster algorithm, with additional utilities in `util.py`

This program takes eight input arguments:
```
python PottsCluster.py <nx> <ny> <nz> <baryons>
                       <gamma> <sweeps> <random_seed> <out_dir>
```

This sets up an `nx`x`ny`x`nz` lattice in the canonical sector with 3x`baryons` quarks.
The Potts model coupling is `gamma`.
`sweeps` update sweeps over the lattice are performed.
Each sweep does each of the following steps  `vol` times:
1. Choose a random site and try to move a full baryon from that site to its neighbor in a random direction.
2. Choose a random site and try to move a single quark from that site to its neighbor in a random direction (which must be in the same cluster).
3. Choose a random {site, direction} and try to change the bond (removing the bond if it's present, adding it if it's not), potentially changing the number of clusters.

The pseudorandom numbers are produced by NumPy's Mersenne Twister generator, initialized with the given `random_seed`.

Output is written to the following files in the output directory `out_dir` (which is created if it does not yet exist):
* `accept.csv` records the average acceptance for each of the three update steps listed above after each sweep
* `action.csv` records the (total and volume-averaged) Potts model action `NB/(1-exp(-gamma))` after each sweep, where NB is the total number of bonds present in the lattice
* `avecluster.csv` records the average size of each cluster after each sweep, in terms of both the number of sites and the fraction of the total volume
* `maxcluster.csv` records the size of the largest cluster after each sweep, in terms of both the number of sites and the fraction of the total volume
* `numbonds.csv` records the number of bonds in the lattice after each sweep, both the total number NB and the fraction of the maximum `Ndim`x`vol`
* `params.txt` records the input parameters and total runtime for reference

Existing files in the output directory are overwritten.\
The `csv` files are formatted as expected by [dygraphs](http://dygraphs.com) dynamical time-series plots.\
All five include a header line for such plots, and all but `accept.csv` also record the initial value before the first sweep.\
Therefore `accept.csv` should have `sweeps`+1 lines while the other four `csv` files should have `sweeps`+2 lines.

TODO:
* Improve performance on larger volumes, especially in the deconfined phase where the clusters can become very large
* Add routines to save and load configurations, appending to output files rather than overwriting them
* Reproduce results in arXiv:1712.07585 (will require additional update steps and/or reweighting)

## Local update algorithm for 3d three-state Potts model

`PottsMRT.py` is the main file for canonical three-dimensional three-state Potts model computations using the Metropolis--Rosenbluth--Teller (MRT) algorithm in the special case of zero baryon density where it does not suffer from a sign problem.  Currently this application uses only one routine from `util.py`.

Here the number of states `Nstate` is hard-coded

This program takes seven input arguments:
```
python PottsMRT.py <nx> <ny> <nz>
                   <gamma> <sweeps> <random_seed> <out_dir>
```

There are only two differences compared to the cluster application.\
First, we only work in the zero-quark canonical sector, to avoid the severe sign problem this algorithm would encounter at non-zero density.\
Second, each sweep over `vol` updates is much simpler.  Each update chooses a random site, sets its spin to a random value (which can be the same as it currently has), and runs the MRT accept/reject test.

As above, output is written to the following files in the output directory `out_dir` (which are created if they don't yet exist, overwritten if they do, and formatted as described above):
* `accept.csv` records the average acceptance for each sweep
* `action.csv` records the (total and volume-averaged) Potts model action (`gamma sum_<ij> delta_{s_i, s_j}`) after each sweep, where the sum is over all nearest-neighbor pairs of sites i and j
* `magnet.csv` records the (total and volume-averaged) magnetization defined by assigning the three Potts states the numerical values {-1, 0, 1}
* `params.txt` records the input parameters and total runtime for reference

TODO:
* Add routines to save and load configurations, appending to output files rather than overwriting them

## Triality cluster algorithm for SU(3) gauge theory

`SU3Cluster.py` is the main file for canonical SU(3) gauge theory computations using the triality cluster algorithm, with additional utilities still in `util.py`

The number of colors `Nc` is hard-coded

This program takes nine input arguments:
```
python SU3Cluster.py <nx> <ny> <nz> <nt> <baryons>
                     <beta> <sweeps> <random_seed> <out_dir>
```

This sets up an `nx`x`ny`x`nz`x`nt` lattice in the canonical sector with 3x`baryons` quarks.
The (inverse) gauge coupling is `beta`.

TODO: To be implemented and filled in...

## Local update algorithm for for SU(3) gauge theory

`SU3MRT.py` is the main file for canonical SU(3) gauge theory computations using the MRT algorithm in the special case of zero baryon density where it does not suffer from a sign problem.  Currently this application uses only one routine from `util.py`.

The number of colors `Nc` is hard-coded

This program takes eight input arguments:
```
python SU3MRT.py <nx> <ny> <nz> <nt>
                 <beta> <sweeps> <random_seed> <out_dir>
```

TODO: To be implemented and filled in...
* Check against pure-gauge over-relaxation algorithm in MILC
