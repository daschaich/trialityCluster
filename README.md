# trialityCluster

Canonical cluster algorithms for the 3d three-state Potts model and heavy-dense SU(3) QCD

The triality cluster algorithmic approach is presented by:
Andrei Alexandru, Georg Bergner, David Schaich and Urs Wenger, *Solution of the sign problem in the Potts model at fixed fermion number*, [arXiv:1712.07585](https://arxiv.org/abs/1712.07585)

Here we're working in python to (hopefully) simplify development and experimentation

Performance optimization is not an immediate priority.
Initial implementations are straightforward and simplistic, and may need significant refinements to enable meaningful computations.

## Triality cluster algorithm for 3d three-state Potts model

`PottsCluster.py` is the main file for canonical three-dimensional three-state Potts model computations using the triality cluster algorithm, with additional utilities in `util.py`.
Numerical Python ([NumPy](https://github.com/numpy/numpy)) is required.

Because we use the bond formulation reviewed in [arXiv:1712.07585](https://arxiv.org/abs/1712.07585), the number of states appears only in some of the update acceptance probabilities; for now these are hard-coded for the three-state model

The code is mostly written in terms of an arbitrary number of dimensions `Ndim`.
To change the dimensionality, the user should only need to modify the definitions of `Ndim` and the lattice volume `vol` along with the input summarized below.

This program takes eight input arguments:
```bash
python PottsCluster.py <nx> <ny> <nz> <baryons>
                       <gamma> <sweeps> <random_seed> <out_dir>
```

This sets up an `nx`x`ny`x`nz` lattice in the canonical sector with 3x`baryons` quarks.
The Potts model coupling is `gamma`.
`sweeps` update sweeps over the lattice are performed.
Each sweep does each of the following steps  `vol` times:
1. Choose a random site and try to move a full baryon from that site to its neighbor in a random direction.
2. Choose a random site and try to move a single quark from that site to its neighbor in a random direction (which must be in the same cluster).
3. Choose a random site and direction and try to change the bond (removing the bond if it's present, adding it if it's not), potentially changing the number of clusters.

The pseudorandom numbers are generated using NumPy's Mersenne Twister generator, initialized with the given `seed`.

Output is written to the following standard files in the output directory `out_dir` (which is created if it does not yet exist):
* `accept.csv` records acceptance (0 vs. 1) for each of the three update steps listed above on each of its `sweeps`x`vol` lines
* `action.csv` records the (total and volume-averaged) Potts model action NB/(1-exp(-gamma)) after each sweep, where NB is the total number of bonds present in the lattice
* `avecluster.csv` records the average size of each cluster after each sweep, in terms of both the number of sites and the fraction of the total volume
* `maxcluster.csv` records the size of the largest cluster after each sweep, in terms of both the number of sites and the fraction of the total volume
* `numbonds.csv` records the number of bonds in the lattice after each sweep, both the total number NB and the fraction of the maximum `Ndim`x`vol`
* `params.csv` records the input parameters and total runtime for reference

Existing files in the output directory are overwritten.
All of `accept.csv`, `action.csv`, `avecluster.csv`, `maxcluster.csv` and `numbonds.csv` also record the initial value of the observables listed above, and therefore should have `sweeps`+1 lines.

TODO:
* Add routines to save and load configurations, appending to output files rather than overwriting them
* Check results against local update algorithm

## Local update algorithm for 3d three-state Potts model

`PottsMRT.py` is the main file for canonical three-dimensional three-state Potts model computations using the Metropolis--Rosenbluth--Teller (MRT) algorithm in the special case of zero baryon density where it does not suffer from a sign problem.

TODO: Details go here...

## Triality cluster algorithm for SU(3) variables

TODO: Work in progress...
