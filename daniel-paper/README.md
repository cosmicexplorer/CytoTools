daniel-paper
============

the current version of the paper is in the doc in this repo

# experiments
## gbm patient
- [biaxial gating](https://irishlab.cytobank.org/cytobank/experiments/22103/illustrations/52353)
- [visne gating](https://irishlab.cytobank.org/cytobank/experiments/22324/illustrations/52355)
## pbmc patient
- [biaxial gating](https://irishlab.cytobank.org/cytobank/experiments/22228/illustrations/52354)
- [visne gating](https://irishlab.cytobank.org/cytobank/experiments/21509/illustrations/49979)

# methodology
- terminal populations chosen, noted in `manual_gating` column
    - this is essentially the "cluster id" but for manual pops
    - terminal populations **partition** the data
- for each cluster from algorithm, find population with most cells in cluster
    - find f-measure *from* cluster *to* that population
    - take mean f-measure for all clusters to get clustering quality
    - rest is in paper

# linx
- https://en.wikipedia.org/wiki/Ensemble_learning
