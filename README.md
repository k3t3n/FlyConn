This repository contains code and data to reproduce the results described in **Circuit Analysis of the Drosophila Brain using Connectivity-based Neuronal Classification Reveals Organization of Key Information Processing Pathways**, by Ketan Mehta, Rebecca F. Goldin, and Giorgio A. Ascoli.

## Contents

**`data`** :
  - nat_flycircuit_neuronlist.rds - scraped from Natverse
  - A_str.rds - strength matrix (same as A_conn, but with corrected names)


**`generate_connectomes`** : Generate stochastic binary connectomes using FlyCircuit Data

**`para_d`** : Automatic picking of the dimenisonalty of embedding parameter using the scree plot

**`clustering`**
  - MBHAC_clustering
  - IVC Consensus Clustering Iterative Voting_clustering
  - final_clustering

**`stat_analysis`**
  - nmi: normalized mutual information to compare different classification schemes
  - ari: Adjusted Rand index 

**`labels`**
  - classification-labels
  - graph-theoretic labels

**`random_walk`** : absorbtion-driftiness

**`spad`**

**`morph_pvec`**
