This repository contains code and data to reproduce the results described in **Circuit Analysis of the Drosophila Brain using Connectivity-based Neuronal Classification Reveals Organization of Key Communication Pathways**, by Ketan Mehta, Rebecca F. Goldin, and Giorgio A. Ascoli.

## Contents

**`data`**
  - Flycircuit neuron list along with other metadata and descriptors
  - Strength connectome from FlyCircuit data

**`generate_connectomes`**
  - Generate stochastic binary connectomes using FlyCircuit Data
  - Automatic picking of the dimenisonalty of embedding parameter using the scree plot

**`clustering`**
  - MBHAC Clustering
  - IVC Consensus Clustering
  - Creating the final connectivity-based classes
  - Obtaining the edge probabilities by paramterizing as a SBM 

**`stat_analysis`**
  - NMI: normalized mutual information to compare different classification schemes
  - ARI: adjusted Rand index 

**`labels`**
  - Classification labels
  - Graph-theoretic labels

**`random_walk`** : Performing random walk on the connectivity-based circuit absorbtion-driftiness

**`spad`** : Spatial distribution of the neurons using FlyCricuit and NatVerse atlas

**`morph_pvec`** : Morphology and persistence vector classification
