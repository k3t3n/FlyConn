## ----------------------------------------------------------------------------
##--Preamble
rm(list=ls())

## ----------------------------------------------------------------------------
##--Sort according to neuron names
neuron_names <- readRDS('/media/WDHDD/clustering/exp/ARGO/P05/MBHAC/data/neuron_names.rds')
morpho23 <- readRDS('/media/WDHDD/clustering/morphology/morpho23feat.rds')

  ##--creates a single empty class vector
  class_vec = rep(0,length(neuron_names))
  ## loop through all neurons
  m_index = NULL
  for(i in 1:length(class_vec))
  {
    m_index = which(morpho23$Neuron_name==neuron_names[i])
    if(length(m_index)==0)
    {
      class_vec[i] = 0
    } else {
      #class_vec[i] = morpho23$Neuron_name[m_index]
      class_vec[i] = m_index
    }
  }

saveRDS(class_vec,file='/media/WDHDD/clustering/morphology/class_vec_morph.rds')
