# Version for Github
source('Step0_Globally_used_Functions_and_Datasets.R')

# Compute the emperical P-values for the three master regulators

Master_Regulators_Cand=c('TP53','VHL','KRAS')
Empipical_Significance=sapply(Master_Regulators_Cand, 
                 function(x) sapply(1:100, function(dummy_var) Randomized_Testing_CRISPR_damage_bias(x) ))

