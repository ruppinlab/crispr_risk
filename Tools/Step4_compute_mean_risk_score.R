# Version for Github
source('Step0_Globally_used_Functions_and_Datasets.R')

#### computing the mean risk score for each master regulator by taking the median differential essentiality value between mutant and WT across CDE+ genes.

Mean_Risk_score<-function(GeneList, MR){
  median(apply(avana[GeneList, ], 1, function(x) ((median(scale(x)[as.logical(Mut_CCLE[MR,])], na.rm = T)
    - median(scale(x)[!as.logical(Mut_CCLE[MR,])], na.rm = T))) ), na.rm = T)
}

CDE_for_volg=readRDS('../Data/DE_posNneg_genes_for_VolgGenes.RDS')

Mean_Score_for_MasterRegulators=sapply(1:length(CDE_for_volg), function(x)
  err_handle(Mean_Risk_score(as.character(CDE_for_volg[[x]][[2]]$GeneName), MR=names(CDE_for_volg)[x])) )

names(Mean_Score_for_MasterRegulators)=names(CDE_for_volg)

Mean_Score_for_MasterRegulators[c('VHL', 'KRAS', 'TP53')]
# a vector comprising median risk of master regulators mutant selection in case their CDE+ genes are CRISPR KO.

