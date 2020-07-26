# Version for Github
source('Tools/Step0_Globally_used_Functions_and_Datasets.R')
# Test the hypothesis that this effect is due to
# incomplete penetrance of exp by shRNA vs CRISPR-KO
# We are planning to achieve this by identifying if indeed these genes
# diff more ess in p53 WT vs Mut in crispr vs shRNA
########################################################################
##Find non-expressing genes
########################################################################
# From Expression profile of cell lines in AVANA::
# we identified genes with mean expression <0.1 (below is the list)
non_expressing_genes=readRDS('Data/non_expressing_genes.RDS')
########################################################################
##Association with TP53
########################################################################
COI=Reduce(intersect, list(colnames(avana),
                           colnames(Mut_CCLE),
                           colnames(achilles)))
GOI=Reduce(intersect, list(rownames(onTarget$avana), rownames(onTarget$mutations_matrix), rownames(onTarget$achilles)))
EffectSize_byMasterREgulator<-function(MR='TP53', 
                                       infunc_mat=avana[GOI, COI],
                                       Mut=Mut_CCLE[GOI, COI]){
  do.call(rbind, mclapply(1:nrow(infunc_mat),
                          function(x) c(wilcox.test(unlist(infunc_mat[x,])~unlist(Mut[MR,]))$p.value,
                                        median(unlist(infunc_mat[x,unlist(Mut[MR,])==0]), na.rm = T) -
                                          median(unlist(infunc_mat[x,unlist(Mut[MR,])==1]), na.rm = T) ),
                          mc.cores = detectCores() ) )
}
## For CRISPR effect SIze
p53_effect=EffectSize_byMasterREgulator('TP53'); rownames(p53_effect)=GOI
## For shRA effect SIze
p53_effect_shRNA=EffectSize_byMasterREgulator('TP53', infunc_mat = onTarget$achilles[GOI, COI])
rownames(p53_effect_shRNA)=GOI

# non-expressing genes are differential more essential in p53 WT vs mutant
# only in crispr screens specifially suggesting that incomplete phenotype 
# difference is not the key contributor of CDE genes identity.
wilcox.test(rowSubset(p53_effect, row_Names = non_expressing_genes)[,2], 
            rowSubset(p53_effect_shRNA, row_Names = non_expressing_genes)[,2],
            alternative='l')
