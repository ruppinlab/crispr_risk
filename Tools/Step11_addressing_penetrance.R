# Version for Github
setwd('/Users/sinhas8/crispr_risk-master/')
source('Tools/Step0_Globally_used_Functions_and_Datasets.R')
source('/Users/sinhas8/my')
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
GOI=Reduce(intersect, list(rownames(avana),
                           rownames(Mut_CCLE),
                           rownames(achilles)))
avana_range01=apply(avana, 1, range01)
avana_range01=t(avana_range01)
achilles_range01=apply(achilles, 1, range01)
achilles_range01=t(achilles_range01)

EffectSize_byMasterREgulator<-function(MR='TP53', 
                                       infunc_mat=avana_range01[GOI, COI],
                                       Mut=Mut_CCLE[GOI, COI]){
  do.call(rbind, mclapply(1:nrow(infunc_mat),
                          function(x) c(wilcox.test(unlist(infunc_mat[x,])~unlist(Mut[MR,]==0))$p.value,
                                        median(unlist(infunc_mat[x,unlist(Mut[MR,])==0]), na.rm = T) -
                                          median(unlist(infunc_mat[x,unlist(Mut[MR,])!=0]), na.rm = T) ),
                          mc.cores = detectCores() ) )
}

## For CRISPR effect SIze
p53_effect=EffectSize_byMasterREgulator('TP53', infunc_mat = avana_range01[GOI, COI]); rownames(p53_effect)=GOI
## For shRA effect SIze
p53_effect_shRNA=EffectSize_byMasterREgulator('TP53', infunc_mat = achilles_range01[GOI, COI])
rownames(p53_effect_shRNA)=GOI

# non-expressing genes are differential more essential in p53 WT vs mutant
# only in crispr screens specifially suggesting that incomplete phenotype 
# difference is not the key contributor of CDE genes identity.
p53_effect_df=data.frame(p=p53_effect[,1],
           eff_size=p53_effect[,2],
           scaled_eff_size=range01(p53_effect[,2]) )
p53_effect_shRNA_df=data.frame(p=p53_effect_shRNA[,1],
                         eff_size=p53_effect_shRNA[,2],
                         scaled_eff_size=range01(p53_effect_shRNA[,2]) )

wilcox.test(rowSubset(p53_effect_df, row_Names = non_expressing_genes)[,3], 
            rowSubset(p53_effect_shRNA_df, row_Names = non_expressing_genes)[,3],
            alternative='g')
