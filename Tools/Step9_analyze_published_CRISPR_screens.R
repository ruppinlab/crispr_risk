# Version for Github
setwd('/Users/sinhas8/crispr_risk-master/')
source('Tools/Step0_Globally_used_Functions_and_Datasets.R')

# Centralized pipeline for All screenings
########################################################################
# Load CDE positive genes identified from avana screens
########################################################################
CDE=read.csv('../Data/CDEs_of_ThreeMaster_Regulators.csv')
CDE=lapply(CDE, function(x) as.character(x[x!='']) )
# ***Publicly available Data***
################################################################
##KRAS DLD1 - shRNA shRNA
################################################################
#shRNA screenings in DLD1 isogenic KRAS mutated and WT pair. 
# from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2768667/
# Ji Luo et al.
df_DLD1_shRNA=read.csv('../Data/Isogenic_Screenings/KRAS_isogenic_shRNA_screening.csv')
df_DLD1_shRNA=df_DLD1_shRNA[!is.na(match(df_DLD1_shRNA$Gene_Symbol, rownames(avana))),]
df_DLD1_shRNA$GeneType=NA
df_DLD1_shRNA$GeneType[!is.na(match(df_DLD1_shRNA$Gene_Symbol, CDE$KRAS.CDE.Pos))]='KRAS CDE Pos'
df_DLD1_shRNA$GeneType[!is.na(match(df_DLD1_shRNA$Gene_Symbol, CDE$KRAS.CDE.Neg))]='KRAS CDE Neg'
df_DLD1_shRNA$GeneType[!(!is.na(match(df_DLD1_shRNA$Gene_Symbol, CDE$KRAS.CDE.Pos)) |
                     !is.na(match(df_DLD1_shRNA$Gene_Symbol, CDE$KRAS.CDE.Neg))) ]='KRAS CDE Neutral'
df_DLD1_shRNA$Rank_WT=rank(df_DLD1_shRNA$KRAS.WT.log2.Mean)
df_DLD1_shRNA$Rank_Mut=rank(df_DLD1_shRNA$KRAS.Mut.log2.mean)
test_differential_essentiality_DLD1_shRNA<-function(cde_type,
                              which_tail='t',
                              Instead_rankMode=F){
  if(Instead_rankMode==F){
    FC_inMutant=df_DLD1_shRNA$KRAS.Mut.log2.mean[df_DLD1_shRNA$GeneType==cde_type]
    FC_inWT=df_DLD1_shRNA$KRAS.WT.log2.Mean[df_DLD1_shRNA$GeneType==cde_type]
  }else {
    FC_inMutant=df_DLD1_shRNA$Rank_Mut[df_DLD1_shRNA$GeneType==cde_type]
    FC_inWT=df_DLD1_shRNA$Rank_WT[df_DLD1_shRNA$GeneType==cde_type]
  }
  wilcox.test(FC_inMutant, FC_inWT, paired = T, alternative = which_tail)$p.value
}

# *Results* : # Testing CDE+/- essentiality direction in shRNA screening
# Using direct FC rank score
test_differential_essentiality_DLD1_shRNA(cde_type="KRAS CDE Neg",
                                    which_tail='l',
                                    Instead_rankMode=T)
test_differential_essentiality_DLD1_shRNA(cde_type="KRAS CDE Pos",
                                    which_tail='g',
                                    Instead_rankMode=T)
# Using direct FC score
test_differential_essentiality_DLD1_shRNA(cde_type="KRAS CDE Neg",
                                          which_tail='l',
                                          Instead_rankMode=F)
test_differential_essentiality_DLD1_shRNA(cde_type="KRAS CDE Pos",
                                          which_tail='g',
                                          Instead_rankMode=F)

################################################################
##KRAS CRSIPR- DLD1
################################################################
# CRISPR-Cas9 screenings in DLD1 isogenic KRAS mutated and WT pair. 
# from https://www.cell.com/cell-reports/pdfExtended/S2211-1247(17)30892-6
# Martine et al 2017

# Load Martin et al Screens for DLD1 cell lines 
df_DLD1_crispr=read.csv('../Data/Isogenic_Screenings/KRAS_genomeWide_CRISPRcas9_screen_DLD1_Martin_.csv')
df_DLD1_crispr=df_DLD1_crispr[!is.na(match(df_DLD1_crispr$Gene, rownames(avana))),]
df_DLD1_crispr$GeneType=NA
df_DLD1_crispr$GeneType[!is.na(match(df_DLD1_crispr$Gene, CDE$KRAS.CDE.Pos))]='KRAS CDE Pos'
df_DLD1_crispr$GeneType[!is.na(match(df_DLD1_crispr$Gene, CDE$KRAS.CDE.Neg))]='KRAS CDE Neg'
df_DLD1_crispr$GeneType[!(!is.na(match(df_DLD1_crispr$Gene, CDE$KRAS.CDE.Pos)) |
                     !is.na(match(df_DLD1_crispr$Gene, CDE$KRAS.CDE.Neg))) ]='KRAS CDE Neutral'

df_DLD1_crispr$Rank_WT=rank(df_DLD1_crispr$WT.log2.fold.change)
df_DLD1_crispr$Rank_Mut=rank(df_DLD1_crispr$KRas.mutant.log2.fold.change)

##Test-Try1:: Method 1
test_direction_krascrispr_DLD1<-function(cde_type=levels(factor(df_DLD1_crispr_kras$GeneType))[1],
                                    which_tail='t', 
                                    compare_FC.ranks_mode=T,
                                    df_DLD1_crispr_kras=df_DLD1_crispr){
  if(compare_FC.ranks_mode==F){
    FC_inMutant=df_DLD1_crispr_kras$WT.log2.fold.change[df_DLD1_crispr_kras$GeneType==cde_type]
    FC_inWT=df_DLD1_crispr_kras$KRas.mutant.log2.fold.change[df_DLD1_crispr_kras$GeneType==cde_type]
  }else {
    FC_inMutant=df_DLD1_crispr_kras$Rank_Mut[df_DLD1_crispr_kras$GeneType==cde_type]
    FC_inWT=df_DLD1_crispr_kras$Rank_WT[df_DLD1_crispr_kras$GeneType==cde_type]
  }
  wilcox.test(FC_inMutant, FC_inWT, paired = T, alternative = which_tail)$p.value
}
# *Results* : Testing CDE+/- essentiality direction in CRISPR DLD1 screening
# Using direct FC score
test_direction_krascrispr_DLD1(cde_type="KRAS CDE Neg",
                          which_tail='l',
                          compare_FC.ranks_mode=F,
                          df_DLD1_crispr_kras=df_DLD1_crispr)
test_direction_krascrispr_DLD1(cde_type="KRAS CDE Pos",
                          which_tail='g',
                          compare_FC.ranks_mode=F,
                          df_DLD1_crispr_kras=df_DLD1_crispr)
# Using direct FC rank score
test_direction_krascrispr_DLD1(cde_type="KRAS CDE Neg",
                          which_tail='l',
                          compare_FC.ranks_mode=T,
                          df_DLD1_crispr_kras=df_DLD1_crispr)
test_direction_krascrispr_DLD1(cde_type="KRAS CDE Pos",
                          which_tail='g',
                          compare_FC.ranks_mode=T,
                          df_DLD1_crispr_kras=df_DLD1_crispr)

################################################################
##KRAS CRSIPR - HCT116
################################################################
df_HCT116_crispr=read.csv('../Data/Isogenic_Screenings/KRAS_Martin_HCT116.csv')
df_HCT116_crispr=df_HCT116_crispr[!is.na(match(df_HCT116_crispr$Gene, rownames(avana))),]
df_HCT116_crispr$GeneType=NA
df_HCT116_crispr$GeneType[!is.na(match(df_HCT116_crispr$Gene, CDE$KRAS.CDE.Pos))]='KRAS CDE Pos'
df_HCT116_crispr$GeneType[!is.na(match(df_HCT116_crispr$Gene, CDE$KRAS.CDE.Neg))]='KRAS CDE Neg'
df_HCT116_crispr$GeneType[!(!is.na(match(df_HCT116_crispr$Gene, CDE$KRAS.CDE.Pos)) |
                            !is.na(match(df_HCT116_crispr$Gene, CDE$KRAS.CDE.Neg))) ]='KRAS CDE Neutral'

df_HCT116_crispr$Rank_WT=rank(df_HCT116_crispr$WT.log2.fold.change)
df_HCT116_crispr$Rank_Mut=rank(df_HCT116_crispr$KRas.mutant.log2.fold.change)

##Test-Try1:: Method 1
test_direction_krascrispr_DLD1<-function(cde_type=levels(factor(df_HCT116_crispr_kras$GeneType))[1],
                                         which_tail='t', 
                                         compare_FC.ranks_mode=T,
                                         df_HCT116_crispr_kras=df_HCT116_crispr){
  if(compare_FC.ranks_mode==F){
    FC_inMutant=df_HCT116_crispr_kras$WT.log2.fold.change[df_HCT116_crispr_kras$GeneType==cde_type]
    FC_inWT=df_HCT116_crispr_kras$KRas.mutant.log2.fold.change[df_HCT116_crispr_kras$GeneType==cde_type]
  }else {
    FC_inMutant=df_HCT116_crispr_kras$Rank_Mut[df_HCT116_crispr_kras$GeneType==cde_type]
    FC_inWT=df_HCT116_crispr_kras$Rank_WT[df_HCT116_crispr_kras$GeneType==cde_type]
  }
  wilcox.test(FC_inMutant, FC_inWT, paired = T, alternative = which_tail)$p.value
}

# *Results* : Testing CDE+/- essentiality direction in CRISPR DLD1 screening
# Using direct FC score
test_direction_krascrispr_DLD1(cde_type="KRAS CDE Neg",
                               which_tail='l',
                               compare_FC.ranks_mode=F,
                               df_HCT116_crispr_kras=df_HCT116_crispr)
test_direction_krascrispr_DLD1(cde_type="KRAS CDE Pos",
                               which_tail='g',
                               compare_FC.ranks_mode=F,
                               df_HCT116_crispr_kras=df_HCT116_crispr)
# Using direct FC rank score
test_direction_krascrispr_DLD1(cde_type="KRAS CDE Neg",
                               which_tail='l',
                               compare_FC.ranks_mode=T,
                               df_HCT116_crispr_kras=df_HCT116_crispr)
test_direction_krascrispr_DLD1(cde_type="KRAS CDE Pos",
                               which_tail='g',
                               compare_FC.ranks_mode=T,
                               df_HCT116_crispr_kras=df_HCT116_crispr)

