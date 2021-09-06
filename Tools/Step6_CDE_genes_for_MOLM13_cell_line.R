# Version for Github
source('Step0_Globally_used_Functions_and_Datasets.R')

# p53/KRAS CDE+/- ranked for the MOLM13 cell line

MOLM=list()
MOLM$p53_CDE_Neg=c(as.character(read.csv('../Data/MOLM_CDE_Genes/P53/CDE_negative/S2A_CDE_neg_MOLM13_specific_overlap.csv')[,1]),
                   as.character(read.csv('../Data/MOLM_CDE_Genes/P53/CDE_negative/S2B_CDE_restof_neg_MOLM13_specific.csv')[,1]))
MOLM$p53_CDE_Pos=as.character(unlist(read.csv('../Data/MOLM_CDE_Genes/P53/CDE_Positive/S1_CDE_Pos.csv')[,1]))


MOLM$KRAS_CDE_Neg=c(as.character(read.csv('../Data/MOLM_CDE_Genes/KRAS/CDE_neg/CDE_Neg_KRAS.csv')[,1]),
                   as.character(read.csv('../Data/MOLM_CDE_Genes/KRAS/CDE_neg/CDE_Neg_rest_KRAS.csv')[,1]))
MOLM$KRAS_CDE_Pos=c(as.character(unlist(read.csv('../Data/MOLM_CDE_Genes/KRAS/CDE_pos/CDE_Pos_KRAS.csv')[,1])),
                   as.character(unlist(read.csv('../Data/MOLM_CDE_Genes/KRAS/CDE_pos/CDE_Pos_rest_KRAS.csv')[,1])))

names(MOLM)= c('P53.CDE.Neg', 'P53.CDE.Pos', 'KRAS.CDE.Neg', 'KRAS.CDE.Pos')

