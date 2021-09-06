# Version for Github
source('Step0_Globally_used_Functions_and_Datasets.R')

# This script is for analyzing the CDE effect (i.e. the imbanlance of the number of DE+/- genes only in CRISPR but not in shRNA screens) of p53
# based on only TP53-LOF mutations

#Mutation matirx::CCLE
sample_info=read.csv('../Data/sample_info.csv')

# Rather than loading all all the variants here, we only load TP53 mutation profile
# To get CCLE_mutations.csv. Please refer to https://depmap.org/portal/
# Mut_raw=read.csv('../Data/CCLE_mutations.csv', sep=',')
# ##Filters to remove the type of mutation one wants:: Currently Silent Mutations are removed.
# Mut_Filtered = Mut_raw[which(!(Mut_raw$Variant_Classification=='Silent')),]
# GeneName='TP53'
# TP53_Mut = Mut_Filtered[which(Mut_Filtered$Hugo_Symbol==GeneName),]
# write.csv(TP53_Mut,'../Data/CCLE_mutations_TP53.csv')
GeneName='TP53'
TP53_Mut=read.csv('../Data/CCLE_mutations_TP53.csv')
avana_colnames_depmapID=sample_info$DepMap_ID[match(gsub('X','',colnames(avana)),
                                                    gsub('X','',sample_info$CCLE_name))]
TP53_Mut_avana=TP53_Mut[na.omit(match(avana_colnames_depmapID,
                                      TP53_Mut$Tumor_Sample_Barcode)),]
TP53_LOF_cellLines= TP53_Mut_avana$Tumor_Sample_Barcode[which(TP53_Mut_avana$isDeleterious)]

# Pre-computed results are saved in the Data folder
# saveRDS(TP53_LOF_cellLines,'../Data/TP53_LOFs.RDS')
# saveRDS(TP53_Mut_avana,'../Data/TP53_Mut_avana.RDS')
TP53_LOFs=readRDS('../Data/TP53_LOFs.RDS')
TP53_Mut_avana=readRDS('../Data/TP53_Mut_avana.RDS')

#Creating P53-LOF status
#initialize as all WT
LOF_Info = rep(0, ncol(avana))
##Make all the cell lines with any mutation NA
Cell_Lines_wd_Mut = TP53_Mut_avana$Tumor_Sample_Barcode
Cell_Lines_wd_Mut_id=match(gsub('X','',sample_info$CCLE_name[
  match(Cell_Lines_wd_Mut, sample_info$DepMap_ID)]), 
  gsub('X','',colnames(avana)))
LOF_Info[Cell_Lines_wd_Mut_id]=NA
Cell_Lines_wd_LOF=match(sample_info$CCLE_name[
  match(TP53_LOF_cellLines, sample_info$DepMap_ID)], colnames(avana))
LOF_Info[Cell_Lines_wd_LOF] = 1

Mut_CCLE_LOF = t(data.frame(TP53=as.numeric(LOF_Info)))

#Comparing WT vs p53-Mut-LOF
p53_DE_usingLOFs = Testing_CRISPR_damage_bias(GeneName, Feature_to_test=Mut_CCLE_LOF) #!!! error, see below
#The above variable comprises DE+/- genes for p53 in crispr and shRNA screens.

# testing for the CDE effect
p53_Contigency=t(sapply(p53_DE_usingLOFs, nrow))
p53_DE_imbalance_prob=significance_test(p53_Contigency)

