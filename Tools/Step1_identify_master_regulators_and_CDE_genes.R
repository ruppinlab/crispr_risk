# Version for Github
source('Step0_Globally_used_Functions_and_Datasets.R')

####################################################################################
##Identifying the DE, CDE genes for TP53 gene
####################################################################################
## This function return four a list of DE+/- for a given master regulator in shrna 
## and crispr with effect size.
p53_DE=Testing_CRISPR_damage_bias(GeneName='TP53')
##Identifying the size of DE genes
#Size of CRISPR DE+ genes
print(paste('Size of CRISPR DE+ genes', nrow(p53_DE$CRISPR_DE_pos)))
#Size of CRISPR DE- genes
print(paste('Size of CRISPR DE- genes', nrow(p53_DE$CRISPR_DE_neg)))
#Size of shrna DE+ genes
print(paste('Size of shrna DE+ genes', nrow(p53_DE$shrna_DE_pos)))
#Size of shrna DE- genes
print(paste('Size of shrna DE- genes', nrow(p53_DE$shrna_DE_neg)))
###Following are the CDE+ and CDE- genes of p53 genes 
p53_CDE_pos= p53_DE$CRISPR_DE_pos$GeneName[p53_DE$CRISPR_DE_pos$GeneName %!in% p53_DE$shrna_DE_pos$GeneName]
p53_CDE_neg= p53_DE$CRISPR_DE_neg$GeneName[p53_DE$CRISPR_DE_neg$GeneName %!in% p53_DE$shrna_DE_neg$GeneName]
length(p53_CDE_pos)
length(p53_CDE_neg)


####################################################################################
##Calculating the probability of imbalance of DE+/- in crispr  in comparison to control (shrna)
####################################################################################
p53_Contigency=t(sapply(p53_DE, nrow))
p53_DE_imbalance_prob=significance_test(p53_Contigency)


####################################################################################
## Using the above framework, we'd identify DE+/- genes for cancer driver genes defined
## in Vogestein et al.
####################################################################################
Vog_Gene_df=get(load('../Data/vogelstein.RData'))
Vog_Genes=Vog_Gene_df$symbol
# the code below is for computing the DE genes for all cancer genes from Vogelstein et al., the result is save in Vog_Genes_DE
# this is saved in ../Data/DE_posNneg_genes_for_VolgGenes.RDS and loaded below to save running time
#Vog_Genes_DE= lapply(Vog_Genes, function(x) Testing_CRISPR_damage_bias(GeneName=x))
#names(Vog_Genes_DE)=Vog_Genes
#saveRDS(Vog_Genes_DE, '../Data/DE_posNneg_genes_for_VolgGenes.RDS')
#To quickly repliacte, instead, one can load Volg_CDE genes for each vog Genes
Vog_Genes_DE=readRDS('../Data/DE_posNneg_genes_for_VolgGenes.RDS')
Contigency=t(sapply(Vog_Genes_DE, function(x) sapply(x, nrow)))
rownames(Contigency)=Vog_Genes
# Significance of imblance of DE+/- in CRISPR and shRNA
Sig=p.adjust(unlist(apply(Contigency, 1, function(x)  err_handle(significance_test(x)) )), method='fdr')
# create a df of log(Sig) and potential Master regulator tested
Log_Sig=-log(Sig, 10)
Log_Sig_df=data.frame(Vog_Genes, Log_Sig)
Log_Sig_df=Log_Sig_df[order(Log_Sig_df$Log_Sig),]
Log_Sig_df$Log_Sig[Log_Sig_df$Log_Sig>300]=300
# Log_Sig_df is a dataframe comprising all the tested potential master regulators names
# with their significnace of imbalance of DE+/- in crispr and shRNA screens


####################################################################################
# write CDEs_of_ThreeMaster_Regulators which are significant in above test:
# TP53, KRAS and VHL
#Creating: CDEs_of_ThreeMaster_Regulators.csv
####################################################################################
CDEs_of_ThreeMaster_Regulators=list(P53.CDE.Pos= as.character(Vog_Genes_DE$'TP53'$CRISPR_DE_pos$GeneName[Vog_Genes_DE$'TP53'$CRISPR_DE_pos$GeneName
                                                                                   %!in% Vog_Genes_DE$'TP53'$shrna_DE_pos$GeneName]),
                                    P53.CDE.Neg= as.character(Vog_Genes_DE$'TP53'$CRISPR_DE_neg$GeneName[Vog_Genes_DE$'TP53'$CRISPR_DE_neg$GeneName
                                                                                   %!in% Vog_Genes_DE$'TP53'$shrna_DE_neg$GeneName]),
                                    VHL.CDE.Pos= as.character(Vog_Genes_DE$'VHL'$CRISPR_DE_pos$GeneName[Vog_Genes_DE$'VHL'$CRISPR_DE_pos$GeneName
                                                                                    %!in% Vog_Genes_DE$'VHL'$shrna_DE_pos$GeneName]),
                                    VHL.CDE.Neg= as.character(Vog_Genes_DE$'VHL'$CRISPR_DE_neg$GeneName[Vog_Genes_DE$'VHL'$CRISPR_DE_neg$GeneName
                                                                                    %!in% Vog_Genes_DE$'VHL'$shrna_DE_neg$GeneName]),
                                    KRAS.CDE.Pos= as.character(Vog_Genes_DE$'KRAS'$CRISPR_DE_pos$GeneName[Vog_Genes_DE$'KRAS'$CRISPR_DE_pos$GeneName
                                                                                    %!in% Vog_Genes_DE$'KRAS'$shrna_DE_pos$GeneName]),
                                    KRAS.CDE.Neg= as.character(Vog_Genes_DE$'KRAS'$CRISPR_DE_neg$GeneName[Vog_Genes_DE$'KRAS'$CRISPR_DE_neg$GeneName
                                                                                    %!in% Vog_Genes_DE$'KRAS'$shrna_DE_neg$GeneName]))
CDEs_of_ThreeMaster_Regulators_df=cbind.fill(CDEs_of_ThreeMaster_Regulators[[1]], CDEs_of_ThreeMaster_Regulators[[2]],
                                          CDEs_of_ThreeMaster_Regulators[[3]], CDEs_of_ThreeMaster_Regulators[[4]],
                                          CDEs_of_ThreeMaster_Regulators[[5]], CDEs_of_ThreeMaster_Regulators[[6]],
                                          fill = '' )
colnames(CDEs_of_ThreeMaster_Regulators_df)=names(CDEs_of_ThreeMaster_Regulators)
write.csv(CDEs_of_ThreeMaster_Regulators_df, '../Data/CDEs_of_ThreeMaster_Regulators.csv', quote=F, row.names = F)


####################################################################################
# get the CDE+/- genes for all the cancer genes: saved in ../Data/CDE_for_volg.RDS
####################################################################################
Vog_Genes_DE=readRDS('../Data/DE_posNneg_genes_for_VolgGenes.RDS')
xs <- names(Vog_Genes_DE)
names(xs) <- xs
CDE_for_volg=lapply(xs, function(K)
  list(CDE.pos= as.character(err_handle(Vog_Genes_DE[[K]]$CRISPR_DE_pos$GeneName[Vog_Genes_DE[[K]]$CRISPR_DE_pos$GeneName
                                                    %!in% Vog_Genes_DE[[K]]$shrna_DE_pos$GeneName])),
       CDE.Neg= as.character(err_handle(Vog_Genes_DE[[K]]$CRISPR_DE_neg$GeneName[Vog_Genes_DE[[K]]$CRISPR_DE_neg$GeneName
                                                    %!in% Vog_Genes_DE[[K]]$shrna_DE_neg$GeneName]))))

names(CDE_for_volg)=Vog_Genes
saveRDS(CDE_for_volg, '../Data/CDE_for_volg.RDS')

