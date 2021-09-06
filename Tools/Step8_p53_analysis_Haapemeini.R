# Version for Github
source('Step0_Globally_used_Functions_and_Datasets.R')

# Here, we are testing the robust enrichment of our CDE+/- in the DE+/- identified
# in Primary cells from Haapemeini et al.
################################################################
##Load Files
################################################################
##Loading CRISPR Screenings from Haapemeini et al
Mut_R1=read.csv('../Data/Isogenic_Screenings/Genome wide screening data RPE/R1-null.gene_summary.txt', sep='\t')
Mut_R2=read.csv('../Data/Isogenic_Screenings/Genome wide screening data RPE/R2-null.gene_summary.txt', sep='\t')
WT_R1=read.csv('../Data/Isogenic_Screenings/Genome wide screening data RPE/R1-WT.gene_summary.txt', sep='\t')
WT_R2=read.csv('../Data/Isogenic_Screenings/Genome wide screening data RPE/R2-WT.gene_summary.txt', sep='\t')
# Our CDE+/- genes of identified from avana
CDE=read.csv('../Data/CDEs_of_ThreeMaster_Regulators.csv')
CDE=lapply(CDE, function(x) as.character(x[x!='']) )
DE_genes=readRDS('../Data/DE_genes.RDS')
################################################################
##Prep: Step 0
################################################################
WT_R1=WT_R1[!is.na(match(as.character(unlist(WT_R1$id)), rownames(avana))), ]
WT_R2=WT_R2[!is.na(match(as.character(unlist(WT_R2$id)), rownames(avana))), ]
Mut_R1=Mut_R1[!is.na(match(as.character(unlist(Mut_R1$id)), rownames(avana))), ]
Mut_R2=Mut_R2[!is.na(match(as.character(unlist(Mut_R2$id)), rownames(avana))), ]
Mut_R2=Mut_R2[match(Mut_R1$id, Mut_R2$id),]
Mut_avgLFC=data.frame(geneName=Mut_R2$id, avgLFC=rowMeans(data.frame(Mut_R2$neg.lfc, Mut_R1$neg.lfc )) )
WT_R2=WT_R2[match(WT_R1$id, WT_R2$id),]
WT_avgLFC=data.frame(geneName=WT_R2$id, avgLFC=rowMeans(data.frame(WT_R2$neg.lfc, WT_R1$neg.lfc)) )
WT_avgLFC=WT_avgLFC[match(Mut_avgLFC$geneName, WT_avgLFC$geneName),]
df_p53Haap=data.frame(WT_avgLFC, Mut=Mut_avgLFC)
# Create FC Rank Rank
df_p53Haap$Rank_WT=rank(df_p53Haap$avgLFC)  
df_p53Haap$Rank_Mut=rank(df_p53Haap$Mut.avgLFC)
#Ranking DE genes
DE_genes$DE_pos_KRAS=DE_genes$DE_pos_KRAS[order(DE_genes$DE_pos_KRAS$Sig..1.....2.),]
DE_genes$DE_neg_KRAS=DE_genes$DE_neg_KRAS[order(DE_genes$DE_neg_KRAS$Sig..2.....2.),]
DE_genes$DE_pos_TP53=DE_genes$DE_pos_TP53[order(DE_genes$DE_pos_TP53$Sig_TP53..1.....2.),]
DE_genes$DE_neg_TP53=DE_genes$DE_neg_TP53[order(DE_genes$DE_neg_TP53$Sig_TP53..2.....2.),]
################################################################
# Comparison of CDE genes at the top/bottom
################################################################
df_p53Haap=df_p53Haap[order(df_p53Haap$avgLFC-df_p53Haap$Mut.avgLFC),]
test_enrich_givenK<-function(flag1=100, flag2=100){
  c(hypergeometric_test_for_twolists(test_list = as.character(unlist(head(df_p53Haap$geneName, flag1))),
                                     base_list = head(as.character(unlist(DE_genes$DE_pos_TP53[,1])), flag2),
                                     global = as.character(unlist(df_p53Haap$geneName)) ),
    hypergeometric_test_for_twolists(test_list = as.character(unlist(tail(df_p53Haap$geneName, flag1))),
                                     base_list = head(as.character(unlist(DE_genes$DE_neg_TP53[,1])), flag2),
                                     global = as.character(unlist(df_p53Haap$geneName)) )
  )
}

# Testing enrichment robustly:
# A) We took top X DE+ genes from Haap screens and top Y DE+ genes from 
# identified previously from avana and calculated an enrichment of their overlap.
# B) We repeated this for X ranging from 50 to size of Haap DE+ genes and Y from
# 50 to size of avana DE+ with an interval of 50.
# C) We finally calculated the proportion of times this enrichemnt is significant
# with FDR<0.2.
# *All the above steps are repeated for DE-/DE- as well
Pos=data.frame(sapply(seq(100, nrow(DE_genes$DE_pos_TP53), 50), function(y)
  sapply(seq(100, nrow(DE_genes$DE_pos_TP53), 50), 
         function(x) test_enrich_givenK(x, y ))[1,]))
#Final robust enrichment score of CDE+
sum(Pos<0.2)/length(Pos<0.2)

Neg=data.frame(sapply(seq(50, nrow(DE_genes$DE_neg_TP53), 50), function(y)
  sapply(seq(10, nrow(DE_genes$DE_neg_TP53), 10), 
             function(x) test_enrich_givenK(x, y ))[2,]))
#Final robust enrichment score of CDE-
sum(Neg<0.2)/length(Neg<0.2)
