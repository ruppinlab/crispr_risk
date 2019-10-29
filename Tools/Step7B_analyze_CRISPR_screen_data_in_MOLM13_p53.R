# Version for Github
source('Step0_Globally_used_Functions_and_Datasets.R')
source('Step6_CDE_genes_for_MOLM13_cell_line.R')
# Processing pooled CRISPR-KO screening using a CRISPR
# library targeting the top p53 CDE+ and CDE- genes
# in isogenic p53 wildtype and p53 mutant pair of the 
# MOLM13 leukemia cell line.

# Here, we tested whether CDE-Positive genes are differentially more ess in P53-WT
# vs Mut in CRISPR-Cas9 based screenings and not in CRISPRi.
# Tested the above test in the opposite direction for CDE-Negative.
############################################################################
#Load and Pre-Processing
############################################################################
# Load Off-Target score for each sgRNAs
offTarget_Scores=readRDS('../Data/sgRNA_Off_Target_Scores.RDS')
offTarget_Scores$gRNA=substring(offTarget_Scores$gRNAsPlusPAM, first=0,last=20)

# Loading Previously identified CDE genes
CDE=read.csv('../Data/CDEs_of_ThreeMaster_Regulators.csv')
############################################################################
# Functions neeeded during the Pre-processing
############################################################################
# The following function takes raw CRISPR-Cas9 screening counts and provides top 10% sgRNAs
# targeting CDE+ & CDE- after removing reads with 1. low Day0 counts,
# 2. present in avana/DepMap dataset 3. Induces more DNA damage.

# P53 - CRSIPR-Cas9
P53_ProcessingD0<-function(removesgRNA_Threshold,    # Remove the sgRNAs with Day0 counts threshold
                           TopXPercent,              # Take top X percent CDE genes
                           whether_paired,           # Take wilcox test is paired or not
                           CDE=MOLM,                 # Consider the MOLM specific ranking
                           only_DepMap=T,            # Only consider sgRNAs in previous Screening
                           include_offTarget=T       # Include Off-Target Score to improve Signal
                           ){
  p53WT=read.csv('../Data/Isogenic_Screenings/p53_CDEScreening_WildType.csv')
  p53Mutant=read.csv('../Data/Isogenic_Screenings/p53_CDEScreening_Mutant.csv')
  colnames(p53Mutant)=colnames(p53WT)
  
  # Step 0: Pre-Processing
  p53WT$avg_D0=rowMeans(p53WT[,c("R1_DO_cpm","R2_DO_cpm")])
  p53WT$avg_D30=rowMeans(p53WT[,c("R1_D3O_cpm","R2_D3O_cpm")])
  p53WT$FC= p53WT$avg_D30/p53WT$avg_D0
  p53Mutant$avg_D0=rowMeans(p53Mutant[,c("R1_DO_cpm","R2_DO_cpm")])
  p53Mutant$avg_D30=rowMeans(p53Mutant[,c("R1_D3O_cpm","R2_D3O_cpm")])
  p53Mutant$FC= p53Mutant$avg_D30/p53Mutant$avg_D0
  
  # Define rank
  p53Mutant$FC_rank=rank(p53Mutant$FC)
  p53WT$FC_rank=rank(p53WT$FC)
  
  p53Mutant$top100OfftargetTotalScore=offTarget_Scores$top100OfftargetTotalScore[
    match(p53Mutant$sgRNA, offTarget_Scores$gRNA)]
  p53Mutant$top50OfftargetTotalScore=offTarget_Scores$top5OfftargetTotalScore[
    match(p53Mutant$sgRNA, offTarget_Scores$gRNA)]
  
  if(only_DepMap){
    p53WT=p53WT[!is.na(p53Mutant$top100OfftargetTotalScore),]
    p53Mutant=p53Mutant[!is.na(p53Mutant$top100OfftargetTotalScore),]
  }
  
  # CDE specific rank
  CDE=lapply(CDE, function(x) as.character(x[x!='']) )
  CDE$P53.CDE.Pos=CDE$P53.CDE.Pos[!is.na(match(CDE$P53.CDE.Pos,
                                                   p53Mutant$Gene[p53Mutant$Gene.group=='p53 CDE Pos']))]
  CDE$P53.CDE.Neg=CDE$`P53.CDE.Neg`[!is.na(match(CDE$P53.CDE.Neg,
                                                   p53Mutant$Gene[p53Mutant$Gene.group=='p53 CDE Neg']))]
  
  # remove sgRNAs with Day-0 READS less than the threshold
  from_mut=p53Mutant$R1_DO_Count<removesgRNA_Threshold | p53Mutant$R2_DO_Count<removesgRNA_Threshold
  from_WT=p53WT$R1_DO_Count<removesgRNA_Threshold | p53WT$R2_DO_Count<removesgRNA_Threshold
  p53WT=p53WT[!(from_WT | from_mut),]
  p53Mutant=p53Mutant[!(from_WT | from_mut),]
  
  # Using MOLM-Specific CDE rank
  p53Mutant$CDE_Rank=NA
  p53Mutant$CDE_Rank[p53Mutant$Gene.group=='p53 CDE Pos']= match(p53Mutant$Gene[p53Mutant$Gene.group=='p53 CDE Pos'],
                                                                 CDE$`P53.CDE.Pos`)
  p53Mutant$CDE_Rank[p53Mutant$Gene.group=='p53 CDE Neg']= match(p53Mutant$Gene[p53Mutant$Gene.group=='p53 CDE Neg'],
                                                                 CDE$`P53.CDE.Neg`)
  # Calculate the topX hits count
  topX_neg=ceiling(length(unique(p53Mutant$Gene[p53Mutant$Gene.group=="p53 CDE Neg"]))*TopXPercent)
  topX_pos=ceiling(length(unique(p53Mutant$Gene[p53Mutant$Gene.group=="p53 CDE Pos"]))*TopXPercent)
  # Considering Off-Target scores of sgRNA to choose sgRNAs inducing more DNA damage
  if (include_offTarget) {
    one_third_point=summary(p53Mutant$top100OfftargetTotalScore)[2]
    two_third_point=summary(p53Mutant$top100OfftargetTotalScore)[3]
    GOI = (p53Mutant$CDE_Rank < topX_neg &  p53Mutant$Gene.group == 'p53 CDE Neg' &
             (is.na(p53Mutant$top100OfftargetTotalScore) |
                p53Mutant$top100OfftargetTotalScore<one_third_point ) ) |
      (p53Mutant$CDE_Rank < topX_pos &  p53Mutant$Gene.group == 'p53 CDE Pos' &
         (is.na(p53Mutant$top100OfftargetTotalScore) |
            p53Mutant$top100OfftargetTotalScore>two_third_point ))
  } else {
    GOI = (p53Mutant$CDE_Rank < topX_neg &  p53Mutant$Gene.group == 'p53 CDE Neg') | (p53Mutant$CDE_Rank < topX_pos &  p53Mutant$Gene.group == 'p53 CDE Pos')
  }
  
  df2plot=rbind(data.frame(geneName=p53WT$Gene[GOI],
                           RANK=p53WT$FC_rank[GOI],
                           type='WT',
                           group=p53Mutant$Gene.group[GOI]),
                data.frame(geneName=p53WT$Gene[GOI],
                           RANK=p53Mutant$FC_rank[GOI],
                           type='Mut',
                           group=p53Mutant$Gene.group[GOI]))
  # This is the one-sided Significances for diff in ess in WT vs Mut
  # in the right directiomn for CDE+ & CDE- in order
  p_value=c(wilcox.test(df2plot$RANK[df2plot$type=='WT' & df2plot$group=='p53 CDE Pos'],
                        df2plot$RANK[df2plot$type=='Mut' & df2plot$group=='p53 CDE Pos'],
                        alternative = 'l', paired=whether_paired)$p.value,
            wilcox.test(df2plot$RANK[df2plot$type=='WT' & df2plot$group=='p53 CDE Neg'],
                        df2plot$RANK[df2plot$type=='Mut' & df2plot$group=='p53 CDE Neg'],
                        alternative = 'g', paired=whether_paired)$p.value)
  
  
  # A=ggplot(df2plot, aes(x=group, fill=type, y=RANK))+geom_boxplot()+
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  #   annotate(geom="text", y=0, x=0.5, label=round(p_value[2], 2), color="blue")+
  #   annotate(geom="text", y=0, x=1.5, label=round(p_value[1], 2), color="blue")
  df2plot_genelevel=do.call(rbind,
                            lapply(split(df2plot, 
                                         list(as.character(df2plot$geneName),
                                              as.character(df2plot$type))),
                                   function(x) data.frame(geneName=unique(x$geneName),
                                                          RANK=mean(x$RANK),
                                                          group=unique(x$group),
                                                          type=unique(x$type)
                                   )))
  p_value_genelevel=c(wilcox.test(df2plot_genelevel$RANK[df2plot_genelevel$type=='WT'
                                                         & df2plot_genelevel$group=='p53 CDE Pos'],
                                  df2plot_genelevel$RANK[df2plot_genelevel$type=='Mut'
                                                         & df2plot_genelevel$group=='p53 CDE Pos'],
                                  alternative = 'l', paired=whether_paired)$p.value,
                      wilcox.test(df2plot_genelevel$RANK[df2plot_genelevel$type=='WT'
                                                         & df2plot_genelevel$group=='p53 CDE Neg'],
                                  df2plot_genelevel$RANK[df2plot_genelevel$type=='Mut'
                                                         & df2plot_genelevel$group=='p53 CDE Neg'],
                                  alternative = 'g', paired=whether_paired)$p.value)
  
  # B=ggplot(df2plot_genelevel, aes(x=group, fill=type, y=RANK))+geom_boxplot()+
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  #   annotate(geom="text", y=0, x=0.5, label=round(p_value_genelevel[2], 2), color="blue")+
  #   annotate(geom="text", y=0, x=1.5, label=round(p_value_genelevel[1], 2), color="blue")
  # grid.arrange(A, B, nrow = 2)
  print(c(sgRNA_level_Significance=p_value, gene_level_Significance=p_value_genelevel))
  df2plot
}

# The following function takes raw CRISPRi screening counts and provides top 10% sgRNAs
# targeting CDE+ & CDE- after removing reads with 1. low Day0 counts.
# iP53 - CRISPRi
iP53_ProcessingD0<-function(removesgRNA_Threshold,
                            TopXPercent,              # Take top X percent CDE genes                            
                            whether_paired=T,
                            CDE_bckp=MOLM){
  ip53WT=read.csv('../Data/Isogenic_Screenings/ip53_WT.csv')
  ip53Mutant=read.csv('../Data/Isogenic_Screenings/ip53_Mut.csv')
  
  # Step 0: Pre-Processing
  ip53WT$avg_D0=rowMeans(ip53WT[,c("cpm_D0_R1","cpm_D0_R2")])
  ip53WT$avg_D30=rowMeans(ip53WT[,c("cpm_D30_R1","cpm_D30_R2")])
  ip53WT$FC= ip53WT$avg_D30/ip53WT$avg_D0
  ip53Mutant$avg_D0=rowMeans(ip53Mutant[,c("cpm_D0_R1","cpm_D0_R2")])
  ip53Mutant$avg_D30=rowMeans(ip53Mutant[,c("cpm_D30_R1","cpm_D30_R2")])
  ip53Mutant$FC= ip53Mutant$avg_D30/ip53Mutant$avg_D0

  # Define rank
  ip53Mutant$FC_rank=rank(ip53Mutant$FC)
  ip53WT$FC_rank=rank(ip53WT$FC)
  
  # MOLM specific rank
  CDE=lapply(CDE, function(x) as.character(x[x!='']) )
  CDE$P53.CDE.Pos=CDE$P53.CDE.Pos[!is.na(match(CDE$P53.CDE.Pos,
                                                   ip53Mutant$Gene[ip53Mutant$CDE.Group=='p53 CDE Pos']))]
  CDE$P53.CDE.Neg=CDE$P53.CDE.Neg[!is.na(match(CDE$P53.CDE.Neg,
                                                   ip53Mutant$Gene[ip53Mutant$CDE.Group=='p53 CDE Neg']))]
  
  # Remove reads less than a threshold at Day0
  from_mut=ip53Mutant$Count_D0_R1<removesgRNA_Threshold | 
    ip53Mutant$Count_D0_R2<removesgRNA_Threshold
  from_WT=ip53WT$Count_D0_R1<removesgRNA_Threshold |
    ip53WT$Count_D0_R2<removesgRNA_Threshold
  ip53WT=ip53WT[!(from_WT | from_mut),]
  ip53Mutant=ip53Mutant[!(from_WT | from_mut),]
  
  # Define Ranks in screening using MOLM specific rank
  ip53Mutant$CDE_Rank=NA
  ip53Mutant$CDE_Rank[ip53Mutant$CDE.Group=='p53 CDE Pos']= 
    match(ip53Mutant$Gene[ip53Mutant$CDE.Group=='p53 CDE Pos'], CDE$P53.CDE.Pos)
  ip53Mutant$CDE_Rank[ip53Mutant$CDE.Group=='p53 CDE Neg']= 
    match(ip53Mutant$Gene[ip53Mutant$CDE.Group=='p53 CDE Neg'], CDE$P53.CDE.Neg)
  
  ###Calculate the topX hits count
  topX_neg=ceiling(length(unique(ip53Mutant$Gene[ip53Mutant$CDE.Group=="p53 CDE Neg"]))*TopXPercent) #!!! I think it's here that p53Mutant is not found, do you mean ip53Mutant?
  topX_pos=ceiling(length(unique(ip53Mutant$Gene[ip53Mutant$CDE.Group=="p53 CDE Pos"]))*TopXPercent) #!!! ditto
  
  # Top X Genes based on the above 
  GOI = (ip53Mutant$CDE_Rank < topX_neg &  ip53Mutant$CDE.Group == 'p53 CDE Neg') |
    (ip53Mutant$CDE_Rank < topX_pos &  ip53Mutant$CDE.Group == 'p53 CDE Pos')
  df2plot=rbind(data.frame(geneName=ip53WT$Gene[GOI],
                           RANK=ip53WT$FC_rank[GOI],
                           type='WT',
                           group=ip53Mutant$CDE.Group[GOI]),
                data.frame(geneName=ip53WT$Gene[GOI],
                           RANK=ip53Mutant$FC_rank[GOI],
                           type='Mut',
                           group=ip53Mutant$CDE.Group[GOI]))
  # This is the one-sided Significances for diff in ess in WT vs Mut
  # in the right directiomn for CDE+ & CDE- in order
  p_value=c(wilcox.test(df2plot$RANK[df2plot$type=='WT' & df2plot$group=='p53 CDE Pos'],
                        df2plot$RANK[df2plot$type=='Mut' & df2plot$group=='p53 CDE Pos'],
                        alternative = 'l', paired=whether_paired)$p.value,
            wilcox.test(df2plot$RANK[df2plot$type=='WT' & df2plot$group=='p53 CDE Neg'],
                        df2plot$RANK[df2plot$type=='Mut' & df2plot$group=='p53 CDE Neg'],
                        alternative = 'g', paired=whether_paired)$p.value)
  
  # A=ggplot(df2plot, aes(x=group, fill=type, y=RANK))+geom_boxplot()+
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  #   annotate(geom="text", y=0, x=0.5, label=round(p_value[2], 2), color="blue")+
  #   annotate(geom="text", y=0, x=1.5, label=round(p_value[1], 2), color="blue")
  df2plot_genelevel=do.call(rbind,
                            lapply(split(df2plot,
                                         list(as.character(df2plot$geneName),
                                              as.character(df2plot$type))),
                                   function(x)
                                     data.frame(geneName=unique(x$geneName),
                                                RANK=mean(x$RANK), 
                                                group=unique(x$group),
                                                type=unique(x$type)
                                                )))
  p_value_genelevel=c(wilcox.test(df2plot_genelevel$RANK[df2plot_genelevel$type=='WT' &
                                                           df2plot_genelevel$group=='p53 CDE Pos'],
                                  df2plot_genelevel$RANK[df2plot_genelevel$type=='Mut' &
                                                           df2plot_genelevel$group=='p53 CDE Pos'],
                                  alternative = 'l', paired=whether_paired)$p.value,
                      wilcox.test(df2plot_genelevel$RANK[df2plot_genelevel$type=='WT' &
                                                           df2plot_genelevel$group=='p53 CDE Neg'],
                                  df2plot_genelevel$RANK[df2plot_genelevel$type=='Mut' &
                                                           df2plot_genelevel$group=='p53 CDE Neg'],
                                  alternative = 'g', paired=whether_paired)$p.value)
  
  # B=ggplot(df2plot_genelevel, aes(x=group, fill=type, y=RANK))+geom_boxplot()+
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  #   annotate(geom="text", y=0, x=0.5, label=round(p_value_genelevel[2], 2), color="blue")+
  #   annotate(geom="text", y=0, x=1.5, label=round(p_value_genelevel[1], 2), color="blue")
  # grid.arrange(A, B, nrow = 2)
  print(c(sgRNA_level_Significance=p_value, 
          gene_level_Significance=p_value_genelevel))
  df2plot
}
############################################################################
#Calling functions
############################################################################
# Following will return a data.frame of top 10% (TopXPercent) processed sgRNA hits
# from CRISPR-Cas9 screenings. In additions, it prints significance of whether the
# difference of 1. CDE+ genes KO viability is greater in Mutant vs WT. 2. CDE- genes
# KO viability is greater in WT vs Mutant. This significance is repeated at sgRNA and
# gene-level.
pooled_crisprCas9_p53=P53_ProcessingD0(removesgRNA_Threshold=20,
                                   TopXPercent=0.1,
                                   whether_paired = T, 
                                   CDE = MOLM,
                                   only_DepMap = T,
                                   include_offTarget = T)
# Following will return top 10% (TopXPercent) processed sgRNA hits from CRISPRi screenings
pooled_crisprI_p53=iP53_ProcessingD0(removesgRNA_Threshold=20,
                                     TopXPercent=0.1,
                                     whether_paired = T,
                                     CDE = MOLM)


