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

# KRAS - CRSIRP-Cas9
KRAS_ProcessingD0<-function(removesgRNA_Threshold,    # Remove the sgRNAs with Day0 counts threshold
                            topX_neg,                 # Take top X percent CDE genes
                            topX_pos,
                            whether_matched=T,
                            CDE,
                            only_DepMap=T,
                            include_offTarget=F){
  KRAS_Mut=read.csv('../Data/Isogenic_Screenings/KRAS_Mut.csv')
  KRAS_WT2=read.csv('../Data/Isogenic_Screenings/KRAS_WT2.csv')
  
  KRAS_Mut$FC = KRAS_Mut$Avg.D30/ KRAS_Mut$Avg.D0
  KRAS_WT2$FC = KRAS_WT2$Avg.D30/KRAS_WT2$Avg.D0
  colnames(KRAS_WT2)=colnames(KRAS_Mut)
  CDE=lapply(CDE, function(x) as.character(x[x!='']) )
  CDE$KRAS.CDE.Pos=CDE$KRAS.CDE.Pos[!is.na(match(CDE$KRAS.CDE.Pos,
                                                     KRAS_Mut$Gene[KRAS_Mut$Gene.group=='KRAS CDE Pos']))]
  CDE$KRAS.CDE.Neg=CDE$KRAS.CDE.Neg[!is.na(match(CDE$KRAS.CDE.Neg,
                                                     KRAS_Mut$Gene[KRAS_Mut$Gene.group=='KRAS CDE Neg']))]
  KRAS_WT2$FC_rank=rank(KRAS_WT2$FC)
  KRAS_Mut$FC_rank=rank(KRAS_Mut$FC)
  
  KRAS_Mut$top100OfftargetTotalScore=offTarget_Scores$top100OfftargetTotalScore[
    match(KRAS_Mut$sgRNA, offTarget_Scores$gRNA)]
  KRAS_Mut$top50OfftargetTotalScore=offTarget_Scores$top5OfftargetTotalScore[
    match(KRAS_Mut$sgRNA, offTarget_Scores$gRNA)]

  if(only_DepMap){
    KRAS_WT2=KRAS_WT2[!is.na(KRAS_Mut$top100OfftargetTotalScore),]
    KRAS_Mut=KRAS_Mut[!is.na(KRAS_Mut$top100OfftargetTotalScore),]
  }
  # remove sgRNAs with Day-0 READS less than the threshold
  from_mut=KRAS_Mut$BD0_G1<removesgRNA_Threshold | KRAS_Mut$BD0_G2<removesgRNA_Threshold
  from_WT=KRAS_WT2$BD0_G1<removesgRNA_Threshold | KRAS_WT2$BD0_G2<removesgRNA_Threshold
  ###Threshold of Raw counts
  from_mut=KRAS_Mut$BD0_G1<removesgRNA_Threshold | KRAS_Mut$BD0_G2<removesgRNA_Threshold
  from_WT=KRAS_WT2$BD0_G1<removesgRNA_Threshold | KRAS_WT2$BD0_G2<removesgRNA_Threshold
  KRAS_WT2=KRAS_WT2[!(from_WT | from_mut),]
  KRAS_Mut=KRAS_Mut[!(from_WT | from_mut),]
  
  # Using general CDE rank
  KRAS_Mut$CDE_Rank=NA
  KRAS_Mut$CDE_Rank[KRAS_Mut$Gene.group=='KRAS CDE Pos']= match(KRAS_Mut$Gene[KRAS_Mut$Gene.group=='KRAS CDE Pos'],
                                                                CDE$KRAS.CDE.Pos)
  KRAS_Mut$CDE_Rank[KRAS_Mut$Gene.group=='KRAS CDE Neg']= match(KRAS_Mut$Gene[KRAS_Mut$Gene.group=='KRAS CDE Neg'],
                                                                CDE$KRAS.CDE.Neg)
  KRAS_WT2=KRAS_WT2[!is.na(KRAS_Mut$CDE_Rank),]
  KRAS_Mut=KRAS_Mut[!is.na(KRAS_Mut$CDE_Rank),]
  
  # Here we take top 5% CDE positive and top 50 CDE negative hits
  # topX_neg=ceiling(length(unique(KRAS_Mut$Gene[KRAS_Mut$Gene.group=="KRAS CDE Neg"]))*TopXPercent)
  # topX_pos=ceiling(length(unique(KRAS_Mut$Gene[KRAS_Mut$Gene.group=="KRAS CDE Pos"]))*TopXPercent)
  
  if(include_offTarget){
    GOI = (KRAS_Mut$CDE_Rank < topX_neg &  KRAS_Mut$Gene.group == 'KRAS CDE Neg' &
             (is.na(KRAS_Mut$top100OfftargetTotalScore) |
                KRAS_Mut$top100OfftargetTotalScore<5.8 ) ) |
      (KRAS_Mut$CDE_Rank < topX_pos &  KRAS_Mut$Gene.group == 'KRAS CDE Pos' &
         (is.na(KRAS_Mut$top100OfftargetTotalScore) |
            KRAS_Mut$top100OfftargetTotalScore>2.9 ))
  } else{
    GOI = (KRAS_Mut$CDE_Rank < topX_neg &  KRAS_Mut$Gene.group == 'KRAS CDE Neg') |
      (KRAS_Mut$CDE_Rank < topX_pos &  KRAS_Mut$Gene.group == 'KRAS CDE Pos')
  }
  
  df2plot=rbind(data.frame(geneName=KRAS_WT2$Gene[GOI],
                           RANK=KRAS_WT2$FC_rank[GOI], 
                           type='WT',
                           group=KRAS_Mut$Gene.group[GOI]),
                data.frame(geneName=KRAS_WT2$Gene[GOI],
                           RANK=KRAS_Mut$FC_rank[GOI],
                           type='Mut',
                           group=KRAS_Mut$Gene.group[GOI]))
  p_value=c(wilcox.test(df2plot$RANK[df2plot$type=='WT' & df2plot$group=='KRAS CDE Pos'],
                        df2plot$RANK[df2plot$type=='Mut' & df2plot$group=='KRAS CDE Pos'],
                        alternative = 'l', paired=whether_matched)$p.value,
            wilcox.test(df2plot$RANK[df2plot$type=='WT' & df2plot$group=='KRAS CDE Neg'],
                        df2plot$RANK[df2plot$type=='Mut' & df2plot$group=='KRAS CDE Neg'],
                        alternative = 'g', paired=whether_matched)$p.value)
  
  # A=ggplot(df2plot, aes(x=group, fill=type, y=RANK))+geom_boxplot()+
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  #   annotate(geom="text", y=0, x=0.5, label=round(p_value[2], 2), color="blue")+
  #   annotate(geom="text", y=0, x=1.5, label=round(p_value[1], 2), color="blue")
  df2plot_genelevel=do.call(rbind, lapply(split(df2plot, list(as.character(df2plot$geneName), as.character(df2plot$type))), function(x) data.frame(geneName=unique(x$geneName),
                                                                                                                                                   RANK=mean(x$RANK),
                                                                                                                                                   group=unique(x$group),
                                                                                                                                                   type=unique(x$type)
  )))
  p_value_genelevel=c(wilcox.test(df2plot_genelevel$RANK[df2plot_genelevel$type=='WT' & df2plot_genelevel$group=='KRAS CDE Pos'],
                                  df2plot_genelevel$RANK[df2plot_genelevel$type=='Mut' & df2plot_genelevel$group=='KRAS CDE Pos'],
                                  alternative = 'l', paired=whether_matched)$p.value,
                      wilcox.test(df2plot_genelevel$RANK[df2plot_genelevel$type=='WT' & df2plot_genelevel$group=='KRAS CDE Neg'],
                                  df2plot_genelevel$RANK[df2plot_genelevel$type=='Mut' & df2plot_genelevel$group=='KRAS CDE Neg'],
                                  alternative = 'g', paired=whether_matched)$p.value)
  
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
# Functions neeeded during the Pre-processing
############################################################################
#Here we provide Top 
pooled_crisprCas9_kras=KRAS_ProcessingD0(removesgRNA_Threshold=150,
                                        topX_neg = 100,
                                        topX_pos = 10,
                                        whether_matched = T,
                                        CDE=CDE,
                                        only_DepMap=F,
                                        include_offTarget = F)
