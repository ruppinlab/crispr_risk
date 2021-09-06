####Centralized pipeline for All screenings
########################################################################
##Load all screenings data
########################################################################
require(ggpubr)
setwd('/Users/sinhas8/Project_CRISPR/2.Data/Isogenic_Screenings/')
CDE=read.csv('/Users/sinhas8/Project_CRISPR/2.Data/CDE_allThree.csv')
names(CDE)[1:2]=c('P53 CDE Pos', 'P53 CDE Neg')
CDE_bckp=CDE
########################################################################
#***Public Data***#
########################################################################
##KRAS shRNA
########################################################################
df_kras=read.csv('KRAS_isogenic_shRNA_screening.csv')
CDE=read.csv('/Users/sinhas8/Project_CRISPR/2.Data/CDE_allThree.csv')
df_kras=df_kras[!is.na(match(df_kras$Gene_Symbol, rownames(mat))),]
df_kras$GeneType=NA
df_kras$GeneType[!is.na(match(df_kras$Gene_Symbol, CDE$KRAS.CDE.))]='KRAS CDE Pos'
df_kras$GeneType[!is.na(match(df_kras$Gene_Symbol, CDE$KRAS.CDE..1))]='KRAS CDE Neg'
df_kras$GeneType[!(!is.na(match(df_kras$Gene_Symbol, CDE$KRAS.CDE.)) | !is.na(match(df_kras$Gene_Symbol, CDE$KRAS.CDE..1))) ]='KRAS CDE Neutral'
df_kras$Rank_WT=rank(df_kras$KRAS.WT.log2.Mean)
df_kras$Rank_Mut=rank(df_kras$KRAS.Mut.log2.mean)
test_direction_kras<-function(cde_type=levels(factor(df_kras$GeneType))[1], which_tail='t', Instead_rankMode=F){
  if(Instead_rankMode==F){
    FC_inMutant=df_kras$KRAS.Mut.log2.mean[df_kras$GeneType==cde_type]
    FC_inWT=df_kras$KRAS.WT.log2.Mean[df_kras$GeneType==cde_type]
  }else {
    FC_inMutant=df_kras$Rank_Mut[df_kras$GeneType==cde_type]
    FC_inWT=df_kras$Rank_WT[df_kras$GeneType==cde_type]
  }
  wilcox.test(FC_inMutant, FC_inWT, paired = T, alternative = which_tail)$p.value
}
##Calling the above function in rank mode
test_direction_kras(cde_type=levels(factor(df_kras$GeneType))[1], which_tail='l', Instead_rankMode=T)
test_direction_kras(cde_type=levels(factor(df_kras$GeneType))[1], which_tail='t', Instead_rankMode=T)
test_direction_kras(cde_type=levels(factor(df_kras$GeneType))[1], which_tail='g', Instead_rankMode=T)

################################################################
##P53 shRNA
################################################################
emsemblmapping=read.csv('mart_export.txt', sep='\t')
df_p53=read.csv('p53_isogenic_shRNA_screening.csv')
CDE=read.csv('/Users/sinhas8/Project_CRISPR/2.Data/CDE_allThree.csv')

###Preprocessing and testing hypothesis
df_p53$Gene_Symbol=emsemblmapping$Gene.name[match(df_p53$Ensemble, emsemblmapping$Gene.stable.ID)]
df_p53=na.omit(df_p53)
df_p53$GeneType=NA
df_p53$GeneType[!is.na(match(df_p53$Gene_Symbol, CDE$p53.CDE.))]='p53 CDE Pos'
df_p53$GeneType[!is.na(match(df_p53$Gene_Symbol, CDE$p53.CDE..1))]='p53 CDE Neg'
df_p53$GeneType[!(!is.na(match(df_p53$Gene_Symbol, CDE$p53.CDE.)) | !is.na(match(df_p53$Gene_Symbol, CDE$p53.CDE..1))) ]='p53 CDE Neutral'
table(df_p53$GeneType)

df_p53$Mut.log2.Mean=rowMeans(df_p53[,c("z..KO.A.","z..KO.B.")])
df_p53$WT.log2.Mean=rowMeans(df_p53[,c("z...WT.A.","z..WT.B.")])

df_p53$Rank_WT=rank(df_p53$WT.log2.Mean)
df_p53$Rank_Mut=rank(df_p53$Mut.log2.Mean)
test_direction_p53RNAi<-function(cde_type=levels(factor(df_p53$GeneType))[1], which_tail='t', Instead_rankMode=F){
  if(Instead_rankMode==F){
    FC_inMutant=df_p53$Mut.log2.Mean[df_p53$GeneType==cde_type]
    FC_inWT=df_p53$WT.log2.Mean[df_p53$GeneType==cde_type]
  }else {
    FC_inMutant=df_p53$Rank_Mut[df_p53$GeneType==cde_type]
    FC_inWT=df_p53$Rank_WT[df_p53$GeneType==cde_type]
  }
  wilcox.test(FC_inMutant, FC_inWT, paired = T, alternative = which_tail)$p.value
}
##Calling the above function in rank mode
test_direction_p53RNAi(cde_type=levels(factor(df_p53$GeneType))[1], which_tail='l', Instead_rankMode=T)
test_direction_p53RNAi(cde_type=levels(factor(df_p53$GeneType))[1], which_tail='t', Instead_rankMode=T)
test_direction_p53RNAi(cde_type=levels(factor(df_p53$GeneType))[1], which_tail='g', Instead_rankMode=T)


################################################################
##KRAS CRSIPR- DLD1
################################################################
##Load Files ##DLD1
CDE=read.csv('/Users/sinhas8/Project_CRISPR/2.Data/CDE_allThree.csv')
df_martin_kras=read.csv('/Users/sinhas8/Project_CRISPR/2.Data/Isogenic_Screenings/KRAS_Martin_HCT116.csv')

df_martin_kras$GeneType=NA
df_martin_kras$GeneType[!is.na(match(df_martin_kras$Gene, CDE$KRAS.CDE.))]='KRAS CDE Pos'
df_martin_kras$GeneType[!is.na(match(df_martin_kras$Gene, CDE$KRAS.CDE..1))]='KRAS CDE Neg'
df_martin_kras$GeneType[!(!is.na(match(df_martin_kras$Gene, CDE$KRAS.CDE.)) | !is.na(match(df_martin_kras$Gene, CDE$KRAS.CDE..1))) ]='KRAS CDE Neutral'

##Compute Ranks
##**Results were based on order func before
# df_martin_kras$Rank_WT=order(df_martin_kras$WT.log2.fold.change)
# df_martin_kras$Rank_Mut=order(df_martin_kras$KRas.mutant.log2.fold.change)
head(df_martin_kras)
df_martin_kras$Rank_WT=rank(df_martin_kras$WT.log2.fold.change)
df_martin_kras$Rank_Mut=rank(df_martin_kras$KRas.mutant.log2.fold.change)
##Test-Try1:: Method 1
colnames(df_martin_kras)
test_direction_krascrispr<-function(cde_type=levels(factor(df_martin_kras$GeneType))[1], which_tail='t', Instead_rankMode=F,
                                    df_martin_kras=df_martin){
  if(Instead_rankMode==F){
    FC_inMutant=df_martin_kras$WT.log2.fold.change[df_martin_kras$GeneType==cde_type]
    FC_inWT=df_martin_kras$KRas.mutant.log2.fold.change[df_martin_kras$GeneType==cde_type]
  }else {
    FC_inMutant=df_martin_kras$Rank_Mut[df_martin_kras$GeneType==cde_type]
    FC_inWT=df_martin_kras$Rank_WT[df_martin_kras$GeneType==cde_type]
  }
  wilcox.test(FC_inMutant, FC_inWT, paired = T, alternative = which_tail)$p.value
}
##

test_direction_krascrispr(cde_type=levels(factor(df_martin_kras$GeneType))[1], which_tail='l', Instead_rankMode=F, df_martin_kras=df_martin_kras)
test_direction_krascrispr(cde_type=levels(factor(df_martin_kras$GeneType))[2], which_tail='t', Instead_rankMode=F, df_martin_kras=df_martin_kras)
test_direction_krascrispr(cde_type=levels(factor(df_martin_kras$GeneType))[3], which_tail='g', Instead_rankMode=F, df_martin_kras=df_martin_kras)

#Rank mode off
test_direction_krascrispr(cde_type=levels(factor(df_martin$GeneType))[1], which_tail='l', Instead_rankMode=T, df_martin_kras=df_martin)
test_direction_krascrispr(cde_type=levels(factor(df_martin$GeneType))[2], which_tail='t', Instead_rankMode=T, df_martin_kras=df_martin)
test_direction_krascrispr(cde_type=levels(factor(df_martin$GeneType))[3], which_tail='g', Instead_rankMode=T, df_martin_kras=df_martin)



################################################################
##KRAS CRSIPR - HCT116
################################################################
df_martin=read.csv('/Users/sinhas8/Project_CRISPR/2.Data/Isogenic_Screenings/KRAS_Martin_HCT116.csv')
df_martin$Rank_WT=order(df_martin$WT.log2.fold.change)
df_martin$Rank_Mut=order(df_martin$KRas.mutant.log2.fold.change)
###
wilcox.test(df_martin$WT.log2.fold.change, df_martin$KRas.mutant.log2.fold.change, alternative='g', paired=T)
wilcox.test(df_martin_neg$WT.log2.fold.change, df_martin_neg$KRas.mutant.log2.fold.change, alternative='g', paired=T)
wilcox.test(df_martin_pos$WT.log2.fold.change, df_martin_pos$KRas.mutant.log2.fold.change, alternative='l', paired=T)
wilcox.test(df_martin_neutral$WT.log2.fold.change, df_martin_neutral$KRas.mutant.log2.fold.change, alternative='g', paired=T)


################################################################
##P53-CRSIPR--RPE1
################################################################
setwd('/Users/sinhas8/Project_CRISPR/Project_CRISPR_bias/Tools/')
source('Globally_used_Functions_and_Datasets.R')
Mut_R1=read.csv('/Users/sinhas8/Project_CRISPR/2.Data/Genome wide screening data RPE/R1-null.gene_summary.txt', sep='\t')
Mut_R2=read.csv('/Users/sinhas8/Project_CRISPR/2.Data/Genome wide screening data RPE/R2-null.gene_summary.txt', sep='\t')
WT_R1=read.csv('/Users/sinhas8/Project_CRISPR/2.Data/Genome wide screening data RPE/R1-WT.gene_summary.txt', sep='\t')
WT_R2=read.csv('/Users/sinhas8/Project_CRISPR/2.Data/Genome wide screening data RPE/R2-WT.gene_summary.txt', sep='\t')
CDE=read.csv('/Users/sinhas8/Project_CRISPR/2.Data/CDE_allThree.csv')
##Prep: Step 0
WT_R1=WT_R1[!is.na(match(as.character(unlist(WT_R1$id)), rownames(avana))), ]
WT_R2=WT_R2[!is.na(match(as.character(unlist(WT_R2$id)), rownames(avana))), ]
Mut_R1=Mut_R1[!is.na(match(as.character(unlist(Mut_R1$id)), rownames(avana))), ]
Mut_R2=Mut_R2[!is.na(match(as.character(unlist(Mut_R2$id)), rownames(avana))), ]
##Prep
Mut_R2=Mut_R2[match(Mut_R1$id, Mut_R2$id),]
Mut_avgLFC=data.frame(geneName=Mut_R2$id, avgLFC=rowMeans(data.frame(Mut_R2$neg.lfc, Mut_R1$neg.lfc )) )
WT_R2=WT_R2[match(WT_R1$id, WT_R2$id),]
WT_avgLFC=data.frame(geneName=WT_R2$id, avgLFC=rowMeans(data.frame(WT_R2$neg.lfc, WT_R1$neg.lfc)) )
WT_avgLFC=WT_avgLFC[match(Mut_avgLFC$geneName, WT_avgLFC$geneName),]
df_p53Haap=data.frame(WT_avgLFC, Mut=Mut_avgLFC)
#####*It was wrong before*
# df_p53Haap$Rank_WT=order(df_p53Haap$avgLFC)
# df_p53Haap$Rank_Mut=order(df_p53Haap$Mut.avgLFC)
df_p53Haap$Rank_WT=rank(df_p53Haap$avgLFC)
df_p53Haap$Rank_Mut=rank(df_p53Haap$Mut.avgLFC)
##Test Hypothesis:: Method 1: CDE medians comparison
df_p53Haap$GeneType=NA
df_p53Haap$GeneType[!is.na(match(df_p53Haap$geneName, CDE$p53.CDE.))]='p53 CDE Pos'
df_p53Haap$GeneType[!is.na(match(df_p53Haap$geneName, CDE$p53.CDE..1))]='p53 CDE Neg'
df_p53Haap$GeneType[!(!is.na(match(df_p53Haap$geneName, CDE$p53.CDE.)) | !is.na(match(df_p53Haap$geneName, CDE$p53.CDE..1))) ]='p53 CDE Neutral'
table(df_p53Haap$GeneType)
colnames(df_p53Haap)
test_direction_p53crispr<-function(cde_type=levels(factor(df_p53Haap$GeneType))[1], which_tail='t', Instead_rankMode=F){
  if(Instead_rankMode==F){
    FC_inMutant=df_p53Haap$avgLFC[df_p53Haap$GeneType==cde_type]
    FC_inWT=df_p53Haap$Mut.avgLFC[df_p53Haap$GeneType==cde_type]
  }else {
    FC_inMutant=df_p53Haap$Rank_Mut[df_p53Haap$GeneType==cde_type]
    FC_inWT=df_p53Haap$Rank_WT[df_p53Haap$GeneType==cde_type]
  }
  wilcox.test(FC_inMutant, FC_inWT, paired = T, alternative = which_tail)$p.value
}
##
levels(factor(df_p53Haap$GeneType))
test_direction_p53crispr(cde_type=levels(factor(df_p53Haap$GeneType))[2], which_tail='t', Instead_rankMode=F)
test_direction_p53crispr(cde_type=levels(factor(df_p53Haap$GeneType))[1], which_tail='l', Instead_rankMode=F)
test_direction_p53crispr(cde_type=levels(factor(df_p53Haap$GeneType))[3], which_tail='g', Instead_rankMode=F)
#Rank mode off
test_direction_p53crispr(cde_type=levels(factor(df_p53Haap$GeneType))[1], which_tail='l', Instead_rankMode=T)
test_direction_p53crispr(cde_type=levels(factor(df_p53Haap$GeneType))[2], which_tail='t', Instead_rankMode=T)
test_direction_p53crispr(cde_type=levels(factor(df_p53Haap$GeneType))[3], which_tail='g', Instead_rankMode=T)


################################################################
#***Our Data***#
################################################################
############################################################################
#P53 CRSIPR
############################################################################
##Test hypothesis
test_direction<-function(cde_type=levexls(p53Mutant$Gene.group)[1],
                         which_tail='t', Instead_rankMode=F, p53Mutant, p53WT){
  if(Instead_rankMode==F){
    FC_inMutant=p53Mutant$FC[p53Mutant$Gene.group==cde_type]
    FC_inWT=p53WT$FC[p53WT$Gene.group==cde_type]
  }else {
    FC_inMutant=p53Mutant$FC_rank[p53Mutant$Gene.group==cde_type]
    FC_inWT=p53WT$FC_rank[p53WT$Gene.group==cde_type]
  }
  wilcox.test(FC_inMutant, FC_inWT, paired = T, alternative = which_tail)$p.value
}

###Define Off-Target score for sgRNAs
setwd('/Users/sinhas8/Project_CRISPR/2.Data')
offTarg_v2=readRDS('results1_v2.RDS')
head(offTarg_v2$offtarget)
head(offTarg_v2$summary)
OT_top100scores=offTarg_v2$summary
head(OT_top100scores)
OT_top100scores$gRNA=substring(OT_top100scores$gRNAsPlusPAM, first=0,last=20)

P53_ProcessingD0<-function(removesgRNA_Threshold=20,  topX_neg=200, topX_pos=200,
                           whether_paired=T, CDE_bckp=MOLM, only_DepMap=F, include_offTarget=T){
  setwd('/Users/sinhas8/Project_CRISPR/2.Data/Isogenic_Screenings/')
  p53WT=read.csv('p53_CDEScreening_WildType.csv')
  p53Mutant=read.csv('p53_CDEScreening_Mutant.csv')
  colnames(p53Mutant)=colnames(p53WT)
  ##Step 0: Pre-Processing
  p53WT$avg_D0=rowMeans(p53WT[,c("R1_DO_cpm","R2_DO_cpm")])
  p53WT$avg_D30=rowMeans(p53WT[,c("R1_D3O_cpm","R2_D3O_cpm")])
  p53WT$FC= p53WT$avg_D30/p53WT$avg_D0
  p53WT$FC_rank=rank(p53WT$FC)
  p53Mutant$avg_D0=rowMeans(p53Mutant[,c("R1_DO_cpm","R2_DO_cpm")])
  p53Mutant$avg_D30=rowMeans(p53Mutant[,c("R1_D3O_cpm","R2_D3O_cpm")])
  p53Mutant$FC= p53Mutant$avg_D30/p53Mutant$avg_D0
  ##Define rank
  p53Mutant$FC_rank=rank(p53Mutant$FC)
  p53WT$FC_rank=rank(p53WT$FC)

  p53Mutant$top100OfftargetTotalScore=OT_top100scores$top100OfftargetTotalScore[match(p53Mutant$sgRNA, OT_top100scores$gRNA)]
  p53Mutant$top50OfftargetTotalScore=OT_top100scores$top5OfftargetTotalScore[match(p53Mutant$sgRNA, OT_top100scores$gRNA)]
  if(only_DepMap){
    p53WT=p53WT[!is.na(p53Mutant$top100OfftargetTotalScore),]
    p53Mutant=p53Mutant[!is.na(p53Mutant$top100OfftargetTotalScore),]
  }
  #CDE specific rank
  CDE=CDE_bckp
  CDE=lapply(CDE, function(x) x)
  CDE$`P53 CDE Pos`=CDE$`P53 CDE Pos`[!is.na(match(CDE$`P53 CDE Pos`,
                                                   p53Mutant$Gene[p53Mutant$Gene.group=='p53 CDE Pos']))]
  CDE$`P53 CDE Neg`=CDE$`P53 CDE Neg`[!is.na(match(CDE$`P53 CDE Neg`,
                                                   p53Mutant$Gene[p53Mutant$Gene.group=='p53 CDE Neg']))]

  #Day-0 READS threshold
  from_mut=p53Mutant$R1_DO_Count<removesgRNA_Threshold | p53Mutant$R2_DO_Count<removesgRNA_Threshold
  from_WT=p53WT$R1_DO_Count<removesgRNA_Threshold | p53WT$R2_DO_Count<removesgRNA_Threshold
  p53WT=p53WT[!(from_WT | from_mut),]
  p53Mutant=p53Mutant[!(from_WT | from_mut),]

  ##Define CDE Rank
  p53Mutant$CDE_Rank=NA
  p53Mutant$CDE_Rank[p53Mutant$Gene.group=='p53 CDE Pos']= match(p53Mutant$Gene[p53Mutant$Gene.group=='p53 CDE Pos'],
                                                                CDE$`P53 CDE Pos`)
  p53Mutant$CDE_Rank[p53Mutant$Gene.group=='p53 CDE Neg']= match(p53Mutant$Gene[p53Mutant$Gene.group=='p53 CDE Neg'],
                                                                CDE$`P53 CDE Neg`)

  #Top X Genes
  if(include_offTarget){
    one_third_point=summary(p53Mutant$top100OfftargetTotalScore)[2]
    two_third_point=summary(p53Mutant$top100OfftargetTotalScore)[3]
    table(p53Mutant$Gene.group)
    GOI = (p53Mutant$CDE_Rank < topX_neg &  p53Mutant$Gene.group == 'p53 CDE Neg' & (is.na(p53Mutant$top100OfftargetTotalScore) | p53Mutant$top100OfftargetTotalScore<one_third_point ) ) |
      (p53Mutant$CDE_Rank < topX_pos &  p53Mutant$Gene.group == 'p53 CDE Pos' & (is.na(p53Mutant$top100OfftargetTotalScore) | p53Mutant$top100OfftargetTotalScore>two_third_point ))
  } else{
    GOI = (p53Mutant$CDE_Rank < topX_neg &  p53Mutant$Gene.group == 'p53 CDE Neg') | (p53Mutant$CDE_Rank < topX_pos &  p53Mutant$Gene.group == 'p53 CDE Pos')
  }

  df2plot=rbind(data.frame(geneName=p53WT$Gene[GOI], RANK=p53WT$FC_rank[GOI], type='WT', group=p53Mutant$Gene.group[GOI]),
                data.frame(geneName=p53WT$Gene[GOI], RANK=p53Mutant$FC_rank[GOI], type='Mut', group=p53Mutant$Gene.group[GOI]))
  p_value=c(wilcox.test(df2plot$RANK[df2plot$type=='WT' & df2plot$group=='p53 CDE Pos'],
                        df2plot$RANK[df2plot$type=='Mut' & df2plot$group=='p53 CDE Pos'],
                        alternative = 'l', paired=whether_paired)$p.value,
            wilcox.test(df2plot$RANK[df2plot$type=='WT' & df2plot$group=='p53 CDE Neg'],
                        df2plot$RANK[df2plot$type=='Mut' & df2plot$group=='p53 CDE Neg'],
                        alternative = 'g', paired=whether_paired)$p.value)

  A=ggplot(df2plot, aes(x=group, fill=type, y=RANK))+geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    annotate(geom="text", y=0, x=0.5, label=round(p_value[2], 2), color="blue")+
    annotate(geom="text", y=0, x=1.5, label=round(p_value[1], 2), color="blue")
  df2plot_genelevel=do.call(rbind, lapply(split(df2plot, list(as.character(df2plot$geneName), as.character(df2plot$type))), function(x) data.frame(geneName=unique(x$geneName),
                                                                                                                                                   RANK=mean(x$RANK),
                                                                                                                                                   group=unique(x$group),
                                                                                                                                                   type=unique(x$type)
  )))
  p_value_genelevel=c(wilcox.test(df2plot_genelevel$RANK[df2plot_genelevel$type=='WT' & df2plot_genelevel$group=='p53 CDE Pos'],
                                  df2plot_genelevel$RANK[df2plot_genelevel$type=='Mut' & df2plot_genelevel$group=='p53 CDE Pos'],
                                  alternative = 'l', paired=whether_paired)$p.value,
                      wilcox.test(df2plot_genelevel$RANK[df2plot_genelevel$type=='WT' & df2plot_genelevel$group=='p53 CDE Neg'],
                                  df2plot_genelevel$RANK[df2plot_genelevel$type=='Mut' & df2plot_genelevel$group=='p53 CDE Neg'],
                                  alternative = 'g', paired=whether_paired)$p.value)

  B=ggplot(df2plot_genelevel, aes(x=group, fill=type, y=RANK))+geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    annotate(geom="text", y=0, x=0.5, label=round(p_value_genelevel[2], 2), color="blue")+
    annotate(geom="text", y=0, x=1.5, label=round(p_value_genelevel[1], 2), color="blue")
  grid.arrange(A, B, nrow = 2)
  c(p_value, p_value_genelevel)
}

##Optimum
P53_ProcessingD0(removesgRNA_Threshold=20,  topX_neg=14, topX_pos=19, whether_paired = T, CDE_bckp = MOLM, only_DepMap = T, include_offTarget = T)

############################################################################
#P53 - CRSIPRi
############################################################################
ip53=read.csv('/Users/sinhas8/Project_CRISPR/2.Data/Isogenic_Screenings/iP53_CDETopGenes.txt', sep='\t')
head(ip53)
############################################################################
#KRAS-CRSIRP -- Without any scaling by NTC
############################################################################
KRAS_ProcessingD0<-function(removesgRNA_Threshold=0, topX_neg=200, topX_pos=200, whether_matched=T, CDE_bckp, only_DepMap=T, include_offTarget=F){
  setwd('/Users/sinhas8/Project_CRISPR/2.Data/Isogenic_Screenings/')
  KRAS_Mut=read.csv('KRAS_Mut.csv')
  KRAS_WT2=read.csv('KRAS_WT2.csv')
  KRAS_Mut$FC = KRAS_Mut$Avg.D30/ KRAS_Mut$Avg.D0
  KRAS_WT2$FC = KRAS_WT2$Avg.D30/KRAS_WT2$Avg.D0
  colnames(KRAS_WT2)=colnames(KRAS_Mut)
  CDE=CDE_bckp
  CDE=lapply(CDE, function(x) x)
  names(CDE)[5:6]=c('KRAS CDE Pos', 'KRAS CDE Neg')
  CDE$`KRAS CDE Pos`=CDE$`KRAS CDE Pos`[!is.na(match(CDE$`KRAS CDE Pos`,
                                                     KRAS_Mut$Gene[KRAS_Mut$Gene.group=='KRAS CDE Pos']))]
  CDE$`KRAS CDE Neg`=CDE$`KRAS CDE Neg`[!is.na(match(CDE$`KRAS CDE Neg`,
                                                     KRAS_Mut$Gene[KRAS_Mut$Gene.group=='KRAS CDE Neg']))]
  KRAS_WT2$FC_rank=rank(KRAS_WT2$FC)
  KRAS_Mut$FC_rank=rank(KRAS_Mut$FC)

  KRAS_Mut$top100OfftargetTotalScore=OT_top100scores$top100OfftargetTotalScore[match(KRAS_Mut$sgRNA, OT_top100scores$gRNA)]
  KRAS_Mut$top50OfftargetTotalScore=OT_top100scores$top5OfftargetTotalScore[match(KRAS_Mut$sgRNA, OT_top100scores$gRNA)]
  dim(KRAS_Mut); dim(KRAS_WT2)
  if(only_DepMap){
    KRAS_WT2=KRAS_WT2[!is.na(KRAS_Mut$top100OfftargetTotalScore),]
    KRAS_Mut=KRAS_Mut[!is.na(KRAS_Mut$top100OfftargetTotalScore),]
  }

  ####Provide ranking

  ###Threshold of Raw counts
  # head(KRAS_Mut)
  from_mut=KRAS_Mut$BD0_G1<removesgRNA_Threshold | KRAS_Mut$BD0_G2<removesgRNA_Threshold
  from_WT=KRAS_WT2$BD0_G1<removesgRNA_Threshold | KRAS_WT2$BD0_G2<removesgRNA_Threshold
  ###Threshold of Raw counts
  from_mut=KRAS_Mut$BD0_G1<removesgRNA_Threshold | KRAS_Mut$BD0_G2<removesgRNA_Threshold
  from_WT=KRAS_WT2$BD0_G1<removesgRNA_Threshold | KRAS_WT2$BD0_G2<removesgRNA_Threshold
  KRAS_WT2=KRAS_WT2[!(from_WT | from_mut),]
  KRAS_Mut=KRAS_Mut[!(from_WT | from_mut),]


  #hist(KRAS_WT2$FC_rank)
  KRAS_Mut$CDE_Rank=NA
  KRAS_Mut$CDE_Rank[KRAS_Mut$Gene.group=='KRAS CDE Pos']= match(KRAS_Mut$Gene[KRAS_Mut$Gene.group=='KRAS CDE Pos'], CDE$`KRAS CDE Pos`)
  KRAS_Mut$CDE_Rank[KRAS_Mut$Gene.group=='KRAS CDE Neg']= match(KRAS_Mut$Gene[KRAS_Mut$Gene.group=='KRAS CDE Neg'], CDE$`KRAS CDE Neg`)
  KRAS_WT2=KRAS_WT2[!is.na(KRAS_Mut$CDE_Rank),]
  KRAS_Mut=KRAS_Mut[!is.na(KRAS_Mut$CDE_Rank),]

  if(include_offTarget){
  GOI = (KRAS_Mut$CDE_Rank < topX_neg &  KRAS_Mut$Gene.group == 'KRAS CDE Neg' & (is.na(KRAS_Mut$top100OfftargetTotalScore) | KRAS_Mut$top100OfftargetTotalScore<5.8 ) ) |
    (KRAS_Mut$CDE_Rank < topX_pos &  KRAS_Mut$Gene.group == 'KRAS CDE Pos' & (is.na(KRAS_Mut$top100OfftargetTotalScore) | KRAS_Mut$top100OfftargetTotalScore>2.9 ))
  } else{
    GOI = (KRAS_Mut$CDE_Rank < topX_neg &  KRAS_Mut$Gene.group == 'KRAS CDE Neg') | (KRAS_Mut$CDE_Rank < topX_pos &  KRAS_Mut$Gene.group == 'KRAS CDE Pos')
  }

  df2plot=rbind(data.frame(geneName=KRAS_WT2$Gene[GOI], RANK=KRAS_WT2$FC_rank[GOI], type='WT', group=KRAS_Mut$Gene.group[GOI]),
                data.frame(geneName=KRAS_WT2$Gene[GOI], RANK=KRAS_Mut$FC_rank[GOI], type='Mut', group=KRAS_Mut$Gene.group[GOI]))
  p_value=c(wilcox.test(df2plot$RANK[df2plot$type=='WT' & df2plot$group=='KRAS CDE Pos'],
                        df2plot$RANK[df2plot$type=='Mut' & df2plot$group=='KRAS CDE Pos'],
                        alternative = 'l', paired=whether_matched)$p.value,
            wilcox.test(df2plot$RANK[df2plot$type=='WT' & df2plot$group=='KRAS CDE Neg'],
                        df2plot$RANK[df2plot$type=='Mut' & df2plot$group=='KRAS CDE Neg'],
                        alternative = 'g', paired=whether_matched)$p.value)

  A=ggplot(df2plot, aes(x=group, fill=type, y=RANK))+geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    annotate(geom="text", y=0, x=0.5, label=round(p_value[2], 2), color="blue")+
    annotate(geom="text", y=0, x=1.5, label=round(p_value[1], 2), color="blue")
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

  B=ggplot(df2plot_genelevel, aes(x=group, fill=type, y=RANK))+geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    annotate(geom="text", y=0, x=0.5, label=round(p_value_genelevel[2], 2), color="blue")+
    annotate(geom="text", y=0, x=1.5, label=round(p_value_genelevel[1], 2), color="blue")
  grid.arrange(A, B, nrow = 2)
}

############################################################################
#Final Settings
############################################################################
P53_ProcessingD0(removesgRNA_Threshold=20,  topX_neg=14, topX_pos=19, whether_paired = T, CDE_bckp = MOLM, only_DepMap = T, include_offTarget = T)
KRAS_ProcessingD0(removesgRNA_Threshold=150, topX_neg = 100, topX_pos = 10, whether_matched = T, CDE_bckp=CDE_bckp, only_DepMap=F, include_offTarget = F)

P53_ProcessingD0(removesgRNA_Threshold=20,  topX_neg=200, topX_pos=200, whether_paired = T, CDE_bckp = MOLM, only_DepMap = F, include_offTarget = F)
iP53_ProcessingD0(removesgRNA_Threshold=20,  topX_neg=200, topX_pos=200, whether_paired = T, CDE_bckp = MOLM)

