########################################################################
#######FIGURE_3B
########################################################################
library(tidyr)
library(Rmisc)
require(ggplot2)
require(ggpubr)
compAssay_Page2=readxl::read_xlsx('/Users/sinhas8/Downloads/Competition_Analysis3.xlsx', sheet = 2)
long_df2plot=gather(compAssay_Page2, Day, Percent_Mutant, D0:D25, factor_key=TRUE)

long_df2plot=long_df2plot[long_df2plot$Gene %in% levels(factor(long_df2plot$Gene))[3:5],]
colnames(long_df2plot)
# long_df2plot=do.call(rbind, lapply(split(long_df2plot, list(long_df2plot$System, long_df2plot$Gene, 
#                          long_df2plot$Day, long_df2plot$Mix)),
#        function(x) err_handle(data.frame(x[,c(1,2,3,4,9)], Percent_Mutant=x$Percent_Mutant))))
# head(long_df2plot)
summarydf2plot = summarySE(long_df2plot, measurevar="Percent_Mutant", groupvars=c("System","Mix", "Day", 'Gene'))


########################################################################
#######Improve Signal of sgRNAs by ranking based on Large screenings
########################################################################
sgRNAseq=read.csv('/Users/sinhas8/Downloads/top_koCDE_P53CompetitionAssay[15512] (1).txt', sep='\t', header=F)
###This is a incomplete thread of code from "Untitled 1-P53 function"
compAssay_Page1$Support_Score=NA
sgRNAseq$sgRNA_id=sapply(sgRNAseq$V1, function(x) strsplit(as.character(x), '-')[[1]][2])
dim(compAssay_Page1)
temp1=lapply(1:nrow(sgRNAseq), function(x) compAssay_Page1[which(compAssay_Page1$Gene==sgRNAseq$GeneName[x] & compAssay_Page1$sgRNA==sgRNAseq$sgRNA_id[x]),])
sgRNAseq[sgRNAseq$GeneName=='NTC',]
temp1_df=do.call(rbind, lapply(1:length(temp1), function(x) err_handle(data.frame(temp1[[x]],sgRNAseq[x,])) ))
NTC_df=cbind(compAssay_Page1[compAssay_Page1$Gene=='NTC',], 
      sgRNAseq[sgRNAseq$GeneName=='NTC',])
colnames(temp1_df)=colnames(NTC_df)
df_CRISPRko=rbind(temp1_df, NTC_df)
#colnames(df_CRISPRko)
dim(compAssay_Page3)
colnames(compAssay_Page3)=colnames(df_CRISPRko)[1:ncol(compAssay_Page3)]
compAssay_Page3=compAssay_Page3[,-ncol(compAssay_Page3)]
colnames(compAssay_Page3)

compAssay_Page3[,colnames(df_CRISPRko)[ncol(compAssay_Page3)+1]]=NA
manual_df=rbind(compAssay_Page3, df_CRISPRko)
df2plot=rbind(data.frame(geneName=p53WT$Gene[GOI], sgRNA=p53WT$sgRNA[GOI], RANK=p53WT$FC_rank[GOI], 
                         type='WT', group=p53Mutant$Gene.group[GOI], avg_D0=p53Mutant$avg_D0[GOI]),
              data.frame(geneName=p53WT$Gene[GOI], sgRNA=p53Mutant$sgRNA[GOI], RANK=p53Mutant$FC_rank[GOI],
                         type='Mut', group=p53Mutant$Gene.group[GOI], avg_D0=p53WT$avg_D0[GOI]))
sgRNAseq$GeneName=sapply(sgRNAseq$V1, function(x) strsplit(as.character(x), '-')[[1]][1])
map_exp12exp2 = lapply(1:nrow(sgRNAseq), function(x) df2plot[which(df2plot$geneName == sgRNAseq$GeneName[x] & 
                                                                   df2plot$sgRNA== as.character(unlist(sgRNAseq$V2[x])) &
                                                                   df2plot$avg_D0 !=  0), ])
map_exp12exp2=map_exp12exp2[sapply(map_exp12exp2, nrow) != 0]
bigExperiemnt_info=do.call(rbind, lapply(map_exp12exp2, function(x) data.frame(unique(x[,2]), diff(x[,3])) ))
sgRNAseq$support = bigExperiemnt_info$diff.x...3..[match(sgRNAseq$V2, bigExperiemnt_info$unique.x...2..)]

########################################################################
##Plot Mixture 1 ###Mixture ratio:: 5/95
########################################################################
###Try1
summarydf2plot$Gene[summarydf2plot$Gene=='NTC']='Control'
summarydf2plot$Day_num= as.numeric(sapply(summarydf2plot$Day, function(x) substring(as.character(x), 2,)))
long_df2plot$Day_num= as.numeric(sapply(long_df2plot$Day, function(x) substring(as.character(x), 2,)))
long_df2plot$Gene[long_df2plot$Gene=='NTC']='Control'
long_df2plot$Day_num=factor(long_df2plot$Day_num)
summarydf2plot$Day_num=factor(summarydf2plot$Day_num)

summarydf2plot$Gene[summarydf2plot$Gene=='Control']='Non-Targeting Control'
long_df2plot$Gene[long_df2plot$Gene=='Control']='Non-Targeting Control'

summarydf2plot$System=factor(summarydf2plot$System, labels = c('CRISPRi', 'CRSIPR-KO'))
long_df2plot$System=factor(long_df2plot$System, labels = c('CRISPRi', 'CRSIPR-KO'))
pdf('/Users/sinhas8/Project_CRISPR/figure3B.pdf', width=3, height=6)
A<-ggplot(summarydf2plot[summarydf2plot$Mix=='5/95',],
       aes(y=Percent_Mutant, x=Day, color=System))+
  geom_point()+
  geom_line(data = summarydf2plot[summarydf2plot$Mix=='5/95',],
            aes(x=as.numeric(Day), y=Percent_Mutant), size=1.3)+
  facet_wrap(Gene~., , nrow=3, scales='free')+
  geom_errorbar(aes(ymin=Percent_Mutant-se, ymax=Percent_Mutant+se),
                 width=.1)+
   theme_classic(base_size = 12)+
   theme(, legend.position = "top")+
  stat_compare_means(data=long_df2plot[long_df2plot$Mix=='50/50',], aes(group = System), 
                     label='p.signif', size=4, label.y = 30)+
  labs(x='Days', y='% p53-mutated cells in co-culture', color='')
dev.off()
ggsave(plot=A, file='/Users/sinhas8/Project_CRISPR/figure3B.pdf', width=3, height=6)
########################################################################
##Try2
########################################################################
Mut=summarydf2plot
WT=summarydf2plot
WT$Percent_Mutant= 100-WT$Percent_Mutant
comb_df=rbind(data.frame(Mut, CellType="Mutant"), 
              data.frame(WT, CellType="WT"))

ggplot(comb_df[comb_df$Mix=='5/95',], aes(y=Percent_Mutant, x=Day, fill=CellType))+
  facet_grid(Gene~System+Mix, scales = 'fixed')+
  geom_bar(stat='identity')+
  # geom_errorbar(aes(ymin=Percent_Mutant-se, ymax=Percent_Mutant+se),
  #               width=.2,                    # Width of the error bars
  #               position=position_dodge(.9))+
  theme_classic(base_size = 15)+
  theme(axis.text.x = element_text(angle = 90))

########################################################################
##Try3:: Manual DF
########################################################################
saveRDS(manual_df,'/Users/sinhas8/Project_CRISPR/2.Data/manual_df.RDS')
manual_df=readRDS('/Users/sinhas8/Project_CRISPR/2.Data/manual_df.RDS')
manual_df$Gene
long_manual_df=gather(manual_df, Day, Percent_Mutant, D5:D25, factor_key=TRUE)
long_manual_df=long_manual_df[long_manual_df$Gene %in% levels(factor(long_manual_df$Gene))[3:5],]
head(long_manual_df)

long_manual_df=do.call(rbind, lapply(split(long_manual_df, list(long_manual_df$System, long_manual_df$Gene, 
                                                                    long_manual_df$Day, long_manual_df$Mix)),
                                           function(x) err_handle(data.frame(x[,c(1,2,3,4,5,9,12, 14)], Percent_Mutant=x$Percent_Mutant ))
))

long_manual_df=long_manual_df[-which(long_manual_df$Gene=='NDUFB10' & long_manual_df$sgRNA=='4'),]


summarymanual_df = summarySE(long_manual_df, measurevar="Percent_Mutant", groupvars=c("System","Mix", "Day", 'Gene'))

Mut=summarydf2plot
WT=summarydf2plot
WT$Percent_Mutant= 100-WT$Percent_Mutant
comb_df=rbind(data.frame(Mut, CellType="Mutant"),
              data.frame(WT, CellType="WT"))
ggplot(comb_df[comb_df$Mix=='5/95',], aes(y=Percent_Mutant, x=Day, fill=CellType))+
  facet_grid(Gene~System+Mix, scales = 'fixed')+
  geom_bar(stat='identity')+
  # geom_errorbar(aes(ymin=Percent_Mutant-se, ymax=Percent_Mutant+se),
  #               width=.2,                    # Width of the error bars
  #               position=position_dodge(.9))+
  theme_classic(base_size = 15)+
  theme(axis.text.x = element_text(angle = 90))


########################################################################
#######Supp Figure S2
########################################################################
library(tidyr)
library(Rmisc)
compAssay_Page2=readxl::read_xlsx('/Users/sinhas8/Downloads/Competition_Analysis3.xlsx', sheet = 2)
head(compAssay_Page2)
compAssay_Page2[compAssay_Page2$Gene=='NDUFB10',]
long_df2plot=gather(compAssay_Page2, Day, Percent_Mutant, D0:D25, factor_key=TRUE)
long_df2plot=do.call(rbind, lapply(split(long_df2plot, list(long_df2plot$System, long_df2plot$Gene, 
                                                            long_df2plot$Day, long_df2plot$Mix)),
                                   function(x) err_handle(data.frame(x[,c(1,2,3,4,9,10)], Percent_Mutant=x$Percent_Mutant ))
))
summarydf2plot = summarySE(long_df2plot, measurevar="Percent_Mutant", groupvars=c("System","Mix", "Day", 'Gene'))
summarydf2plot_bckp=summarydf2plot
summarydf2plot=summarydf2plot[summarydf2plot$Gene=='ANKRD49' | summarydf2plot$Gene=='FANCG'| summarydf2plot$Gene=='TFAP4',]
summarydf2plot$Gene=factor(as.character(summarydf2plot$Gene))

tiff('/Users/sinhas8/Project_CRISPR/figureS2.tiff', width=800, height=300)
ggplot(summarydf2plot[summarydf2plot$Mix=='5/95',], aes(y=Percent_Mutant, x=Day, color=System))+
  geom_point()+
  geom_line(data = summarydf2plot[summarydf2plot$Mix=='5/95',],
            aes(x=as.numeric(Day), y=Percent_Mutant), size=1.3)+
  facet_wrap(Gene~., , nrow=1, scales='free')+
  geom_errorbar(aes(ymin=Percent_Mutant-se, ymax=Percent_Mutant+se),
                width=.1)+
  theme_classic(base_size = 20)+
  theme(, legend.position = "top")+
  labs(x='Days', y='% p53-mutated cells \n in co-culture', color='')+
  lims(y=c(0,50))
dev.off()
