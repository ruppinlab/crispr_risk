# Version for Github
source('Step0_Globally_used_Functions_and_Datasets.R')
source('Step7B_analyze_CRISPR_screen_data_in_MOLM13_p53.R') #!!! since this script currently has problem, I have not tested the p53 plotting

# Here, we plto Figure 2B and 4G
########################################################################
##Pre-processing - TP53
########################################################################
df2plot_CRISPRcas9=pooled_crisprCas9_p53
df2plot_CRISPRi=pooled_crisprI_p53
df2plot_CRISPRcas9$`Gene Silencing Methods`='CRISPR-Cas9'
df2plot_CRISPRi$`Gene Silencing Methods`='CRISPR-interference'
df2plot_combined=rbind(df2plot_CRISPRi, df2plot_CRISPRcas9)
levels(df2plot_combined$type)=c('WildType', 'Mutated')
########################################################################
##Calculate Significances for the plots to follow
########################################################################
p_values=c(wilcox.test(df2plot_combined$RANK[df2plot_combined$group=='p53 CDE Pos' &
                                               df2plot_combined$`Gene Silencing Methods`=='CRISPR-Cas9'] ~
                         df2plot_combined$type[df2plot_combined$group=='p53 CDE Pos' &
                                                 df2plot_combined$`Gene Silencing Methods`=='CRISPR-Cas9'],
                       alternative='l', paired=T)$p.value,
           wilcox.test(df2plot_combined$RANK[df2plot_combined$group=='p53 CDE Neg' &
                                               df2plot_combined$`Gene Silencing Methods`=='CRISPR-Cas9'] ~
                         df2plot_combined$type[df2plot_combined$group=='p53 CDE Neg' &
                                                 df2plot_combined$`Gene Silencing Methods`=='CRISPR-Cas9'],
                       alternative='g', paired=T)$p.value,
           wilcox.test(df2plot_combined$RANK[df2plot_combined$group=='p53 CDE Pos' &
                                               df2plot_combined$`Gene Silencing Methods`=='CRISPR-interference'] ~
                         df2plot_combined$type[df2plot_combined$group=='p53 CDE Pos' &
                                                 df2plot_combined$`Gene Silencing Methods`=='CRISPR-interference'],
                       alternative='l', paired=T)$p.value,
           wilcox.test(df2plot_combined$RANK[df2plot_combined$group=='p53 CDE Neg' &
                                               df2plot_combined$`Gene Silencing Methods`=='CRISPR-interference'] ~
                         df2plot_combined$type[df2plot_combined$group=='p53 CDE Neg' &
                                                 df2plot_combined$`Gene Silencing Methods`=='CRISPR-interference'],
                       alternative='g', paired=T)$p.value
)

########################################################################
##Plotting sgRNA-Level
########################################################################
df2plot_combined=na.omit(df2plot_combined)
pvalues_df=data.frame(round(p_values, 3),
                      `Gene Silencing Methods`=rep(c('CRISPR-Cas9','CRISPR-interference'), each=2),
                      group=c(rep(c('p53 CDE Pos', 'p53 CDE Neg'), 2)),
                      ypos=c(4500),
                      xpos=rep(c('CRISPR-Cas9','CRISPR-interference'), each=2),
                      type=rep(c('WildType', 'Mutated'),2)
)

tiff('../Plots/isogenic_screen_p53.tiff')
ggplot(df2plot_combined, aes(x=`Gene Silencing Methods`, fill=type, y=RANK))+
  facet_wrap(group~., scales = 'free', nrow = 2)+
  geom_boxplot()+
  labs(fill = "Status of P53", y='sgRNA Essentiality Rank', x='Silencing Method Used')+
  theme_classic(base_size = 15)+
  geom_text(data=pvalues_df, aes(label=round(p_values, 3), x=xpos, y=ypos))
dev.off()

##OVERlapping Arguemnets
#**Plot KRAS figure : 4G

########################################################################
##PreProcessing
########################################################################
source('Step7C_analyze_CRISPR_screen_data_in_MOLM13_KRAS.R') #!!! originally here source the file for Step7B (p53), you mean KRAS? this part is supposed to be for KRAS?
df2plot=pooled_crisprCas9_kras
df2plot$`Gene Silencing Methods`='CRISPR-Cas9'
whether_matched=TRUE
p_value_sgRNAlevel=c(wilcox.test(df2plot$RANK[df2plot$type=='WT' & df2plot$group=='KRAS CDE Pos'],
                                 df2plot$RANK[df2plot$type=='Mut' & df2plot$group=='KRAS CDE Pos'],
                                 alternative = 'l', paired=whether_matched)$p.value,
                     wilcox.test(df2plot$RANK[df2plot$type=='WT' & df2plot$group=='KRAS CDE Neg'],
                                 df2plot$RANK[df2plot$type=='Mut' & df2plot$group=='KRAS CDE Neg'],
                                 alternative = 'g', paired=whether_matched)$p.value) #!!! whether_matched not found
levels(df2plot$type)=c('WildType', 'Mutated')
########################################################################
#At sgRNA level : 
########################################################################
pvalues_df=data.frame(p_values=p_value_sgRNAlevel,
                      `Gene Silencing Methods`=rep(c('CRISPR-Cas9'), each=2),
                      group=c(rep(c('KRAS CDE Pos', 'KRAS CDE Neg'), 1)),
                      ypos=c(5500),
                      xpos=rep(c('CRISPR-Cas9'), each=2),
                      type=rep(c('WildType', 'Mutated'),1)
)

tiff('../Plots/isogenic_screen_KRAS.tiff', height=300)
ggplot(df2plot, aes(x=`Gene Silencing Methods`, fill=type, y=RANK))+
  facet_wrap(group~., scales = 'free', nrow = 1)+
  geom_boxplot()+
  labs(fill = "Status of KRAS", y='Gene Essentiality Rank', x='Silencing Method Used')+
  theme_classic(base_size = 15)+
  geom_text(data=pvalues_df, aes(label=round(p_values, 3), x=xpos, y=ypos))
dev.off()

