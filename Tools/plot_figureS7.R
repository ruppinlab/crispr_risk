# Version for Github
source('Tools/Step0_Globally_used_Functions_and_Datasets.R')
require(PRROC);   # Library to generate and analyze PR curves
require(rafalib);   # This library is used to adjust labeling on graphs

# Cas9 inducing signature difference by KRAS status
############################################################
# Step 0:: Preprocessing data
############################################################
all_screens=readRDS('Data/seven_RPE_crisprcas9_screen.RDS')
names(all_screens)=c('p53 WT fg1_zimmerman',
                     'p53 WT fg2_unp',
                     'p53 WT fg3_unp', 
                     'p53 WT fg4_hart',
                     'p53 Null fg5_hart',
                     'p53 Null fg6_haapameini',
                     'p53 WT fg7_haapameini')
common_genes=Reduce(intersect, lapply(all_screens, names))
CDE=read.csv('Data/CDE_allThree.csv')
sapply(1:7, function(x) wilcox.test(vectorSubset(all_screens[[6]], CDE$p53.CDE.),
                                    vectorSubset(all_screens[[x]], CDE$p53.CDE.),
                                    alternative='l')$p.value)
all_screens_CDE_scores=sapply(all_screens, function(x) x[match(CDE$p53.CDE.[CDE$p53.CDE.!=''], names(x))] )
colnames(all_screens_CDE_scores)=c('p53 WT fg1_zimmerman',
                                   'p53 WT fg2_unp',
                                   'p53 WT fg3_unp', 
                                   'p53 WT fg4_hart',
                                   'p53 Null fg5_hart',
                                   'p53 Null fg6_haapameini',
                                   'p53 WT fg7_haapameini')
all_screens_CDE_scores_long=melt(all_screens_CDE_scores, id.vars=c('p53 WT fg1_zimmerman',
                                                                   'p53 WT fg2_unp',
                                                                   'p53 WT fg3_unp', 
                                                                   'p53 WT fg4_hart',
                                                                   'p53 Null fg5_hart',
                                                                   'p53 Null fg6_haapameini',
                                                                   'p53 WT fg7_haapameini'))
all_screens_CDE_scores_long$p53='WT'
all_screens_CDE_scores_long$p53[grep('Null',all_screens_CDE_scores_long$Var2)]='Mut'
all_screens_CDE_scores_long$Var2=factor(all_screens_CDE_scores_long$Var2)
all_screens_CDE_scores_long$Var2 = factor(all_screens_CDE_scores_long$Var2,
                                          levels(all_screens_CDE_scores_long$Var2)[
                                            order(aggregate(value ~ Var2,all_screens_CDE_scores_long, median)[,2],
                                                  decreasing = T)])
comparisons_needed=combn(levels(all_screens_CDE_scores_long$Var2), 2, simplify = F)
comparisons_needed=comparisons_needed[sapply(comparisons_needed, function(x) sum(grepl('WT',x))<2)]

tiff('Plots/figureS7.tiff')
ggplot(all_screens_CDE_scores_long, aes(x=reorder(Var2, value, FUN = median),
                                        y=value,
                                        fill=p53))+
  geom_boxplot()+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x='Screens', y='CDE+ genes essentiality')+
  stat_compare_means(method='wilcox',
                     label = 'p',
                     comparisons = comparisons_needed)
dev.off()