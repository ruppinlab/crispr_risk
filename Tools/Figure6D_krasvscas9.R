# Scatter plot for showing the rankca9 vs KRAS
##################################################################
# Test stable cas9 extent in WT vs Mut for each of the cancer gene
##################################################################
# Test cas9 activity by KRAS
CCLE_cellLineName=sapply(as.character(onTarget$annotation$CCLE_ID),
                         function(x) strsplit(x, '_')[[1]][1])
matched_mutmatrix=onTarget$mutations_matrix[,match(onTarget$annotation$depMapID, colnames(onTarget$mutations_matrix))]
load('/Users/sinhas8/Project_CRISPR/Project_CRISPR.Risk_Github/Data/vogelstein.RData')
cas9Activity=readxl::read_xlsx('/Users/sinhas8/Project_CRISPR/cas9_Activity_inCellLines.xlsx')

cas9Activity=cbind(cas9Activity, 
                   t(matched_mutmatrix[match(vogelstein$symbol, rownames(matched_mutmatrix)),
                                       match(cas9Activity$`Official Cell Line Name`,
                                             CCLE_cellLineName)]))
myhead(cas9Activity)

Volg_DifferentialCas9=t(apply(cas9Activity[,-(1:3)], 2, 
                              function(x) c(p_value=err_handle(wilcox.test(as.numeric(cas9Activity$`%GFP positive-cells`) ~ x,
                                                                           alternative='g')$p.value),
                                            Effect_size=err_handle(median(100-as.numeric(cas9Activity$`%GFP positive-cells`)[which(x==1)], na.rm=T)-
                                                                     median(100-as.numeric(cas9Activity$`%GFP positive-cells`)[which(x==0)], na.rm=T)))))
summary(lm(as.numeric(cas9Activity$`%GFP positive-cells`) ~ cas9Activity$KRAS+cas9Activity$TP53))$coefficients
median(as.numeric(cas9Activity$`%GFP positive-cells`)[cas9Activity$`TP53 status`=='mut'])
median(as.numeric(cas9Activity$`%GFP positive-cells`)[cas9Activity$`TP53 status`=='WT'])
##################################################################
# Plot cas9 vs KRAS
##################################################################
install.packages('tidyverse')
library(cowplot)
cas9Activity=readxl::read_xlsx('/Users/sinhas8/Project_CRISPR/cas9_Activity_inCellLines.xlsx')
cas9Activity=na.omit(cas9Activity)

pdf('/Users/sinhas8/Project_CRISPR/cas9byKRAS_figure4xx.pdf')
ggplot(cas9Activity, aes(x=factor(KRAS, labels = c('WT', 'Mutated')),
                         y=100-as.numeric(cas9Activity$`%GFP positive-cells`),
                         fill=factor(KRAS)))+
  geom_boxplot()+
  theme_classic(base_size = 20)+
  labs(x='KRAS Status', y='Stable Cas9 Activity')+
  theme(legend.position = 'none')+
  scale_fill_manual(values = c('grey', 'indianred2'))
dev.off()

table(cas9Activity$KRAS)