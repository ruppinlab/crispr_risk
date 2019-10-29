# Version for Github
source('Step0_Globally_used_Functions_and_Datasets.R')

## Plotting Panel4_chart_hist_4th_aspectratio.pdf
FANC_Genes=as.character(read.csv('../Data/FANC_geneset.txt', header=T, skip=1)[,1])
CEll_Cycle_Genes=as.character(read.csv('../Data/Cell-Cycle_geneset.txt')[,1])
CDE_Pos_Neg=readRDS('../Data/CDE_for_volg.RDS')$TP53

#Similar plots could be created for other master Regulators
master_regulators='TP53'
median_diff_avana=apply(avana, 1, function(x)
  median(scale(x)[unlist(as.logical(Mut_CCLE[master_regulators,]))], na.rm=T) -
    median(scale(x)[unlist(!as.logical(Mut_CCLE[master_regulators,]))], na.rm=T))

##
sig=rep(1, nrow(avana))
sig[match(CDE_Pos_Neg[[1]], rownames(avana))]=2
sig[match(CDE_Pos_Neg[[2]],rownames(avana))]=3
sig=as.factor(sig)


##labelling
Gene=rep('',nrow(avana))
Top_Pos_to_mention = CDE_Pos_Neg[[1]][na.omit(match(FANC_Genes, CDE_Pos_Neg[[1]]))]
Gene[na.omit(match(Top_Pos_to_mention, rownames(avana)))]=rownames(avana)[na.omit(match(Top_Pos_to_mention, rownames(avana)))]

Top_Neg_to_mention = CDE_Pos_Neg[[2]][na.omit(match(CEll_Cycle_Genes, CDE_Pos_Neg[[2]]))]
Gene[match(Top_Neg_to_mention,rownames(avana))]=rownames(avana)[match(Top_Neg_to_mention,rownames(avana))]

##
levels(sig)=c('Non-Significant Risk','CDE+', 'CDE-')
test=data.frame(Risk_Score=median_diff_avana, sig, Gene)
test=test[order(test$Risk_Score, decreasing=T),]
test$x=seq(length(sig))
test1=test[test$Gene != '',]

##Plotting the risk plot
tt=ggplot(data=test, aes(x=x, y=Risk_Score, fill=sig))+ 
  geom_bar(stat='identity', aes(fill=sig))+
  labs(x='Genes at Risk', y="Risk Score")+theme_classic()+
  theme(axis.text=element_text(size=10), axis.title=element_text(size=15),
        legend.text=element_text(size=7),
        legend.title=element_text(size=0),
        legend.position="top", plot.margin=margin(10,10,30,10))+
  geom_text_repel(
    data = subset(test, Gene != ''),
    aes(label = Gene),
    size = 2.5,
    #	    box.padding = unit(0.35, "lines"),
    #			    point.padding = unit(0.3, "lines"),
    min.segment.length = 0,nudge_y=ifelse(test1$Risk_Score>0, 0.5, -0.5), 
    nudge_x=ifelse(test1$Risk_Score>0, 5000, -5000), segment.alpha=0.2, segment.size=0.25
  )+
  scale_fill_manual(values=c('grey','red','blue'))

##
tiff('../Plots/figure1D.tiff')
plot(tt)
dev.off()

