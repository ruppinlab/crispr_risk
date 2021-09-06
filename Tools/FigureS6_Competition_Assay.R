# Version for Draft May 10
source('Tools/Step0_Globally_used_Functions_and_Datasets.R')

########################################################################
# Plot Supp Figure S2
########################################################################
compAssay_Page2=readxl::read_xlsx('Data/Competition_Analysis3.xlsx', sheet = 2)
long_df2plot=gather(compAssay_Page2, Day, Percent_Mutant, D0:D25, factor_key=TRUE)
long_df2plot=as.data.frame(long_df2plot)
summarydf2plot = summarySE(long_df2plot, measurevar="Percent_Mutant", groupvars=c("System","Mix", "Day", 'Gene'))
summarydf2plot=summarydf2plot[summarydf2plot$Gene=='ANKRD49' | summarydf2plot$Gene=='FANCG'| summarydf2plot$Gene=='TFAP4',]
summarydf2plot$Gene=factor(as.character(summarydf2plot$Gene))

tiff('Plots/figureS2.tiff', width=800, height=300)
ggplot(summarydf2plot[summarydf2plot$Mix=='5/95',], aes(y=Percent_Mutant, x=Day, color=System))+
  geom_point()+
  geom_line(data = summarydf2plot[summarydf2plot$Mix=='5/95',],
            aes(x=as.numeric(Day), y=Percent_Mutant), size=1.3)+
  facet_wrap(Gene~.,nrow=1, scales='free')+
  geom_errorbar(aes(ymin=Percent_Mutant-se, ymax=Percent_Mutant+se),
                width=.1)+
  theme_classic(base_size = 20)+
  theme(legend.position = "top")+
  labs(x='Days', y='% p53-mutated cells \n in co-culture', color='')+
  lims(y=c(0,45))
dev.off()
