# Version for Github
source('Step0_Globally_used_Functions_and_Datasets.R')

#Plot Figure 3B
########################################################################
#######Load and Preprocess the data
########################################################################
compAssay_Page2=readxl::read_xlsx('../Data/Competition_Analysis3.xlsx', sheet = 2)
long_df2plot=gather(compAssay_Page2, Day, Percent_Mutant, D0:D25, factor_key=TRUE)
long_df2plot=long_df2plot[long_df2plot$Gene %in% levels(factor(long_df2plot$Gene))[3:5],]
##No Day 31 score available
long_df2plot=long_df2plot[,-7]
long_df2plot=as.data.frame(long_df2plot)
#Data.Table to Data frame
long_df2plot=do.call(rbind, lapply(split(long_df2plot,
                                         list(long_df2plot$System,
                                              long_df2plot$Gene, 
                                              long_df2plot$Day,
                                              long_df2plot$Mix)),
                                   function(x) 
                                     err_handle(data.frame(x[,c(1:7)],
                                                           Percent_Mutant=x$Percent_Mutant ))
))

summarydf2plot = summarySE(long_df2plot, measurevar="Percent_Mutant", groupvars=c("System","Mix", "Day", 'Gene'))
summarydf2plot$Gene[summarydf2plot$Gene=='NTC']='Non-Targeting Control'
summarydf2plot$Day_num= as.numeric(sapply(summarydf2plot$Day, function(x) substring(as.character(x), 2)))
summarydf2plot$Day_num=factor(summarydf2plot$Day_num)
summarydf2plot$System=factor(summarydf2plot$System, labels = c('CRISPRi', 'CRSIPR-KO'))
long_df2plot$Gene[long_df2plot$Gene=='NTC']='Non-Targeting Control'
long_df2plot$Day_num= as.numeric(sapply(long_df2plot$Day, function(x) substring(as.character(x), 2)))
long_df2plot$Day_num=factor(long_df2plot$Day_num)
long_df2plot$System=factor(long_df2plot$System, labels = c('CRISPRi', 'CRSIPR-KO'))
########################################################################
##Only Plotting for 5/95 mixture
########################################################################
######Plot
tiff('../Plots/figure3B.tiff', width=300, height=600)
ggplot(summarydf2plot[summarydf2plot$Mix=='5/95',], aes(y=Percent_Mutant, x=Day, color=System))+
  geom_point()+
  geom_line(data = summarydf2plot[summarydf2plot$Mix=='5/95',],
            aes(x=as.numeric(Day), y=Percent_Mutant), size=1.3)+
  facet_wrap(Gene~., , nrow=3, scales='free')+
  geom_errorbar(aes(ymin=Percent_Mutant-se, ymax=Percent_Mutant+se),
                 width=.1)+
   theme_classic(base_size = 12)+
   theme(, legend.position = "top")+
  stat_compare_means(data=long_df2plot, aes(group = System), 
                     label='p.signif', size=4, label.y.npc = 0.85)+
  labs(x='Days', y='% p53-mutated cells in co-culture', color='')
dev.off()
########################################################################
#######Supp Figure S2
########################################################################
compAssay_Page2=readxl::read_xlsx('../Data/Competition_Analysis3.xlsx', sheet = 2)
long_df2plot=gather(compAssay_Page2, Day, Percent_Mutant, D0:D25, factor_key=TRUE)
long_df2plot=as.data.frame(long_df2plot)
summarydf2plot = summarySE(long_df2plot, measurevar="Percent_Mutant", groupvars=c("System","Mix", "Day", 'Gene'))
summarydf2plot=summarydf2plot[summarydf2plot$Gene=='ANKRD49' | summarydf2plot$Gene=='FANCG'| summarydf2plot$Gene=='TFAP4',]
summarydf2plot$Gene=factor(as.character(summarydf2plot$Gene))

tiff('../Plots/figureS2.tiff', width=800, height=300)
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

