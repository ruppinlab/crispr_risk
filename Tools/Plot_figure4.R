# Version for Github
setwd('/Users/sinhas8/crispr_risk-master/')
source('Tools/Step0_Globally_used_Functions_and_Datasets.R')
####################################################################################
# Plot Figure 4A
####################################################################################
# Import output of TIDE assay measuring RNP efficiency after CRISPR-KO in p53 WT vs Mutant
figure4_data=readxl::read_xlsx('Data/RNP_CT_TIDE.xlsx', sheet = 1)
melted <- melt(figure4_data[,c(1,2,3,4,6)], id.vars=c("P53 Status",
                                                        "Timepoint",
                                                        'sgRNA','Cell line'))
melted$`P53 Status`=factor(melted$`P53 Status`)
levels(melted$`P53 Status`)=c('p53-Mutant', 'p53-Wildtype')
df2plot=ddply(melted, c("`P53 Status`","Timepoint",'sgRNA','variable','`Cell line`'),
              summarise,
              mean = mean(value), sd = sd(value),
              sem = sd(value)/sqrt(length(value)))
df2plot$`P53 Status`=factor(df2plot$`P53 Status`)
levels(df2plot$`P53 Status`)=c('p53-Mutant', 'p53-Wildtype')
# manual position to show p-value 
p_df=data.frame(p_value=c(0.2, 0.1, 0.077, 0.059),
                `P53 Status`=rep(c('p53-Mutant', 'p53-Wildtype'), 2),
                `Cell line`=rep(c('MOLM13', 'RPE1'), each=2))
p_df$Timepoint='Day 0'
# x-y position of p-values
p_df$x_value=1.5
p_df$y_value=c(25, 40, 12.5, 6)
colnames(p_df)[2:3]=c('P53 Status', 'Cell line')
p_df$FC=aggregate(mean ~ `P53 Status`+`Cell line`,df2plot, function(x) x[1]-x[2])[,3]

tiff('Plots/figure4A.tiff', height=400, width=400)
ggplot(df2plot, aes(x=Timepoint,fill=Timepoint,y=mean))+
  facet_wrap(`Cell line`~`P53 Status`, scales = 'free')+
  geom_bar(stat='identity')+
  geom_point(data = melted, aes(y=value, x=Timepoint))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  theme_bw(base_size = 13)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_text(data=p_df, aes(x=x_value,y=y_value,label=paste('p=',p_value,
                                                           ';',
                                                           'fc=',round(FC,2))))
dev.off()
########################################################################
# Figure 4B
########################################################################
# Import output of TIDE assay measuring RNP efficiency after CRISPR-KO in p53 WT vs Mutant
figure4_dataB=readxl::read_xlsx('Data/RNP_CT_TIDE.xlsx', sheet = 2)
figure4_dataB_dt=data.table(figure4_dataB)
figure4_dataB_dt=figure4_dataB_dt[!(figure4_dataB_dt$Timepoint=='D0' |
                                      figure4_dataB_dt$Timepoint=='D2'|
                                      figure4_dataB_dt$Timepoint=='D4'),]
# Dataframe to add p-values
p_dataframe=rbind(data.frame(pvalue=sapply(unique(figure4_dataB_dt$Timepoint),
       function(x) t.test(figure4_dataB_dt$`APC Mean Fluorescence Intensity (MFI)`[figure4_dataB_dt$`P53 Status`=='Wildtype'
                                                                                        & figure4_dataB_dt$sgRNA=='NDUFB6'
                                                                                        & figure4_dataB_dt$Timepoint==x],
                               figure4_dataB_dt$`APC Mean Fluorescence Intensity (MFI)`[figure4_dataB_dt$`P53 Status`=='Wildtype'
                                                                                        & figure4_dataB_dt$sgRNA=='NTC'
                                                                                        & figure4_dataB_dt$Timepoint==x])$p.value),
       `P53 Status`='Wildtype'),
      data.frame(pvalue=sapply(unique(figure4_dataB_dt$Timepoint),
                               function(x) t.test(figure4_dataB_dt$`APC Mean Fluorescence Intensity (MFI)`[figure4_dataB_dt$`P53 Status`=='Mutant'
                                                                                                           & figure4_dataB_dt$sgRNA=='NDUFB6'
                                                                                                           & figure4_dataB_dt$Timepoint==x],
                                                  figure4_dataB_dt$`APC Mean Fluorescence Intensity (MFI)`[figure4_dataB_dt$`P53 Status`=='Mutant'
                                                                                                           & figure4_dataB_dt$sgRNA=='NTC'
                                                                                                           & figure4_dataB_dt$Timepoint==x])$p.value),
                 `P53 Status`='Mutant'))
p_dataframe$Timepoint=rep(unique(figure4_dataB_dt$Timepoint), 2)
melted <- melt(figure4_dataB_dt[,c(2,3,4,6)], id.vars=c("P53 Status",
                                        "Timepoint",
                                        'sgRNA'))
df2plot=ddply(melted, c("`P53 Status`","Timepoint",'sgRNA','variable'),
      summarise,
      mean = mean(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)))
df2plot$day=as.numeric(gsub('D','',df2plot$Timepoint))
df2plot$Timepoint=factor(df2plot$Timepoint)
df2plot$Timepoint=factor(df2plot$Timepoint, 
                         levels=levels(df2plot$Timepoint)[
                           order(as.numeric(gsub('D','',levels(df2plot$Timepoint))))])
colnames(p_dataframe)[2]=c('P53 Status')
p_dataframe$sgRNA='NTC'
FC=aggregate(mean ~ Timepoint+`P53 Status`,df2plot, function(x) x[1]/x[2])
p_dataframe=p_dataframe[order(paste(p_dataframe$Timepoint, p_dataframe$`P53 Status`),
                  decreasing =T),]
p_dataframe$FC=FC[order(paste(FC$Timepoint, FC$`P53 Status`),
         decreasing =T),3]
p_dataframe$`P53 Status`=factor(p_dataframe$`P53 Status`)
levels(p_dataframe$`P53 Status`)=c('p53-Mutant', 'p53-Wildtype')
df2plot$`P53 Status`=factor(df2plot$`P53 Status`)
levels(df2plot$`P53 Status`)=c('p53-Mutant', 'p53-Wildtype')
y_value=aggregate(mean ~ Timepoint+`P53 Status`,df2plot, max)
p_dataframe$y_value=y_value[order(paste(y_value$Timepoint, FC$`P53 Status`),
                        decreasing =T),3]+1000
p_dataframe$Timepoint=factor(p_dataframe$Timepoint)
p_dataframe$Timepoint=factor(p_dataframe$Timepoint,
                             levels=levels(p_dataframe$Timepoint)[
                               order(as.numeric(gsub('D','',levels(p_dataframe$Timepoint))))])

figure4_dataB_dt$`P53 Status`=factor(figure4_dataB_dt$`P53 Status`)
levels(figure4_dataB_dt$`P53 Status`)=c('p53-Mutant', 'p53-Wildtype')


tiff('Plots/figure4B.tiff', height=400, width=500)
ggplot(df2plot, aes(x=sgRNA,y=mean, fill=sgRNA)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_point(data = figure4_dataB_dt,
             aes(y=`APC Mean Fluorescence Intensity (MFI)`, fill=sgRNA))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  facet_grid(`P53 Status` ~ Timepoint,  scales = 'free', shrink = T)+
  theme_bw(base_size = 12)+
  theme(legend.position="top")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_brewer(palette = 3)+
  geom_text(data=p_dataframe, aes(x=1.5,y=y_value,
                                  label=paste('\np=',format.pval(pvalue, 2),
                                              ';\n','fc=',round(FC, 2),
                                              sep='')))+
  labs(y='APC Mean Fluorescence Intensity (MFI)', x='Timepoint')
dev.off()
