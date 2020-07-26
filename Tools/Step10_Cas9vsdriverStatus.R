# Version for Github
source('Tools/Step0_Globally_used_Functions_and_Datasets.R')

####################################################################################
##Identifying the cas9-activity in WT vs mutant for cancer drivers
####################################################################################
# GFP reporter assay
setwd('/Users/sinhas8/crispr_risk-master/')
cas9_activity=readRDS('Data/cas9_activity.RDS')
summary(lm(as.numeric(as.character(cas9_activity$cas9_activity)) ~ 
             cas9_activity[,'KRAS']+cas9_activity[,'TP53']))

Volg_DifferentialCas9=fdrcorr(apply(gfp_wdMut[,-(1:23)], 2, 
                              function(x) c(p_value=err_handle(
                                wilcox.test(as.numeric(gfp_wdMut[,'cas9_activity']) ~ x,
                                            alternative='l')$p.value))))
