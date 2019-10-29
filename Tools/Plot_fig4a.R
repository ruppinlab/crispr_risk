# Version for Github
source('Step0_Globally_used_Functions_and_Datasets.R')

####################################################################################
# Plot 4A
# Using CDE genes positive and negative to identify the Master Regualtors
####################################################################################
##Variable comprising saved CDE genes for each genes
CDE=readRDS('../Data/CDE_for_volg.RDS')
CDE_pos_size=sapply(CDE, function(x) length(x[[1]]))
####################################################################################
## Calcualte the significnace of imbalanmce of DE+/- genes
####################################################################################
Vog_Gene_df=get(load('../Data/vogelstein.RData'))
Vog_Genes=Vog_Gene_df$symbol
mat=readRDS('../Data/DE_posNneg_genes_for_VolgGenes.RDS')
Contigency=t(sapply(mat, function(x) sapply(x, nrow)))
rownames(Contigency)=Vog_Genes
# Significance of imblance of DE+/- in CRISPR and shRNA
Sig=p.adjust(unlist(apply(Contigency, 1, function(x)  significance_test(x))), method='fdr')
# df of log(Sig) and potential Master regulator tested
Log_Sig=-log(Sig, 10)
Log_Sig_df=data.frame(Vog_Genes, Log_Sig)
Log_Sig_df=Log_Sig_df[order(Log_Sig_df$Log_Sig),]
Log_Sig_df$Log_Sig[Log_Sig_df$Log_Sig>300]=300

Log_Sig_df$CDE_pos=CDE_pos_size[match(Log_Sig_df$Vog_Genes, names(CDE_pos_size))]
Log_Sig_df$col=''
Log_Sig_df$col[Log_Sig_df$CDE_pos>100]='Hits of Interest'

##plotting CDEsize vs significance
plot1=ggplot(Log_Sig_df, aes(x=Log_Sig, y=CDE_pos, color=col))+
  geom_point(aes(size=ifelse(CDE_pos>1 | Log_Sig>50 , 3, 2)))+
  labs(x=expression(-log[10](p)), y='Number of Genes\nat Significant Risk')+
  geom_text_repel(
    data = subset(Log_Sig_df, CDE_pos > 1 | Log_Sig>50 ), nudge_y=100,
    aes(label = Vog_Genes),
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )+
  scale_color_manual(values=c("black", "#56B4E9"), guide=FALSE)+
  theme_classic()+
  scale_size(guide = 'none')
plot2 <- ggplot(Log_Sig_df, aes(Log_Sig)) + 
  geom_density(alpha=.5)+theme_classic()+
  labs(x=expression(-log[10](p)), y='Density')

joined_plot=grid.arrange(plot2, plot1, 
                         ncol=1, nrow=2, heights=c(1,3))
tiff('../Plots/figure4A.tiff', 480, 220)
plot(joined_plot)
dev.off()

