# Load GSP and GSN standards
setwd('/Users/sinhas8/Downloads/msb188679-sup-0001-Data_and_Code')
nonEss <- read.csv("training_nonessential.txt",header=T,sep="\t",stringsAsFactors=F);
ess <- read.csv("core-v2-gene-hgnc",header=F,sep="\t",stringsAsFactors=F);

###
# Figure 1E - PR Curves
###
tkov3 <- read.csv("wt_rpe1_screens_tkov3.txt",header=T,sep="\t",stringsAsFactors=F);
# Load data file from Hart et al Cell 2015 paper
mat2 <- read.csv("Hart_et_al_Cell2015_guideLevel_logFC_matrix.txt",header=T,sep="\t",stringsAsFactors=F);

# Compute gene-level logFC matrix
mat2.avg <- aggregate( mat2[,-1:-2],by=list(mat2$GENE),mean);
mat2.avg <- mat2.avg[,c("Group.1","DLD1","HeLa","RPE1","HCT116","GBM")];

# Compute PR curves for datasets
GOI=rownames(avana)
fg1_zimmerman <- -unlist(tkov3$Zimmerman.logFC[ which(tkov3$Gene.symbol %in% GOI) ]);
names(fg1_zimmerman)=tkov3$Gene.symbol[tkov3$Gene.symbol %in% GOI]
fg2_unp <- -unlist(tkov3$RPE1_002.logFC[ tkov3$Gene.symbol %in% GOI ]);
names(fg2_unp)=tkov3$Gene.symbol[tkov3$Gene.symbol %in% GOI]
fg3_unp <- -unlist(tkov3$RPE1_006.logFC[ tkov3$Gene.symbol %in% GOI ]);
names(fg3_unp)=tkov3$Gene.symbol[tkov3$Gene.symbol %in% GOI]

fg4_hart <- -unlist(mat2.avg$RPE1[ mat2.avg$Group.1 %in% GOI ]);
names(fg4_hart)=mat2.avg$Group.1[mat2.avg$Group.1 %in% GOI]

# Get Hart et al p53-null RPE1 screen
# RPE1
p53null <- read.csv("matrix-rpe1_p53Null.txt",header=T,sep="\t",stringsAsFactors=F);
p53null <- p53null[ !apply( p53null[,-1:-2] == 0,1,all ), ];
p53null.norm <- sweep( p53null[,-1:-2],2,apply( p53null[,-1:-2],2,sum)/1000000,FUN="/")
p53null.norm <- as.data.frame( cbind( p53null[,1:2], p53null.norm ));
p53null.norm$RPE1.logFC <- log2( (rowMeans( p53null.norm[,c("Y_RPE1_WT_T18A_S3_count","Y_RPE1_WT_T18B_S4_count")] )+0.1) / (p53null.norm$Y_RPE1_WT_T0_S1_count +0.1) );
p53null.avg <- aggregate( p53null.norm[,c("RPE1.logFC")],by=list( p53null.norm$geneA ),mean );
colnames(p53null.avg) <- c("Gene.symbol","logFC");

fg5_hart_p53Null <- -unlist(p53null.avg$logFC[ p53null.avg$Gene.symbol %in% GOI ]);
names(fg5_hart_p53Null)=p53null.avg$Gene.symbol[p53null.avg$Gene.symbol %in% GOI]

# Finally, redo the Haapaniemi data
data <- read.csv("haapaniemi_NatMed_June2018_rawReadCounts.txt",header=T,sep="\t",stringsAsFactors=F);
data.norm <- data;
data.norm[,-1:-2] <- sweep(data[,-1:-2],2,apply( data[,-1:-2],2,sum)/1000000,FUN="/");
data.norm <- data.norm[ !apply(data.norm[,-1:-2] == 0,1,all), ];
# Compute mean log2FC and plot out GSP and GSN sets
data.norm$RPEnull.logFC <- log2( rowMeans( data.norm[,c("RPEnull_R1_D28","RPE1null_R2_D28")]+0.1 )) - log2( rowMeans( data.norm[,c("RPEnull_R1_D4","RPEnull_R2_D4")]+0.1));
data.norm$RPEwt.logFC <- log2( rowMeans( data.norm[,c("RPEWT_R1_D28","RPEWT_R2_D28")]+0.1 )) - log2( rowMeans( data.norm[,c("RPEWT_R1_D4","RPEWT_R2_D4")]+0.1));
data.norm.geneLevel <- aggregate( data.norm[,c("RPEnull.logFC","RPEwt.logFC")],by=list( data.norm$Gene ),mean );

fg_haap_p53null <- -unlist(data.norm.geneLevel$RPEnull.logFC[ data.norm.geneLevel$Group.1 %in% GOI ])
names(fg_haap_p53null)=data.norm.geneLevel$Group.1[data.norm.geneLevel$Group.1 %in% GOI]
fg_haap_p53WT <- -unlist(data.norm.geneLevel$RPEwt.logFC[ data.norm.geneLevel$Group.1 %in% GOI ]);
names(fg_haap_p53WT)=data.norm.geneLevel$Group.1[data.norm.geneLevel$Group.1 %in% GOI]

all_screens=list(fg1_zimmerman, fg2_unp, fg3_unp, 
                        fg4_hart, fg5_hart_p53Null,
                        fg_haap_p53null,fg_haap_p53WT)

saveRDS(all_screens, '/Users/sinhas8/Project_CRISPR/2.Data/seven_RPE_crisprcas9_screen.RDS')
names(all_screens)=c('p53 WT fg1_zimmerman',
                     'p53 WT fg2_unp',
                     'p53 WT fg3_unp', 
                     'p53 WT fg4_hart',
                     'p53 Null fg5_hart',
                     'p53 Null fg6_haapameini',
                     'p53 WT fg7_haapameini')
common_genes=Reduce(intersect, lapply(list(fg1_zimmerman, fg2_unp, fg3_unp, fg4_hart, fg5_hart_p53Null,
     fg_haap_p53null,fg_haap_p53WT), names))
CDE=read.csv('/Users/sinhas8/Project_CRISPR/2.Data/cde')
sapply(1:7, function(x) wilcox.test(vectorSubset(all_screens[[6]], CDE$p53.CDE.),
                                    vectorSubset(all_screens[[x]], CDE$p53.CDE.),
                                    alternative='l')$p.value)
# length(na.omit(unique(hart_genes$`ConstitutiveCoreEssentials(CCE)`)))
# length(hart_genes$`Nonessential Genes (NE)`)
# 217/(length(na.omit(unique(hart_genes$`ConstitutiveCoreEssentials(CCE)`)))+
#        length(hart_genes$`Nonessential Genes (NE)`))


all_screens_CDE_scores=sapply(all_screens, function(x) x[match(CDE$p53.CDE.[CDE$p53.CDE.!=''], names(x))] )

colnames(all_screens_CDE_scores)=c('p53 WT fg1_zimmerman',
                                   'p53 WT fg2_unp',
                                   'p53 WT fg3_unp', 
                                   'p53 WT fg4_hart',
                                   'p53 Null fg5_hart',
                                   'p53 Null fg6_haapameini',
                                   'p53 WT fg7_haapameini')
require(reshape2)
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

# Currently Figure S5
tiff('/Users/sinhas8/Project_CRISPR/figureS7.tiff')
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