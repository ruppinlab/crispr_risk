# Version for Github
source('Tools/Step0_Globally_used_Functions_and_Datasets.R')

# Cas9 inducing signature difference by KRAS status
############################################################
# Step 0:: Preprocessing data
############################################################
# hallmark genes
hallmark_genes=read.gmt('Data/h.all.v7.1.symbols.gmt')
# previously identified CDE genes
CDE=read.csv('Data/CDE_allThree.csv')
CDE=lapply(CDE, function(x) as.character(x[x!='']))
CDE=CDE[c(1, 2, 5, 6)]
exp=readRDS('Data/exp_cas9_ParentalVSCas9Induced.RDS')
annotation=exp[,1:4]; exp=exp[,-(1:4)]
exp_cas9=exp[,grep('cas9',colnames(exp), ignore.case = T)]
exp_Parental=exp[,-grep('cas9',colnames(exp), ignore.case = T)]
colnames(exp_cas9)=colnames(exp_Parental)
# Kras mutation profile
kras_mut=readRDS('Data/kras_mut.RDS')
Diff_inWT=t(sapply(1:nrow(exp), function(x)
  c(pvalue=wilcox.test(as.numeric(exp_cas9[x,which(kras_mut==0)]),
                       as.numeric(exp_Parental[x,which(kras_mut==0)]))$p.value,
    FC_wtbymut=median(as.numeric(exp_cas9[x,which(kras_mut==0)])/
      as.numeric(exp_Parental[x,which(kras_mut==0)]), na.rm = T))))
Diff_inMut=t(sapply(1:nrow(exp), function(x)
  c(pvalue=wilcox.test(as.numeric(exp_cas9[x,which(kras_mut==1)]),
                       as.numeric(exp_Parental[x,which(kras_mut==1)]))$p.value,
    FC_wtbymut=median(as.numeric(exp_cas9[x,which(kras_mut==1)])/
      as.numeric(exp_Parental[x,which(kras_mut==1)]), na.rm = T))))

############################################################
# hallmark Pathway - Differnetially regulated by KRAS
############################################################
diffscore_inMut=Diff_inMut[,2]
names(diffscore_inMut)=annotation$Gene_Symbol
mut_pathways_enrichment=enriched_pathway(score = diffscore_inMut,
                                         geneList = names(diffscore_inMut),
                                         pathways_list=hallmark_genes)
diffscore_inWT=Diff_inWT[,2]
names(diffscore_inWT)=annotation$Gene_Symbol
WT_pathways_enrichment=enriched_pathway(score = diffscore_inWT,
                                        geneList = names(diffscore_inWT),
                                        pathways_list=hallmark_genes)

df2plot=rbind(data.frame(mut_pathways_enrichment, status='Mut'),
              data.frame(WT_pathways_enrichment, status='WT'))
df2plot$pathway=substring(df2plot$pathway, 10, )
df2plot$POI=''
df2plot$POI[df2plot$pathway %in% df2plot[df2plot$padj<0.05,]$pathway]=
df2plot$pathway[df2plot$pathway %in% df2plot[df2plot$padj<0.05,]$pathway]
df2plot$flipped=FALSE
POI=intersect(df2plot$pathway[df2plot$padj<0.05 & df2plot$NES>0 &
                            df2plot$status=='WT'],
          df2plot$pathway[df2plot$padj<0.05 & df2plot$NES<0 & 
                            df2plot$status=='Mut'])
POI=c(POI, intersect(df2plot$pathway[df2plot$padj<0.05 &
                            df2plot$status=='WT'],
          df2plot$pathway[df2plot$padj>0.05 & 
                            df2plot$status=='Mut']))
POI=c(POI, intersect(df2plot$pathway[df2plot$padj>0.05 &
                                       df2plot$status=='WT'],
                     df2plot$pathway[df2plot$padj<0.05 & 
                                       df2plot$status=='Mut']))
df2plot$flipped[df2plot$pathway%in% POI]=TRUE
levels(df2plot$status)=c('KRAS Mut', 'KRAS WT')
df2plot=df2plot[order(df2plot$pathway),]
tiff('Plots/FigureS9A.tiff',
     width = 1200, height = 800)
ggplot(df2plot, aes(x=NES, y=-log(padj), label=POI, color= !flipped))+
  geom_point()+
  facet_grid(~status)+
  geom_label_repel()+
  geom_hline(yintercept =  - log(0.05), color='red', linetype='longdash')+
  geom_vline(xintercept =  0, color='black', size=2)+
  theme_bw(base_size = 20)+
  theme(legend.position = 'none')+
  labs( y= '-log(FDR)', x='effect size(NES)',
        title = 'Differential Expression Cas9 vs Parental CCL')
dev.off()
############################################################
# Diff in NES
############################################################
df2plot2=data.frame(pathway=unique(df2plot$pathway),
           NES_diff=split(df2plot$NES, df2plot$status)[[1]] -
             split(df2plot$NES, df2plot$status)[[2]])
df2plot2$`Differential regulation direction\n by KRAS Status`=FALSE
df2plot2$`Differential regulation direction\n by KRAS Status`[df2plot2$pathway %in% df2plot$POI[df2plot$flipped]]=TRUE
df2plot2$`Differential regulation direction\n by KRAS Status`=factor(df2plot2$`Differential regulation direction\n by KRAS Status`)
df2plot2$`Differential regulation direction\n by KRAS Status`=factor(df2plot2$`Differential regulation direction\n by KRAS Status`, levels = c('TRUE', 'FALSE'))
tiff('Plots/FigureS9B', height=600, width=800)
ggplot(df2plot2, aes(y=reorder(pathway, -abs(NES_diff)),
                     x=NES_diff,
                     fill=`Differential regulation direction\n by KRAS Status`))+
  geom_bar(stat='identity')+
  theme_bw(base_size = 12)+
  labs(y='Hallmark Pathways', x='NES (WT - mutant)')+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
