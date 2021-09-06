# Scatter plot for showing the rank of driver genes.
CDE=readRDS('/Users/sinhas8/Project_CRISPR/Project_CRISPR.Risk_Github/Data/CDE_for_volg.RDS')

CDEbased_Test=data.frame(p_value=10^(-mr.by.pval),
                         CDEsize=sapply(CDE, function(x) length(x[[1]]))[names(mr.by.pval)]
                         )
colnames(Volg_Mutselection)=c('Sig', 'Effect Size')
colnames(Volg_DifferentialCas9)=c('Sig', 'Effect Size')
fdrcorr(Volg_DifferentialCas9[,1])['TP53']
colnames(CDEbased_Test)=c('Sig', 'Effect Size')
Volg_DifferentialCas9=Volg_DifferentialCas9[rownames(Volg_Mutselection),]
CDEbased_Test=CDEbased_Test[rownames(Volg_Mutselection),]

Volg_Mutselection[order(Volg_Mutselection[,1]),]
df2plot=rbind(data.frame(Volg_DifferentialCas9, Type='DiffCas9'),
              data.frame(Volg_Mutselection, Type='MutSelection'),
              data.frame(CDEbased_Test, Type='CDE'))

df2plot=na.omit(df2plot)
levels(df2plot$Type)=c('Differnetially lower cas9',
                       'Mutation Selection after cas9',
                       'Skewness of CDE genes')
saveRDS(df2plot, '/Users/sinhas8/Project_CRISPR/df2plot_figure4_scatterplot.RDS')

tiff('/Users/sinhas8/Project_CRISPR/Threemethods_forMR_figure4A.tiff',
     width = 900, height=300)
ggplot(df2plot, aes(y=Effect.Size,
                    x= -log(Sig, 10),
                    color=Sig<0.05))+
  geom_point()+
  facet_wrap(~Type, scales = 'free')+
  theme_bw(base_size = 15)+
  geom_text_repel(aes(label=ifelse((Sig<0.05 & Type=='Differnetially lower cas9') |
                                     (fdrcorr(Sig)<0.001 & Effect.Size>100 & Type=='Skewness of CDE genes')|
                                     (Sig<0.1 & Type=='Mutation Selection after cas9' | 
                                        rownames(Volg_Mutselection) =='TP53'), rownames(Volg_Mutselection),'')))+
  scale_color_manual(values = c('Black', 'Red'))+
  geom_vline(xintercept = 1.3, linetype='dashed', color='red')+
  theme(legend.position = 'none')
dev.off()  
  


fdrcorr(sapply(rownames(Volg_DifferentialCas9), function(x) 
  Mut_Assoc_given_a_mutmatrix('TP53', x, mut_matrix = cas9Activity))
)