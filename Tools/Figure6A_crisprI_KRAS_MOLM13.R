# CRISPRi KRAS MOLM
df2plot=rbind(data.frame(geneName=KRAS_WT2$Gene[GOI],
                         RANK=KRAS_WT2$FC_rank[GOI], 
                         type='WT',
                         group=KRAS_Mut$Gene.group[GOI],
                         sgRNA=KRAS_WT2$sgRNA[GOI]),
              data.frame(geneName=KRAS_WT2$Gene[GOI],
                         RANK=KRAS_Mut$FC_rank[GOI],
                         type='Mut',
                         group=KRAS_Mut$Gene.group[GOI],
                         sgRNA=KRAS_WT2$sgRNA[GOI]))
iKRAS_WT=readxl::read_xlsx('/Users/sinhas8/Project_CRISPR/Summary_iKRASCDERuns (2).xlsx')
iKRAS_WT$FC_rank=rank(iKRAS_WT$`Fold change D30/D0`)
iKRAS_WT_matched1=iKRAS_WT[!is.na(match(iKRAS_WT$Gene, df2plot$geneName)),]
iKRAS_Mut=readxl::read_xlsx('/Users/sinhas8/Project_CRISPR/Summary_iKRASCDERuns (2).xlsx', sheet = 2)
iKRAS_Mut$FC_rank=rank(iKRAS_Mut$`Fold change D30/D0`)
iKRAS_Mut_matched1=iKRAS_Mut[!is.na(match(iKRAS_Mut$Gene, df2plot$geneName)),]

df2plot=rbind(data.frame(system='CRISPRi',
                         type='Mut',
                         FC=iKRAS_Mut_matched1$`Fold change D30/D0`,
                         FC_rank=iKRAS_Mut_matched1$FC_rank,
                         Group=iKRAS_Mut_matched1$`Gene group`),
              data.frame(system='CRISPRi',
                         type='WT',
                         FC=iKRAS_WT_matched1$`Fold change D30/D0`,
                         FC_rank=iKRAS_WT_matched1$FC_rank,
                         Group=iKRAS_WT_matched1$`Gene group`))

saveRDS(df2plot, '/Users/sinhas8/Project_CRISPR/df2plot_iKRAS_matched_cas9_genes.RDS')
df2plot=readRDS('/Users/sinhas8/Project_CRISPR/df2plot_iKRAS_matched_cas9_genes.RDS')
require(ggpubr)
df2plot
tiff('/Users/sinhas8/Project_CRISPR/Figure4G_Panel2.tiff', width = 400, height = 400)
ggplot(df2plot, aes(x=Group, fill=type, y=FC_rank))+geom_boxplot()+
  stat_compare_means(method='wilcox', paired = T)+
  facet_grid(~system)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
    # annotate(geom="text", y=0, x=0.5, label=round(p_value[2], 2), color="blue")+
    # annotate(geom="text", y=0, x=1.5, label=round(p_value[1], 2), color="blue")
