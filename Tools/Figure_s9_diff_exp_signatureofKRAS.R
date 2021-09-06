# Test diff exp genes
# # Identify Pathways

saveRDS(kras_mut, 'Data/kras_mut.RDS')
# Define function to idnetify a differnetial signature
differential_signature_wtvsMut<-function(geneName='KRAS'){
  table(kras_mut)
  names(kras_mut)=colnames(exp_cas9)
  Diff_inWT=t(sapply(1:nrow(exp), function(x)
    c(pvalue=wilcox.test(as.numeric(exp_cas9[x,which(kras_mut==0)]),
                         as.numeric(exp_Parental[x,which(kras_mut==0)]))$p.value,
      FC_wtbymut=median(as.numeric(exp_cas9[x,which(kras_mut==0)]), na.rm = T)/
        median(as.numeric(exp_Parental[x,which(kras_mut==0)]), na.rm = T))))
  # Diff=t(sapply(1:nrow(exp), function(x)
  #   c(pvalue=wilcox.test(as.numeric(exp_cas9[x,]),
  #                        as.numeric(exp_Parental[x,]))$p.value,
  #     FC_wtbymut=median(as.numeric(exp_cas9[x,]), na.rm = T)/
  #       median(as.numeric(exp_Parental[x,]), na.rm = T))))
  
  Diff_inMut=t(sapply(1:nrow(exp), function(x)
    c(pvalue=wilcox.test(as.numeric(exp_cas9[x,which(kras_mut==1)]),
                         as.numeric(exp_Parental[x,which(kras_mut==1)]))$p.value,
      FC_wtbymut=median(as.numeric(exp_cas9[x,which(kras_mut==1)]), na.rm = T)/
        median(as.numeric(exp_Parental[x,which(kras_mut==1)]), na.rm = T))))
  # Mut vs WT
  diffscore=Diff_inMut[,2]
  names(diffscore)=annotation$Gene_Symbol
  mut_pathways_enrichment=enriched_pathway()
  diffscore=Diff_inWT[,2]
  names(diffscore)=annotation$Gene_Symbol
  WT_pathways_enrichment=enriched_pathway()
  WT_pathways_enrichment=WT_pathways_enrichment[match(mut_pathways_enrichment$pathway,WT_pathways_enrichment$pathway),]
  df2plot=data.frame(pathway=WT_pathways_enrichment$pathway,
                     wt_nes=WT_pathways_enrichment$NES, 
                     mut_nes=mut_pathways_enrichment$NES)

  # Plot the above
  tiff(paste('/Users/sinhas8/Project_CRISPR/pathway_difference_',geneName,'.tiff'), width = 800, height=800)
  ggplot(df2plot, aes(x=wt_nes, y=mut_nes, label=pathway))+
    geom_text_repel(size=3.5)+
    geom_point()+
    theme_bw(base_size = 15)+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)
  dev.off()
}
differential_signature_wtvsMut<-function(geneName='TP53'){
  kras_mut=matched_mutmatrix[geneName, match(colnames(exp_cas9),
                                           CCLE_cellLineName)]
  
  names(kras_mut)=colnames(exp_cas9)
  Diff_inWT=t(sapply(1:nrow(exp), function(x)
    c(pvalue=wilcox.test(as.numeric(exp_cas9[x,which(kras_mut==0)]),
                         as.numeric(exp_cas9[x,which(kras_mut==1)]))$p.value,
      FC_wtbymut=median(as.numeric(exp_cas9[x,which(kras_mut==0)]), na.rm = T)/
        median(as.numeric(exp_cas9[x,which(kras_mut==1)]), na.rm = T))))
  Diff_inMut=t(sapply(1:nrow(exp), function(x)
    c(pvalue=wilcox.test(as.numeric(exp_Parental[x,which(kras_mut==0)]),
                         as.numeric(exp_Parental[x,which(kras_mut==1)]))$p.value,
      FC_wtbymut=median(as.numeric(exp_Parental[x,which(kras_mut==0)]), na.rm = T)/
        median(as.numeric(exp_Parental[x,which(kras_mut==1)]), na.rm = T))))
  
  # Mut vs WT
  diffscore=Diff_inMut[,2]
  names(diffscore)=annotation$Gene_Symbol
  mut_pathways_enrichment=enriched_pathway(pathways_list = DDR_Geneset_List)
  diffscore=Diff_inWT[,2]
  names(diffscore)=annotation$Gene_Symbol
  WT_pathways_enrichment=enriched_pathway(pathways_list = DDR_Geneset_List)
  WT_pathways_enrichment=WT_pathways_enrichment[match(mut_pathways_enrichment$pathway,WT_pathways_enrichment$pathway),]
  head(WT_pathways_enrichment)
  df2plot=data.frame(pathway=WT_pathways_enrichment$pathway,
                     wt_nes=WT_pathways_enrichment$NES, 
                     mut_nes=mut_pathways_enrichment$NES)
  head(df2plot[order(df2plot$wt_nes - df2plot$mut_nes, decreasing = T),], 20)
  # Plot the above
  tiff(paste('/Users/sinhas8/Project_CRISPR/pathway_difference_cas9vsParental_DDR_',geneName,'.tiff'))
  ggplot(df2plot, aes(x=wt_nes, y=mut_nes, label=pathway))+
    geom_text_repel(size=3.5)+
    geom_point()+
    theme_bw(base_size = 15)+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    labs(x='NES cas9 (p53 WT vs Mut)', y='NES parental (p53 WT vs Mut)')
  dev.off()
  
}

# Test for each pathway - How many samples is it upregulated in
