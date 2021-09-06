# Version for Github
source('Step0_Globally_used_Functions_and_Datasets.R')
# Plotting figure S11
####################################################################################
##Functions required to plot
####################################################################################
Panel3<-function(base_dataset=achilles,
                 which_hit=1,
                 Dataset='shRNA \n Viability',
                 is.Mut=FALSE){
  #	Sig_Score=significance_scores[[1]][[1]]
  #	names(Sig_Score)=rownames(avana)
  #	WTgreaterMT_genes=rownames(avana)[which(significance_scores[[1]][[1]]<0.1)]
  #	WTgreaterMT_genes_control=rownames(avana)[which(significance_scores[[1]][[3]]<0.1)]
  #	CRISPR_specific=which(significance_scores[[1]][[1]]<0.1)[which(significance_scores[[1]][[1]]<0.1) %!in% which(significance_scores[[1]][[3]]<0.1)]
  #	CRISPR_Hits_in_Sig.order=(Sig_Score[CRISPR_specific])[order(Sig_Score[CRISPR_specific])]
  
  Hit=CDE_Pos_Neg[[1]][which_hit]
  baseorder=order(avana[Hit,])
  forboxes=xtile(unlist(avana[Hit, baseorder]), 8)
  df=data.frame(base=unlist(base_dataset[Hit, baseorder]), forboxes)
  if(!is.Mut & Dataset!='CNV')
    ggplot(df, aes(x=as.factor(forboxes), y=base, group=forboxes))+
    geom_boxplot()+
    labs(x='', y=paste(Dataset))+
    theme(axis.text=element_text(size=14))
  else if(Dataset=='CNV'){
    ggplot(df, aes(x=seq(nrow(df)), y=base))+
      geom_smooth()+
      geom_line()+
      labs(x='CellLines', y='Copy Number')+
      ylim(0, 4)+
      theme(axis.text=element_text(size=14))
  }
  else{
    TP53=base_dataset['TP53', baseorder]
    TP53=factor(TP53, labels=c('WT', 'MT'))
    df=data.frame(TP53, forboxes)
    ggplot(df)+
      geom_bar(aes(x=as.factor(forboxes), fill=TP53))+
      labs(x='', y='TP53 Mutation')+
      theme(axis.text=element_text(size=14),
            legend.text=element_text(size=14),
            legend.title=element_text(size=14),
            legend.position="top",
            legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(-10,-10,-10,-10))
  }
}

####################################################################################
##calling funcs for the respective master regulator of interest
####################################################################################
Master_Regulator_ofInterest='TP53'
CDE_Pos_Neg=readRDS('../Data/CDE_for_volg.RDS')[[Master_Regulator_ofInterest]]
Mut_Prep=t(apply(Mut_CCLE,1, function(x) as.numeric(as.logical(x))))

## Various Panels of the figure are created sep to combine
which_hit=1
g1 <- ggplotGrob(Panel3(achilles, which_hit,'shRNA \n Viability'))
g2 <- ggplotGrob(Panel3(avana, which_hit,'CRSIPR Viability'))
g3 <- ggplotGrob(Panel3(Mut_Prep, which_hit,'TP53 Mutation Status',is.Mut=TRUE))
g4 <- ggplotGrob(Panel3((2*(2^(CNV))), which_hit,'CNV'))

## Combining and aligning the above part
g <- rbind(g1, g2, g3, g4, size = "last")
g$widths <- unit.pmax(g1$widths, g4$widths, g2$widths, g3$widths)

#Save in a File
tiff('../Plots/figureS11.tiff')
grid.newpage()
grid.draw(g)
dev.off()

