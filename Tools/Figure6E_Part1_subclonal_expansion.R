# Expansion of sub-clonal mutation
# onTarget=readRDS('/Users/sinhas8/Project_OffTarget/2.Data/onTarget_v2.RDS')
source('/Users/sinhas8/myCustom_functions.R')
library(statar)

comp_natGen=readxl::read_xlsx('/Users/sinhas8/Project_CRISPR/comp_natGen_Enache.xlsx')
# colnames(comp_natGen)
remove_silent=comp_natGen$Canonical_Variant_Classification!='synonymous'
comp_natGen=comp_natGen[remove_silent,]
comp_natGen$variant_ID= paste(comp_natGen$Cell_Line_Origin,
                  comp_natGen$Hugo_Symbol,
                  comp_natGen$Position, 
                  comp_natGen$Canonical_Variant_Classification, 
                  comp_natGen$Category,
                  comp_natGen$Canonical_cDNA_Change,
                  comp_natGen$Canonical_Protein_Change,
                  sep='_')
comp_natGen$variant_ID[grep('KRAS',comp_natGen$variant_ID)]
sum(sort(table(comp_natGen$Canonical_Variant_Classification)))
# both_in_cas9vsparent=names(which(table(comp_natGen$variant_ID)==2))
# Identify the cas9 high cell lines
cas9Activity=readxl::read_xlsx('/Users/sinhas8/Project_CRISPR/cas9_Activity_inCellLines.xlsx')
# lowCas9_CellLines=cas9Activity$`Official Cell Line Name`[which(as.numeric(cas9Activity$`%GFP positive-cells`)<20)]
lowCas9_CellLines=cas9Activity$`Official Cell Line Name`[xtile(as.numeric(cas9Activity$`%GFP positive-cells`), 4)==4]
CDE=readRDS('/Users/sinhas8/Project_CRISPR/Project_CRISPR.Risk_Github/Data/CDE_for_volg.RDS')
test_enrichment<-function(infunc_geneName='KRAS'){
  kras_mut=comp_natGen[comp_natGen$Hugo_Symbol==infunc_geneName 
                       &
                         comp_natGen$Canonical_Variant_Classification!='intron'
                       ,]
  kras_mut$Cell_Line_Origin
  kras_mut_subset=do.call(rbind, lapply(split(kras_mut, as.character(kras_mut$variant_ID)),
                                        function(x) c(x$Allele_Fraction[match('CAS9',x$`Parental/Cas9`)],
                                                      x$Allele_Fraction[match('PARENTAL', x$`Parental/Cas9`)])))
  colnames(kras_mut_subset)=c('CAS9', 'PARENTAL')
  kras_mut_subset[is.na(kras_mut_subset)]=0
  not_all_mut=which(kras_mut_subset[,2]!=1)
  cellLinesName=sapply(rownames(kras_mut_subset), function(x) strsplit(x, '_')[[1]][1])
  highCas9_CellLines=which(is.na(match(cellLinesName, lowCas9_CellLines)))
  
  c(wilcox.test(kras_mut_subset[intersect(highCas9_CellLines, not_all_mut),1],
              kras_mut_subset[intersect(highCas9_CellLines, not_all_mut),2],
              paired = T,
              alternative = 'g')$p.value,
    median(kras_mut_subset[intersect(highCas9_CellLines, not_all_mut),1] -
             kras_mut_subset[intersect(highCas9_CellLines, not_all_mut),2], na.rm=T)
    )
}

Volg_Mutselection=do.call(rbind, sapply(names(CDE), function(x) err_handle(test_enrichment(x)) ))
Volg_Mutselection=Volg_Mutselection[order(Volg_Mutselection[,1]),]
Volg_Mutselection[which(rownames(Volg_Mutselection[order(Volg_Mutselection[,1]),])=='KRAS'),]
length(rownames(Volg_Mutselection[order(Volg_Mutselection[,1]),])=='KRAS')

##################################################################
# Test stable cas9 extent in WT vs Mut for each of the cancer gene
##################################################################
# Test cas9 activity by KRAS
CCLE_cellLineName=sapply(as.character(onTarget$annotation$CCLE_ID),
       function(x) strsplit(x, '_')[[1]][1])
matched_mutmatrix=onTarget$mutations_matrix[,match(onTarget$annotation$depMapID, colnames(onTarget$mutations_matrix))]
load('/Users/sinhas8/Project_CRISPR/Project_CRISPR.Risk_Github/Data/vogelstein.RData')
cas9Activity=cbind(cas9Activity, 
                   t(matched_mutmatrix[match(names(CDE), rownames(matched_mutmatrix)),
                                              match(cas9Activity$`Official Cell Line Name`,
                                                            CCLE_cellLineName)]))

cas9Activity
Volg_DifferentialCas9=t(apply(cas9Activity[,-(1:3)], 2, 
                             function(x) c(p_value=err_handle(wilcox.test(as.numeric(cas9Activity$`%GFP positive-cells`) ~ x,
                                                                alternative='g')$p.value),
                                           Effect_size=err_handle(median(100-as.numeric(cas9Activity$`%GFP positive-cells`)[which(x==1)], na.rm=T)-
                                                           median(100-as.numeric(cas9Activity$`%GFP positive-cells`)[which(x==0)], na.rm=T)))))
summary(lm(as.numeric(cas9Activity$`%GFP positive-cells`) ~ cas9Activity$KRAS+cas9Activity$TP53))$coefficients

median(as.numeric(cas9Activity$`%GFP positive-cells`)[cas9Activity$`TP53 status`=='mut'])
median(as.numeric(cas9Activity$`%GFP positive-cells`)[cas9Activity$`TP53 status`=='WT'])
##################################################################
# Test enrichment
##################################################################
mr.by.cde=readRDS('/Users/sinhas8/Downloads/mr.by.cde.RDS')
mr.by.de.in.crispr.RDS=readRDS('/Users/sinhas8/Downloads/mr.by.de.in.crispr.RDS')
mr.by.pval=readRDS('/Users/sinhas8/Downloads/mr.by.pval.RDS')
topK=10
hypergeometric_test_for_twolists(tail(names(sort(mr.by.pval)), topK), 
                                 head(names(Selection_sig_ALLvogl), topK),
                                 global = names(Selection_sig_ALLvogl))
hypergeometric_test_for_twolists(tail(names(sort(mr.by.cde)), topK), 
                                 head(names(Selection_sig_ALLvogl), topK),
                                 global = names(Selection_sig_ALLvogl))
hypergeometric_test_for_twolists(tail(names(sort(mr.by.de.in.crispr.RDS)), topK), 
                                 head(names(Selection_sig_ALLvogl), topK),
                                 global = names(Selection_sig_ALLvogl))

topK=10
hypergeometric_test_for_twolists(tail(names(sort(mr.by.pval)), topK), 
                                 head(names(sig_forall_cancergenes), topK),
                                 global = names(Selection_sig_ALLvogl))
hypergeometric_test_for_twolists(tail(names(sort(mr.by.cde)), topK), 
                                 head(names(sig_forall_cancergenes), topK),
                                 global = names(Selection_sig_ALLvogl))
hypergeometric_test_for_twolists(tail(names(sort(mr.by.de.in.crispr.RDS)), topK), 
                                 head(names(sig_forall_cancergenes), topK),
                                 global = names(Selection_sig_ALLvogl))
##################################################################
# Test enrichment
##################################################################
sort(Selection_sig_ALLvogl[match(vogelstein$symbol, names(Selection_sig_ALLvogl))])
df2plot=data.frame(MutSelection=Selection_sig_ALLvogl[match(vogelstein$symbol, names(Selection_sig_ALLvogl))],
                   bycas9=sig_forall_cancergenes[match(vogelstein$symbol, names(sig_forall_cancergenes))],
                   mr.by.pval=mr.by.pval[match(vogelstein$symbol, names(mr.by.pval))],
                   mr.by.cde=mr.by.cde[match(vogelstein$symbol, names(mr.by.cde))],
                   mr.by.de.in.crispr.RDS=mr.by.de.in.crispr.RDS[match(vogelstein$symbol, names(mr.by.de.in.crispr.RDS))]
                   )
rownames(df2plot)=vogelstein$symbol
df2plot=na.omit(df2plot)
sapply(1:ncol(df2plot), function(x) cor.test(df2plot[,2], df2plot[,x])[c(3, 4)] )
rownames(df2plot)[order(df2plot$mr.by.cde, decreasing = T)]


onTarget$annotation