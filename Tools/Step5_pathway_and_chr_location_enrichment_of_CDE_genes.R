# Version for Github
source('Step0_Globally_used_Functions_and_Datasets.R')

##### This script reproduces the results on the pathway and chromosomal location enrichments of the p53/KRAS/VHL CDE+/- genes

# load necessary files
cdes <- readRDS("../Data/CDE_for_volg.RDS") # list of CDE+/- genes
bg.genes <- readRDS("../Data/Background_Genes.RDS") # all the genes present in the data, used as background gene set for enrichment tests
load("../Data/reactome.RData") # the Reactome pathway/gene set data


####################################
######## Required Functions ########
####################################

make.confus.mat <- function(qset, refset, uset) {
  # create a contingency table/cofusion matrix
  uset <- unique(uset)
  qsetl <- factor(uset %in% qset, levels=c(TRUE, FALSE), labels=c("Positive","Negative"))
  refsetl <- factor(uset %in% refset, levels=c(TRUE, FALSE), labels=c("Positive","Negative"))
  table(`Query/Prediction`=qsetl, `Reference/Actual`=refsetl)
}

enrich.test <- function(qset, refset, uset) {
  # hypergeometric test
  conf <- make.confus.mat(qset, refset, uset)
  res <- fisher.test(conf, alternative="greater")
  res$table <- addmargins(conf)
  res
}

enrich.gset0 <- function(fg, gset, bg) {
  # test the enrichment of a set of fg genes for a given gene set (gset), with bg as the background gene set
  res <- enrich.test(qset=fg, refset=gset, uset=bg)
  data.table(overlap.gene.cnt=res$table[1,1], gset.gene.cnt=res$table[3,1], odds.ratio=res$estimate, p=res$p.value)
}

enrich.gsets <- function(fg, gsets=reactome, bg=bg.genes, padj.cutoff=0.1) {
  # a wrapper function around enrich.gset0 for testing for multiple gene sets and do FDR correction
  res <- lapply(gsets, enrich.gset0, fg=fg, bg=bg)
  res <- rbindlist(res, idcol="name")
  ncde <- length(fg)
  res <- res[overlap.gene.cnt>=ifelse(ncde>50, 5, ncde/10)]
  res[, p.adj:=p.adjust(p, method="BH")]
  res[order(p)][p.adj<=padj.cutoff]
}

inv.norm <- function(vec) {
  # inverse normal transformation
  qnorm(qrank(vec))
}

trans4m <- function(dat, method) {
  # transform a matrix
  res <- apply(dat, 1, method)
  res <- t(res)
  if (is.matrix(dat)) {
    rownames(res) <- rownames(dat)
    colnames(res) <- colnames(dat)
  } else if (is.data.frame(dat)) {
    rownames(res) <- row.names(dat)
    colnames(res) <- names(dat)
  }
  res
}

wilt <- function(x, group) {
  # wilcoxon test for differential essentiality
  c(p.pos=wilcox.test(x~group, alternative="less")$p.value,
    p.neg=wilcox.test(x~group, alternative="greater")$p.value)
}

p.dess <- function(gn) {
  # compute P values for differential essentiality, for each gene in both CRISPR and shRNA screens
  grp <- factor(mut.info[gn,]!=0)
  sh <- apply(rnai, 1, wilt, group=grp)
  cr <- apply(crispr, 1, wilt, group=grp)
  p <- as.data.table(t(rbind(cr, sh)), keep.rownames=TRUE)
  setnames(p, c("gene", "p.cr.pos","p.cr.neg","p.sh.pos","p.sh.neg"))
  p[, c("padj.cr.pos","padj.cr.neg","padj.sh.pos","padj.sh.neg"):=lapply(.SD, p.adjust, method="BH"), .SDcols=c("p.cr.pos","p.cr.neg","p.sh.pos","p.sh.neg")]
}

dess <- function(gn, func=median) {
  # compute differential essentiality scores, for each gene in both CRISPR and shRNA screens
  x <- mut.info[gn,]
  dgi.cr <- apply(crispr.norm[,x!=0],1,func) - apply(crispr.norm[,x==0],1,func) # crispr, direction as DE+ (i.e. the more positive, the more DE+)
  dgi.sh <- apply(rnai.norm[,x!=0],1,func) - apply(rnai.norm[,x==0],1,func) # rnai, ditto
  ddgi <- dgi.cr - dgi.sh
  data.table(gene=names(dgi.cr), dess.cr=dgi.cr, dess.sh=dgi.sh, ddess=ddgi)
}

gsea0 <- function(dat, pathways=reactome, minSize=10, nperm=10000) {
  # function for GSEA analysis
  set.seed(1)
  res <- fgsea(pathways=pathways, stats=dat, minSize=minSize, nperm=nperm)
  res[, leadingEdge:=sapply(leadingEdge, paste0, collapse=",")]
  res[order(padj, pval)]
}

gsea <- function(mr) {
  # a wrapper function around gsea0 to obtain significantly enriched pathways
  gns <- dess.p[[mr]][padj.sh.pos>=0.1 & padj.sh.neg>=0.1][, gene]
  gvec <- dess[[mr]][gene %in% gns, .(gene, dess.cr)]
  tmp <- gvec$gene
  gvec <- gvec$dess.cr
  names(gvec) <- tmp
  cr <- gsea0(gvec)
  gns <- dess.p[[mr]][padj.cr.pos>=0.1 & padj.cr.neg>=0.1][, gene]
  gvec <- dess[[mr]][gene %in% gns, .(gene, dess.sh)]
  tmp <- gvec$gene
  gvec <- gvec$dess.sh
  names(gvec) <- tmp
  sh <- gsea0(gvec)
  cr[(ES>0 & padj<0.1 & !pathway %in% sh[ES>0 & padj<0.1,pathway]) | (ES<0 & padj<0.1 & !pathway %in% sh[ES<0 & padj<0.1,pathway])]
}


################################################################
######## 1. Pathway enrichment analysis of CDE+/- genes ########
################################################################

## Method 1: based on hypergeometric test
p53.enr <- list(CDE_pos=enrich.gsets(cdes$TP53[[1]]), CDE_neg=enrich.gsets(cdes$TP53[[2]]))
kras.enr <- list(CDE_pos=enrich.gsets(cdes$KRAS[[1]]), CDE_neg=enrich.gsets(cdes$KRAS[[2]]))
vhl.enr <- list(CDE_pos=enrich.gsets(cdes$VHL[[1]]), CDE_neg=enrich.gsets(cdes$VHL[[2]]))

## Method 2: based on GSEA

# these are the code to compute the differential essentiality scores and P values, used later for GSEA analysis; the results are contained in the objects `dess` and `dess.p`, respectively, and are saved as ../Data/de_score.RData and ../Data/de_pval.RData, these files have been pre-loaded.
#rnai <- achilles
#crispr <- avana
## the above two are in the same order
#mut.info <- Mut_CCLE
#which((names(rnai)!=colnames(mut.info))) # 1 2 3, I checked that they are actually the same cell lines
#colnames(mut.info)[1:3] <- names(rnai)[1:3]
## normalize data
#rnai.norm <- trans4m(rnai, inv.norm)
#crispr.norm <- trans4m(crispr, inv.norm)
## results: differential essentiality scores
#dess <- list(TP53=dess("TP53"), VHL=dess("VHL"), KRAS=dess("KRAS"))
#save(dess, file="../Data/de_score.RData")
## results: P value of differential essentiality
#dess.p <- list(TP53=p.dess("TP53"), VHL=p.dess("VHL"), KRAS=p.dess("KRAS"))
#save(dess.p, file="../Data/de_pval.RData")
load("../Data/de_score.RData")
load("../Data/de_pval.RData")

# run GSEA
p53.gsea <- gsea("TP53")
kras.gsea <- gsea("KRAS")
vhl.gsea <- gsea("VHL")


######################################################################################
# 2. The enrichment of p53/VHL CDE+ genes in common chromosomal fragile sites (CFSs) #
######################################################################################

# the lines below are used to obtain the chromosomal bands of all the genes in the data (`bands` saved in ../Data/chr_bands.RData), loaded below
#library(biomaRt) # need to have the biomaRt package installed from bioconductor
#mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#bands <- getBM(attributes=c("hgnc_symbol", "chromosome_name", "band"), filters="hgnc_symbol", values=bg.genes, mart=mart)
#bands <- as.data.table(bands)
#bands[, chr:=paste0(chromosome_name, band)]
#bands[, chr.arm:=paste0(chromosome_name, stringr::str_extract(band, "p|q"))]
load("../Data/chr_bands.RData")
# chromosomal bands of CFSs, obtained from Lukusa and Fryns, Biochim Biophys Acta, 2008 and HumCFS database
load("../Data/cfs.RData")
cfs <- fragile.sites[frequency=="common", unique(band)]
cfs.genes <- list(cfs=bands[chr %in% cfs, unique(hgnc_symbol)])

# result for p53 CDE+ genes
p53.cfs.res <- enrich.gsets(fg=cdes$TP53[[1]], gsets=cfs.genes)
# result for VHL CDE+ genes
vhl.cfs.res <- enrich.gsets(fg=cdes$VHL[[1]], gsets=cfs.genes)


#############################################################################
####### 3. The enrichment of p53 CDE+ genes in open chromatin regions #######
#############################################################################
# Functions Required
# Method- Via sgRNA target position overlapping with open chromatin regions
# Function which identifies whether a gene is present in a Open Chromatin region
whetherin_OpenCHR_sgRNA<-function(genomic_pos=sgRNA[1,6], chr=sgRNA[1,5]){
  OC_in_respective_Chr=openChr_list[[chr]]
  genomic_pos=as.numeric(genomic_pos)
  sum(genomic_pos>OC_in_respective_Chr$V2 & genomic_pos<OC_in_respective_Chr$V3)
}
CDE_allgenes=readRDS('../Data/CDE_for_volg.RDS')
p53_CDE_pos=CDE_allgenes$TP53[[1]]
openChr=read.csv('../Data/Ubiquitous.bed', sep='\t', header=F)
geneAnnotation=read.csv('../Data/mart_export.txt', sep='\t')
geneAnnotation=geneAnnotation[match(rownames(avana), geneAnnotation$Gene.name),]
global_list=rownames(avana)[!is.na(match(rownames(avana), geneAnnotation$Gene.name))]
geneAnnotation=na.omit(geneAnnotation)
openCHR_counts=mclapply(1:nrow(geneAnnotation), function(x) 
  err_handle(whetherin_OpenCHR_gene(geneAnnotation[x,])))
openCHR_counts=unlist(openCHR_counts)

sgRNA=read.csv('../Data/sgRNA_with_GeneMapping.csv')
sgRNA$chr = sapply(sgRNA$genome_alignment, function(x) strsplit(as.character(x), '_')[[1]][1])
sgRNA$Position = sapply(sgRNA$genome_alignment, function(x) strsplit(as.character(x), '_')[[1]][2])
#break by CHR
sgRNA_list=split(sgRNA, sgRNA$chr)
openChr_list=split(openChr, openChr$V1)
openChr_list=openChr_list[match(names(sgRNA_list),names(openChr_list))]
OpenRegions_Score=apply(sgRNA, 1, function(x) whetherin_OpenCHR_sgRNA(genomic_pos=x[6], chr=x[5]))
sgRNA=cbind(sgRNA, OpenRegions_Score)
geneLevel_Score=aggregate(sgRNA$OpenRegions_Score~sgRNA$gene, data=sgRNA, function(x) sum(x))
genes_withatleastonesgRNAinopenCHR=geneLevel_Score[which(geneLevel_Score$`sgRNA$OpenRegions_Score`>0),]
genes_withatleastonesgRNAinopenCHR=genes_withatleastonesgRNAinopenCHR[order(genes_withatleastonesgRNAinopenCHR$`sgRNA$OpenRegions_Score`, decreasing = T),]
genes_withatleastonesgRNAinopenCHR$GeneName=sapply(genes_withatleastonesgRNAinopenCHR$`sgRNA$gene`,
                                                   function(x) strsplit(as.character(x), ' ')[[1]][1])

##Hypergeometric test for the enrichment for CDE+ genes in top sgRNAs in open Regions
hits_in_openRegion_Threshold=1
hypergeometric_test_for_twolists(test_list = p53_CDE_pos, #!!! p53_CDE_pos not found
                                 base_list = unique(genes_withatleastonesgRNAinopenCHR$GeneName[
                                   genes_withatleastonesgRNAinopenCHR$`sgRNA$OpenRegions_Score`>hits_in_openRegion_Threshold]),
                                 global=global_list)
## Here, we have varied the above threshold from 1:3 to identify test the robustness of the results
hits_in_openRegion_Threshold=2
hypergeometric_test_for_twolists(test_list = p53_CDE_pos, #!!! p53_CDE_pos not found
                                 base_list = unique(genes_withatleastonesgRNAinopenCHR$GeneName[
                                   genes_withatleastonesgRNAinopenCHR$`sgRNA$OpenRegions_Score`>hits_in_openRegion_Threshold]),
                                 global=global_list)
hits_in_openRegion_Threshold=3
hypergeometric_test_for_twolists(test_list = p53_CDE_pos, #!!! p53_CDE_pos not found
                                 base_list = unique(genes_withatleastonesgRNAinopenCHR$GeneName[
                                   genes_withatleastonesgRNAinopenCHR$`sgRNA$OpenRegions_Score`>hits_in_openRegion_Threshold]),
                                 global=global_list)

