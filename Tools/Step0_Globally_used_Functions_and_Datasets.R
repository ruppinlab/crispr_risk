########################################################################
###My Commonly Used functions
########################################################################
# Below function takes an input of list of packages and loads a Package,
# if available or Installs&loads, if not.
installORload<-function(packages){
  package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
}
installORload.bioc<-function(packages){
  package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x)
      library(x, character.only = TRUE)
    }
  })
}
packages_required<-c('parallel', 'grid', 'ggplot2', 'gridExtra',
                     'cowplot', 'ggpubr', 'ggrepel', 'readxl',
                     'reshape2','data.table', 'statar', 'tidyr', 'Rmisc',
                     'seqinr', 'matrixStats','tictoc','BiocManager','rowr')
installORload(packages_required)
installORload.bioc("fgsea")

##Defining not included
'%!in%' <- function(x,y)!('%in%'(x,y))
myhead<-function(x){
  x[1:min(5, nrow(x)), 1:min(5, ncol(x))]
}
err_handle<-function(x){ tryCatch(x, error=function(e){NA}) }
hypergeometric_test_for_twolists<-function(test_list, base_list, global, lowertail=FALSE) {
  #If lowertail=FALSE - we calculate the probability for enrichment
  length(base_list)
  adj_base_list=global[na.omit(match(base_list, global))]
  Matched_list=test_list[!is.na(match(test_list, adj_base_list))]
  phyper(length(Matched_list)-1, length(adj_base_list),
         length(global)- length(adj_base_list), length(test_list), lower.tail=lowertail)
}
fdrcorr<-function(test_list){p.adjust(test_list, method = 'fdr')}
#Subsetting a set of columns
colSubset<-function(mat, column_Names){
  mat[,na.omit(match(column_Names, colnames(mat)))]
}
#Subsetting a set of rows
rowSubset<-function(mat, row_Names){
  mat[na.omit(match(row_Names, rownames(mat))),]
}
########################################################################
###Project-specific functions
########################################################################
##This function return four a list of DE+/- for a given master regulator in shrna and crispr
Testing_CRISPR_damage_bias<-function(GeneName='TP53',
                                     Feature_to_test=Mut_CCLE, 
                                     base_dataset=avana, 
                                     reference_dataset=achilles,
                                     confounding_variable_dataset=CNV, 
                                     cores2use=3){
	##Creating MR contigency
	print(paste('The Candidate gene is ', GeneName))
	CRISPR_left =mclapply(1:nrow(base_dataset), function(x) 
	  wilcox.test(unlist(base_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]),
	              alternative='l' )$p.value, mc.cores=cores2use)
	GeneCount_inCRISPR_pos_affected =p.adjust(CRISPR_left, method='fdr')
	print(paste("Step 1 Completed", sum(GeneCount_inCRISPR_pos_affected<0.1) ))

	CRISPR_right =mclapply(1:nrow(base_dataset), function(x)
	  wilcox.test(unlist(base_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]),
	              alternative='g' )$p.value, mc.cores=cores2use)
	GeneCount_inCRISPR_neg_affected =p.adjust(CRISPR_right, method='fdr')
	print(paste("Step 2 Completed", sum(GeneCount_inCRISPR_neg_affected<0.1) ))

	shrna_left=mclapply(1:nrow(reference_dataset), function(x) 
	  wilcox.test(unlist(reference_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]), 
	              alternative='l' )$p.value, mc.cores=cores2use)
	GeneCount_inshrna_pos_affected =p.adjust(shrna_left, method='fdr')
	print(paste("Step 3 Completed", sum(GeneCount_inshrna_pos_affected<0.1) ))

  shrna_right=mclapply(1:nrow(reference_dataset), function(x)
    wilcox.test(unlist(reference_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]), 
                alternative='g' )$p.value, mc.cores=cores2use)
	GeneCount_inshrna_neg_affected =p.adjust(shrna_right, method='fdr')
	print(paste("Step 4 Completed. ", sum(GeneCount_inshrna_neg_affected<0.1) ))

	##Controlling for Copy Number Effect
	CNV_association_sig =mclapply(1:nrow(confounding_variable_dataset),
	                              function(x) tryCatch(wilcox.test(unlist(confounding_variable_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]))$p.value, 
	                                                   error=function(err){NA}), mc.cores=detectCores())
	FDR_CNV_association_sig =p.adjust(CNV_association_sig, method='fdr')
	print(paste("Step 1 Completed", sum(FDR_CNV_association_sig<0.1, na.rm=T) ))
	
	#Removing the genes whose CNV and Ess are assocaited.
	CRISPR_neg_id=which(GeneCount_inCRISPR_neg_affected<0.1)[which(GeneCount_inCRISPR_neg_affected<0.1) %!in% which(FDR_CNV_association_sig<0.1)]
	CRISPR_pos_id=which(GeneCount_inCRISPR_pos_affected<0.1)[which(GeneCount_inCRISPR_pos_affected<0.1) %!in% which(FDR_CNV_association_sig<0.1)]
	shrna_neg_id =which(GeneCount_inshrna_neg_affected<0.1)
	shrna_pos_id =which(GeneCount_inshrna_pos_affected<0.1)
	
	#Processing and variables to return a df
	CRISPR_neg=rownames(avana)[CRISPR_neg_id]
	CRISPR_pos=rownames(avana)[CRISPR_pos_id]
	shrna_neg=rownames(avana)[shrna_neg_id]
	shrna_pos=rownames(avana)[shrna_pos_id]
	
	#Number of DE+ DE- in CRISPR
	print(list(length(CRISPR_neg), length(CRISPR_pos)))
	#df of DE+ DE- in CRISPR and shRNA with their respective Effect Size
	CRISPR_DE_neg=data.frame(GeneName=CRISPR_neg,
	                         Effect_Size=GeneCount_inCRISPR_neg_affected[CRISPR_neg_id])
	CRISPR_DE_pos=data.frame(GeneName=CRISPR_pos, 
	                         Effect_Size=GeneCount_inCRISPR_pos_affected[CRISPR_pos_id])
	shrna_DE_neg=data.frame(GeneName=shrna_neg, 
	                        Effect_Size=GeneCount_inshrna_neg_affected[shrna_neg_id])
	shrna_DE_pos=data.frame(GeneName=shrna_pos, 
	                        Effect_Size=GeneCount_inshrna_pos_affected[shrna_pos_id])

	#Variable to return
	list(CRISPR_DE_neg=CRISPR_DE_neg[order(CRISPR_DE_neg$Effect_Size,decreasing = F),],
	     CRISPR_DE_pos=CRISPR_DE_pos[order(CRISPR_DE_pos$Effect_Size,decreasing = F),],
	     shrna_DE_neg=shrna_DE_neg[order(shrna_DE_neg$Effect_Size,decreasing = F),],
	     shrna_DE_pos=shrna_DE_pos[order(shrna_DE_pos$Effect_Size,decreasing = F),])
}
#Function to get CDE+/- genes in an ordered form corrected for Copy number
CDE_Genes<- function(Gene_Contigency){
	CDE_pos=Gene_Contigency[[1]][Gene_Contigency[[1]][,1]  %!in% Gene_Contigency[[3]][,1],]
	CDE_neg=Gene_Contigency[[2]][Gene_Contigency[[2]][,1]  %!in% Gene_Contigency[[4]][,1],]
	list(CDE_pos[order(CDE_pos[,2]),1], CDE_neg[order(CDE_neg[,2]),1])
}
#func2#Our aim is to remove the genes whose mutation profile is assocaited with TP53
Mut_Assoc<-function(gene1='TP53', gene2){
	none=sum(!as.logical(Mut_CCLE[gene1,]) & !as.logical(Mut_CCLE[gene2,]) )
	base_only=sum(as.logical(Mut_CCLE[gene1,]) & !as.logical(Mut_CCLE[gene2,]) )
	reference_only=sum(!as.logical(Mut_CCLE[gene1,]) & as.logical(Mut_CCLE[gene2,]) )
	together=sum(as.logical(Mut_CCLE[gene1,]) & as.logical(Mut_CCLE[gene2,]) )
	fisher.test(matrix(c(none, base_only, reference_only, together),2,2), alternative='g')$p.value	
}
#func1#Significance test## 30 events minimun to choose a Chisq test
significance_test<-function(x){
	##We assume that if all the counts>30, chisq would be approximately close to Fisher's test.
  if(x[1]> 30 & x[2] >30 & x[3]>30 & x[4]>30)
		err_handle(chisq.test(matrix(c(x[1],x[2]
		                             ,nrow(avana)-(x[1]+x[2]),x[3],x[4],
		                             nrow(avana)-(x[3]+x[4])), 3, 2))[[3]])
	else
		# tryCatch(fisher.test(matrix(c(x[1],x[2],
		#                               nrow(avana)-(x[1]+x[2]),x[3],x[4],
		#                               nrow(avana)-(x[3]+x[4])), 3, 2),
		#                      simulate.p.value=TRUE, B=10000)[[1]][1], error=function(err){NA})
	  err_handle(fisher.test(matrix(c(x[1],x[2],
	                                nrow(avana)-(x[1]+x[2]),
	                                x[3],x[4],
	                                nrow(avana)-(x[3]+x[4])), 3, 2))$p.value)
}
##DE+/- Complete set for an enrichment
for_GSEA<-function(GeneName='TP53', Feature_to_test=Mut_CCLE,
                   base_dataset=avana, reference_dataset=achilles, 
                   confounding_variable_dataset=CNV, cores2use=detectCores()){
	##Creating MR contigency
	print(paste('The Candidate gene is ', GeneName))
	CRISPR_left =mclapply(1:nrow(base_dataset), function(x) wilcox.test(unlist(base_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]), alternative='l' )$p.value, mc.cores=cores2use)
	GeneCount_inCRISPR_neg_affected =p.adjust(CRISPR_left, method='fdr')
	print(paste("Step 1 Completed", sum(GeneCount_inCRISPR_neg_affected<0.1) ))

	CRISPR_right =mclapply(1:nrow(base_dataset), function(x) wilcox.test(unlist(base_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]), alternative='g' )$p.value, mc.cores=cores2use)
	GeneCount_inCRISPR_pos_affected =p.adjust(CRISPR_right, method='fdr')
	print(paste("Step 2 Completed", sum(GeneCount_inCRISPR_pos_affected<0.1) ))

	shrna_left=mclapply(1:nrow(reference_dataset), function(x) wilcox.test(unlist(reference_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]), alternative='l' )$p.value, mc.cores=cores2use)
	GeneCount_inshrna_neg_affected =p.adjust(shrna_left, method='fdr')
	print(paste("Step 3 Completed", sum(GeneCount_inshrna_neg_affected<0.1) ))

       shrna_right=mclapply(1:nrow(reference_dataset), function(x) wilcox.test(unlist(reference_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]), alternative='g' )$p.value, mc.cores=cores2use)
	GeneCount_inshrna_pos_affected =p.adjust(shrna_right, method='fdr')
	print(paste("Step 4 Completed. ", sum(GeneCount_inshrna_pos_affected<0.1) ))

	##Controlling for Copy Number Effect
	CNV_association_sig =mclapply(1:nrow(confounding_variable_dataset), function(x) tryCatch(wilcox.test(unlist(confounding_variable_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]))$p.value, error=function(err){NA}), mc.cores=detectCores())
	FDR_CNV_association_sig =p.adjust(CNV_association_sig, method='fdr')
	print(paste("Step 1 Completed", sum(FDR_CNV_association_sig<0.1, na.rm=T) ))
		median_diff = apply(avana, 1, function(x) median(x[!as.logical(Mut_CCLE[GeneName,])]) - median(x[as.logical(Mut_CCLE[GeneName,])]) )
	if(length(union(which(FDR_CNV_association_sig<0.1), which(GeneCount_inshrna_pos_affected<0.1)))>0 ){
		DE_pos_ranked_Score =  median_diff[-union(which(FDR_CNV_association_sig<0.1), which(GeneCount_inshrna_pos_affected	<0.1))]
		DE_pos_ranked_GeneName = rownames(base_dataset)[-(union(which(FDR_CNV_association_sig<0.1), which(GeneCount_inshrna_pos_affected<0.1)))]
	}
	else{
		DE_pos_ranked_Score = median_diff
		DE_pos_ranked_GeneName = rownames(base_dataset)
	}
	DE_neg_ranked_Score= -median_diff[-union(which(FDR_CNV_association_sig<0.1), which(GeneCount_inshrna_neg_affected<0.1))]
	DE_neg_ranked_GeneName = rownames(base_dataset)[-union(which(FDR_CNV_association_sig<0.1), which(GeneCount_inshrna_neg_affected<0.1))]

	DE_pos=data.frame(DE_pos_ranked_GeneName, DE_pos_ranked_Score)
	DE_neg=data.frame(DE_neg_ranked_GeneName, DE_neg_ranked_Score)
	list(DE_pos, DE_neg)
}
Randomized_Testing_CRISPR_damage_bias<-function(GeneName='TP53',
                                     Feature_to_test=Mut_CCLE, 
                                     base_dataset=avana, 
                                     reference_dataset=achilles,
                                     confounding_variable_dataset=CNV, 
                                     cores2use=detectCores()){
  ##Creating MR contigency
  print('This tests shuffles the samples annotation/cols and provides p-value and
        could be used to provide emperical p-value')
  Feature_to_test = Feature_to_test[,sample(1:ncol(Feature_to_test))]
  ##Creating MR contigency
  print(paste('The Candidate gene is ', GeneName))
  CRISPR_left =mclapply(1:nrow(base_dataset), function(x) 
    wilcox.test(unlist(base_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]),
                alternative='l' )$p.value, mc.cores=cores2use)
  GeneCount_inCRISPR_pos_affected =p.adjust(CRISPR_left, method='fdr')
  print(paste("Step 1 Completed", sum(GeneCount_inCRISPR_pos_affected<0.1) ))
  
  CRISPR_right =mclapply(1:nrow(base_dataset), function(x)
    wilcox.test(unlist(base_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]),
                alternative='g' )$p.value, mc.cores=cores2use)
  GeneCount_inCRISPR_neg_affected =p.adjust(CRISPR_right, method='fdr')
  print(paste("Step 2 Completed", sum(GeneCount_inCRISPR_neg_affected<0.1) ))
  
  shrna_left=mclapply(1:nrow(reference_dataset), function(x) 
    wilcox.test(unlist(reference_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]), 
                alternative='l' )$p.value, mc.cores=cores2use)
  GeneCount_inshrna_pos_affected =p.adjust(shrna_left, method='fdr')
  print(paste("Step 3 Completed", sum(GeneCount_inshrna_pos_affected<0.1) ))
  
  shrna_right=mclapply(1:nrow(reference_dataset), function(x)
    wilcox.test(unlist(reference_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]), 
                alternative='g' )$p.value, mc.cores=cores2use)
  GeneCount_inshrna_neg_affected =p.adjust(shrna_right, method='fdr')
  print(paste("Step 4 Completed. ", sum(GeneCount_inshrna_neg_affected<0.1) ))
  
  ##Controlling for Copy Number Effect
  CNV_association_sig =mclapply(1:nrow(confounding_variable_dataset),
                                function(x) tryCatch(wilcox.test(unlist(confounding_variable_dataset[x,]) ~ as.logical(Feature_to_test[GeneName,]))$p.value, 
                                                     error=function(err){NA}), mc.cores=detectCores())
  FDR_CNV_association_sig =p.adjust(CNV_association_sig, method='fdr')
  print(paste("Step 1 Completed", sum(FDR_CNV_association_sig<0.1, na.rm=T) ))
  
  #Removing the genes whose CNV and Ess are assocaited.
  CRISPR_neg_id=which(GeneCount_inCRISPR_neg_affected<0.1)[which(GeneCount_inCRISPR_neg_affected<0.1) %!in% which(FDR_CNV_association_sig<0.1)]
  CRISPR_pos_id=which(GeneCount_inCRISPR_pos_affected<0.1)[which(GeneCount_inCRISPR_pos_affected<0.1) %!in% which(FDR_CNV_association_sig<0.1)]
  shrna_neg_id =which(GeneCount_inshrna_neg_affected<0.1)
  shrna_pos_id =which(GeneCount_inshrna_pos_affected<0.1)
  
  #Processing and variables to return a df
  CRISPR_neg=rownames(avana)[CRISPR_neg_id]
  CRISPR_pos=rownames(avana)[CRISPR_pos_id]
  shrna_neg=rownames(avana)[shrna_neg_id]
  shrna_pos=rownames(avana)[shrna_pos_id]
  
  #Number of DE+ DE- in CRISPR
  print(list(length(CRISPR_neg), length(CRISPR_pos)))
  #df of DE+ DE- in CRISPR and shRNA with their respective Effect Size
  CRISPR_DE_neg=data.frame(GeneName=CRISPR_neg,
                           Effect_Size=GeneCount_inCRISPR_neg_affected[CRISPR_neg_id])
  CRISPR_DE_pos=data.frame(GeneName=CRISPR_pos, 
                           Effect_Size=GeneCount_inCRISPR_pos_affected[CRISPR_pos_id])
  shrna_DE_neg=data.frame(GeneName=shrna_neg, 
                          Effect_Size=GeneCount_inshrna_neg_affected[shrna_neg_id])
  shrna_DE_pos=data.frame(GeneName=shrna_pos, 
                          Effect_Size=GeneCount_inshrna_pos_affected[shrna_pos_id])
  
  #Variable to return
  list(CRISPR_DE_neg=CRISPR_DE_neg[order(CRISPR_DE_neg$Effect_Size,decreasing = F),],
       CRISPR_DE_pos=CRISPR_DE_pos[order(CRISPR_DE_pos$Effect_Size,decreasing = F),],
       shrna_DE_neg=shrna_DE_neg[order(shrna_DE_neg$Effect_Size,decreasing = F),],
       shrna_DE_pos=shrna_DE_pos[order(shrna_DE_pos$Effect_Size,decreasing = F),])
}

########################################################################
###Files required in the Project
########################################################################
# loading matrices needed
# CRISPR screening version : https://depmap.org/portal/download/
avana=readRDS('../Data/avana.RDS')
# avana screening : https://depmap.org/portal/download/
achilles=readRDS('../Data/achilles.RDS')
# Mutation profiles of respective cell lines
Mut_CCLE=readRDS('../Data/Mut_CCLE.RDS')
# gene-level CNV profile
CNV=readRDS('../Data/CNV.RDS')
# Expression profile
Exp=readRDS('../Data/Exp.RDS')
