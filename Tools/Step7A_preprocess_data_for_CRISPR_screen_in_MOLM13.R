# Version for Github
source('Step0_Globally_used_Functions_and_Datasets.R')

# This script prepares our CRISPR screening data in isogenic
# MOLM13 cell lines for downstream analyses with the scripts 
# Step7B and Step7C
####################################################################################
#### CRSIPRi-Screenings Preprocessing
####################################################################################
ip53=read.csv('../Data/Isogenic_Screenings/ip53_karina.csv')

ip53_Mut=ip53[c(1:5, 9:10, 15:16, 21:22)]
ip53_WT=ip53[c(1:5, 9:10, 15:16, 21:22)+24]
colnames(ip53_Mut)=as.character(unlist(ip53_Mut[1,]))
colnames(ip53_Mut)[c(3:4, 6:11)]=c('Count_D0_R1', 'cpm_D0_R1', 'Count_D0_R2',  'cpm_D0_R2', 'Count_D30_R1',  'cpm_D30_R1', 'Count_D30_R2',  'cpm_D30_R2' )
colnames(ip53_WT)=colnames(ip53_Mut)
ip53_WT=ip53_WT[-1,]
ip53_Mut=ip53_Mut[-1,]

write.csv(ip53_WT, '../Data/Isogenic_Screenings/ip53_WT.csv')
write.csv(ip53_Mut, '../Data/Isogenic_Screenings/ip53_Mut.csv')


####################################################################################
#### Processing the sgRNA sequence data and compute off-target scores
####################################################################################

# *PLEASE READ* - This part could take 6-7 days to finish with a CPUs=64, Mem=64G #
# The results are saved as sgRNA_Off_Target_Scores.RDS in the Data folder, and will be loaded for the downstream analysis
# Write sgRNA sequences in FASTA
#sgRNA=read.csv('../Data/sgRNA_with_GeneMapping.csv')
#sgRNA_list=mclapply(sgRNA$sgrna,  
#                    function(x) paste(as.character(x), 'NGG', sep=''), mc.cores = detectCores())

#K=length(sgRNA_list)
#write.fasta(sequences = sgRNA_list[1:K], names = make.names(sgRNA$gene, unique = T)[1:K],
#            file.out = "../Data/sgRNA_sequences_avana_foroffTarget_scores.fa")

#if (!requirenamespace("biocmanager", quietly = true))
#  install.packages("biocmanager")
#biocmanager::install("crisprseek")
#biocmanager::install("bsgenome.hsapiens.ucsc.hg19")
#biocmanager::install("txdb.hsapiens.ucsc.hg19.knowngene")
#biocmanager::install("org.hs.eg.db")
#lib_usedhere=c('crisprseek', 'bsgenome.hsapiens.ucsc.hg19', 'txdb.hsapiens.ucsc.hg19.knowngene', 'org.hs.eg.db')
#installorload(lib_usedhere)

#outputdir <- "../Data"
#inputfilepath <- '../Data/sgrna_sequences_avana_forofftarget_scores.fa'
#repatternfile <- system.file('extdata', 'nebenzymes.fa', package = 'crisprseek')
#grnafilepath <- system.file('extdata', 'testhsap_gata1_ex2_grna1.fa', package = 'crisprseek')
# please check the manual for arguements reference
#results1 <- offtargetanalysis(inputfilepath = inputfilepath,
#                             findgrnaswithrecutonly = false,
#                             #repatternfile = repatternfile,
#                             findpairedgrnaonly = false,
#                             findgrnas = false,
#                             bsgenomename = hsapiens,
#                             #chromtosearch = 'chrx',
#                             txdb = txdb.hsapiens.ucsc.hg19.knowngene,
#                             organn = org.hs.egsymbol,
#                             max.mismatch = 3, 
#                             outputdir = outputdir, 
#                             overwrite = true,
#                             enable.multicore=true,
#                             n.cores.max=detectcores(),
#                             annotateexon = t,
#                             topn.offtargettotalscore=100,
#                             fetchsequence=f,
#                             scoring.method='cfdscore')
#saveRDS(results1$summary, '../Data/sgrna_off_target_scores.RDS')
