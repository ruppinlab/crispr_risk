# A systematic genome-wide mapping of oncogenic mutation selection during CRISPR-Cas9 genome editing

This repository contains code used in our study on the potential selection of specific cancer drivers mutation during CRISPR-Cas9 gene editing.

The Data folder contains the data needed for analysis.
The Tools folder contains various R scripts for performing the analysis (file names starting with "Step") and plotting the results (file names starting with "Plot").
Specifically:

* Step0_Globally_used_Functions_and_Datasets.R: defining some functions used for the analysis

* Step1_identify_master_regulators_and_CDE_genes.R: identification of the CDE+/- genes, as well as the potential CRISPR-selected cancer genes (CCDs)

* Step2_emperical_p_value_calculations.R: computing empricial P value for the CCD identification analysis

* Step3_p53_analysis_with_only_LOF_mutations.R: confirming the robustness of the p53 results, considering only p53 LOF mutations

* Step4_compute_mean_risk_score.R: computing the score representing the effect size for mutant selection for each CCD gene

* Step5_pathway_and_chr_location_enrichment_of_CDE_genes.R: pathway and chromosomal location enrichment analysis of CDE+/- genes

* Step6_CDE_genes_for_MOLM13_cell_line.R: CDE+/- genes identification specifically for the MOLM13 cell line

* Step7A_preprocess_data_for_CRISPR_screen_in_MOLM13.R: preprocessing the CRISPR-KO/CRISPRi genetic screen data in MOLM13 cell line

* Step7B_analyze_CRISPR_screen_data_in_MOLM13_p53.R: nalyzing the CRISPR-KO/CRISPRi genetic screen data in p53-isogenic MOLM13 cell line

* Step7C_analyze_CRISPR_screen_data_in_MOLM13_KRAS.R: analyzing the CRISPR-KO/CRISPRi genetic screen data in KRAS-isogenic MOLM13 cell line

* Step8_p53_analysis_Haapemeini.R: re-analysis of the genetic screen data from Haapemeini et al.

* Step9_analyze_published_CRISPR_screens.R: analysis of published CRISPR screen in p53 or KRAS-isogenic cells

* Step10_Cas9vsdriverStatus.R: analysis of Cas9 activity difference by the mutation status of cancer genes

* Step11_addressing_penetrance.R: addressing the potential confounding factor of the different penetrance associated with CRISPR-KO and RNAi

