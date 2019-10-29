##ReadME for both the CDE_Negative files::

Our two objectives are following:
1. To test the robustness of the p53 CDE effect in MOLM13 and THP1 cell lines, respectively.
2. Provide ranked lists of CDE+ genes for MOLM13 and THP1 cell lines, respectively. 

Analysis for MOLM13 cell line:

To test the robustness, we repeated the same analysis as in the manuscript focusing on MOLM13 cell line. Since we do not have p53-mutant MOLM13 cell line in the CCLE collection, we used two different controls: (i) p53-mutant haematopoietic cell lines (n=10) and (ii) all p53-mutant cell lines (n=173). For both cases (i) and (ii), we find the CRISPR-KO has significantly larger number of DE+ while shRNA-KD is balanced (Chi-square test P<4.8E-66 and P<3.0E-43 in respective order), which confirms the robustness of our selection for p53-mutant.  
Ranking CDE+ genes. (i) Using p53-mutant haematopoietic cell lines as control, we identified 756 CDE+ genes, and these genes are significantly enriched with our original CDE+ genes in the manuscript (Hypergeometric p-value<1.7E-12; Overlap: 94 genes). (ii) Using all p53-mutant cell lines as control, we identified 1970 CDE+ genes are enriched with our original CDE+ genes (Hypergeometric p-value< 6.4E-57; Overlap: 291 genes). We took a union of the above two CDE+ genesets, resulting in 305 genes, and ranked them in the order of median essentiality difference between the MOLM13 and p53-mutant haematopoietic cell lines (i.e. based on (i)).
 
Above analysis is repeated for THP1 cell line, where p53 is mutated.

Testing the robustness. (i) Comparing THP1 with haematopoietic cell lines with WT p53 (n=9), no genes are significantly differentially essential with respect to p53 mutation status in the CRISPR-KO screens, and very few genes for shRNA screens. (ii) Comparing THP1 profile to all of the p53 wild type cell lines of different cancer types, we find the CRISPR-KO has significantly larger number of DE+ while shRNA-KD is balanced (Chi-square test P<8.15E-17). CDE+ identified from this analysis is also enriched for our previously identified CDE+ genes. (Hypergeometric p<2.1E-43; overlap)


Based on the above basic analysis, we suggest the use of MOLM13 over THP1. We have also provided a list of 305 genes presented in the order of median essentiality difference in supplementary table 1.(S1)
