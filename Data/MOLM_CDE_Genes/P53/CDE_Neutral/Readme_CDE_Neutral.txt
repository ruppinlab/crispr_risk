##ReadME for both the CDE_Neutral files::

Herewith we have provided the methods to curate, MOLM13 specific CDE_Neutral.

**Definition**

1. CDE_Neutral: Genes which doesn't show a significant differential essentiality in any of the two screens, CRISPR or shRNA.

2. Complementary of a Set A: (Universal - Set A),
which is denoted by  the first four initials "Comp".

Method:
For a set of mutated and wild type cell lines, CDE_Neutral genes are curated by taking an intersection of Comp(DE+ and DE-) of CRISPR screenings and Comp(DE+ and DE-) of shRNA screenings.

*S3 table contains a list curated by taking a intersection of CDE_neutral genes found using control (i), (ii) and manuscript files.


NOTE: CDE_Neutral is ordered in increasing order as we want our top hits to show no differential essentaility.
