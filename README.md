# nilas_haploimputer
## Overview of nilas_haploimputer:

NILAS_haploimputer was written in Python 3.6; data frames constructed using Pandas v0.19.2 (McKinney 2010). NILAS_haploimputer selects the donor and recurrent parent genotype matrices from consensus matrix and then parses the corresponding NILAS lines from the genotype extracted VCF; concatenating the founder and NILAS progeny into one genotype matrix in order to:  i) filter heterozygous and non-polymorphic markers from donor/recurrent parents ii) encode NILAS line genotypes according to parent-of-origin based on the remaining homozygous polymorphic markers iii) aggregate consecutive encoded genotypes into haplotype blocks iv) impute missing marker and haplotype blocks based on adjacent contig congruence v) call second alternate alleles vi)
exports genomic coordinates and encoded/imputed genotype matrices for visualization in FlapJack v1.8.0 (Milne et al. 2010) and vii) convert haplotype block marker data into physical coordinates for characterization.

## Installation:
Using SSH protocol

`$ git clone git@github.com:maizeatlas/nilas_haploimputer.git`

