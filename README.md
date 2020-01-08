# nilas_haploimputer
## Overview of nilas_haploimputer:

NILAS_haploimputer was written in Python 3.6; data frames constructed using Pandas v0.25.0 (McKinney 2010). NILAS_haploimputer selects the donor and recurrent parent genotype matrices from consensus matrix and then parses the corresponding NILAS lines from the genotype extracted VCF; concatenating the founder and NILAS progeny into one genotype matrix in order to:  i) filter heterozygous and non-polymorphic markers from donor/recurrent parents ii) encode NILAS line genotypes according to parent-of-origin based on the remaining homozygous polymorphic markers iii) aggregate consecutive encoded genotypes into haplotype blocks iv) impute missing marker and haplotype blocks based on adjacent contig congruence v) call second alternate alleles vi)
exports genomic coordinates and encoded/imputed genotype matrices for visualization in FlapJack v1.8.0 (Milne et al. 2010) and vii) convert haplotype block marker data into physical coordinates for characterization.


## Dependencies:

* pybedtools <https://daler.github.io/pybedtools/>  >= v0.7.10
* numPY <https://numpy.org/> >= v1.17.0
* pandas <https://pandas.pydata.org/> >= v0.25.0 

## Installation:

`$ git clone git@github.com:maizeatlas/nilas_haploimputer.git`

## Input:
1. Recurrent Parent Consensus Genotype - ex. 2369C
2. Donor Parent Consensus Genotype - ex. CML277C
3. Group Number - ex. g31
4. Directory containing:  
    -NILAS Genotype Extracted VCF file - ex. NILAS_g31.GT  
    -NewCon.txt: Founder Genotype Consenus File  
    -isdigB73v4_pid96_convert.bed: In silico digestion bed file  
    -multimode.txt: multimodal marker sites excluded in Founder genotype consensus file  

**All input files must be tab seperated**

## Output:
1. Group bed file containing marker-in silico-contig intersection coordinates and identifer
2. Encoded genotype matrix for NILAS group 
3. Imputed genotype matrix for NILAS group
4. Genotype matrices per subgroup for introgression and complement lines
5. Genotype Composition per subgroup files for introgression and complement lines
6. Subgroup summary files 

## Usage:

`$ python NILAS_haploimputer.py 2369C CML277C g31 '/Users/Todd/Dropbox/NILAS/' NILAS_g31.GT`
