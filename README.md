# nilas_haploimputer
## Overview of nilas_haploimputer:

NILAS_haploimputer was written in Python 3.6; data frames constructed using Pandas v0.25.0 (McKinney 2010). NILAS_haploimputer selects the donor and recurrent parent genotypes (stored in genotype matrix of all parent lines) and then parses the corresponding NILAS lines for a given group (stored in a genontype matrix of all NILAS lines), concatenating the founder and NILAS lines into one genotype matrix in order to: i) filter heterozygous and non-polymorphic markers from donor/recurrent parents; ii) encode NILAS line genotypes according to parent-of-origin based on the remaining homozygous polymorphic markers; iii) aggregate consecutive encoded genotypes into haplotype blocks; iv) impute missing marker and haplotype blocks based on adjacent read contig (local haplotype) congruence; v) infer second alternate alleles (presumed to originate from residual heterozygosity in the parents); vi)
export genomic coordinates along with encoded and imputed genotype matrices for visualization in FlapJack v1.8.0 (Milne et al. 2010); and vii) convert haplotype block marker data into inferred breakpoint data for whole genome summaries.

### Parent-of-Origin Encoding:
NILAS_haploimputer.py was written to encode each marker as a recurrent, donor, heterozygous, or missing state using a series of conditions. The script merges the genotype data on the parental lines from a consensus genotype matrix with genotype data for the NILAS lines into one genotype matrix. The genotype matrix coordinates are filtered to exclude markers that are heterozygous, missing or non-polymorphic in the parents. The subsetted homozygous-polymorphic parental markers are then used to encode parent-of-origin information for markers scored on the NILAS lines.
### Genotype Imputation:  
A haplotype block is a section of the genome inherited from one of the two parents. Once markers have been encoded by NILAS_haploimputer.py, the script then defines haplotype blocks by iterating through the encoded genotypes. Consecutive markers of the same encoded parent-of-origin are aggregated into a haplotype block. A change in parent-of-origin between two adjacent markers is considered a crossover site. Once haplotype blocks are defined, NILAS_haploimputer.py applies a series of conditions to impute: i) missing data; ii) haplotype blocks that are supported by markers within and across read contigs; iii) haplotype blocks called as secondary alternate alleles. 
### Missing Data Imputation:  
The first condition applied by NILAS_haploimputer.py imputes missing marker data by assigning flanking markers to the 5’ and 3’ of each haplotype block. Haplotype blocks defined as missing are then subject to imputation by comparing the encoded genotypes of the flanking markers. If the markers match then the haplotype block is imputed as the genotype of adjacent haplotype block.  
### Marker Coverage Imputation:  
In the VCF filtering stage of analysis markers are selected by identity to in silico restriction enzyme digest cut sites. Pybedtools is utilized to assign RE site identifiers to NILAS markers at intersecting coordinates. The number of unique in silico sites, informs the number of contigs supporting the haplotype block. This contig number is weighted in imputation; haplotype blocks with less than two supporting contigs are imputed to the nearest neighboring haplotype block genotype. 
### Secondary Alternate Allele Imputation:
Secondary alternative haplotypes from recurrent parents were distinguished from donor parents based on their frequency of occurrence across different lines. Outside of the target interval, donor parent alleles are unlikely to occur by chance at the same site multiple times, as the breeding scheme involved selection of separate individuals for each introgression and included enrichment of the genetic background toward the recurrent parent. Therefore, alleles that did not match the recurrent parent but occured at a frequency of >20% within a given NILAS group (sets included 96 lines belonging to the same donor parent by recurrent parent cross, including the 48 introgression lines the 48 complement lines) were encoded as alternate recurrent parent alleles (RP2).  

Residual heterozygosity was also observed in the donor parent lines, leading to the expectation that multiple donor alleles might have been captured within NILAS. A similar function was employed as described above, but was restricted to evaluation of the ZmPR target interval per NILAS group. This group specificity required an adjustment of allele frequency to account for the change in population, 24 lines (12 introgression and complement lines). A frequency of >95% observed recurrent parent genotypes was set per marker in ZmPR target intervals. The allele frequency condition was contingent on flanking marker logic so that alternate donor genotypes were only called when flanked by haplotype blocks with donor genotypes. A tertiary condition was then applied to the length of the alternate donor haplotype block. The alternate donor allele (DP2) would only be assigned if the distance between adjacent haplotype blocks was less than 10 Mb -- assuming a double recombination event unlikely in this genomic space.  

With the imputation rules applied by NILAS_haploimputer.py, less than 0.1% of the genotype data remains as missing data. The imputed genotypes allow for annotation of introgression lines, providing valuable information on genotype composition, introgression number, recombinant breakpoint resolution, introgression size and position.  
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

      - NILAS Genotype Extracted VCF file - ex. NILAS_g31.GT  
      - NewCon.txt: Founder Genotype Consenus File  
      - isdigB73v4_pid96_convert.bed: In silico digestion bed file  
      - multimode.txt: multimodal marker sites excluded in Founder genotype consensus file  

**All input files must be tab separated**

## Output:

### FlapJack Visualization Files:
1. Marker Coordinate file:[Group]_Coordinates.txt 
2. Encoded genotype matrix for NILAS group: [Group]_Encoded_Genotypes.txt
3. Imputed genotype matrix for NILAS group: [Group]_Imputed_Genotypes.txt
4. ZmPR Genotype matrices per subgroup for introgression and complement lines: [Group]_[ZMPR]_fig.txt/COMP_fig.txt  
### Genotype Composition Summary Files: 
1. Group bed file containing marker-in silico-contig intersection coordinates and identifer:[Group].bed
2. Genotype Composition per subgroup files for introgression and complement lines: [Group]_[ZMPR]_Geno_Sum.txt
   - contains percentages of Foreground/Background Donor/Recurrent/Heterzygous/Missing Data
3. Subgroup haplotype block summary files: [Group]_[ZMPR]_SUM.txt/SUMC.txt
   - contains data per haplotype block including: genotype, identifier, marker number, 5/3' breakpoint density, 5/3' breakpoint coordinate, physical genomic length, ZmPR status, foreground/background-- percentage/ 5/3' breakpoint density
   

## Usage:

`$ python NILAS_haploimputer.py 2369C CML277C g31 '/Users/Todd/Dropbox/NILAS/' NILAS_g31.GT`
