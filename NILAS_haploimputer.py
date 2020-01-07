#!/usr/bin/python

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
# Author:	Todd D. Yoder
# 
# Date:	September 2018
# 
# NILAS_haploimputer was written in Python 3.6; data frames constructed using Pandas v0.19.2 (McKinney 2010). NILAS_haploimputer selects the donor
# and recurrent parent genotype matrices from consensus matrix and then parses the corresponding NILAS lines from the genotype extracted VCF; concatenating
# the founder and NILAS progeny into one genotype matrix in order to:  i) filter heterozygous and non-polymorphic markers from donor/recurrent parents
# ii) encode NILAS line genotypes according to parent-of-origin based on the remaining homozygous polymorphic markers iii) aggregate consecutive encoded
# genotypes into haplotype blocks iv) impute missing marker and haplotype blocks based on adjacent contig congruence v) call second alternate alleles vi)
# exports genomic coordinates and encoded/imputed genotype matrices for visualization in FlapJack v1.8.0 (Milne et al. 2010) and vii) convert haplotype block
# marker data into physical coordinates for characterization.


#Modules:--------------------------------------------------------------------------------------------------------------

import re, sys
import numpy as np
import pandas as pd
import os
import pybedtools



#User Variables:--------------------------------------------------------------------------------------------------------------

Recurrent_Parent = sys.argv[1]
Donor_Parent = sys.argv[2]
NILAS_Crossing_Scheme = sys.argv[3]
dir = sys.argv[4]
NILAS_gt = sys.argv[5]

path = dir + '/' + sys.argv[1] + "x" + sys.argv[2] + '/'

new_dir = path.replace('/','//')

if os.path.isdir(new_dir) != True:
    
    print(new_dir)
    NILASdir = os.makedirs(new_dir)
    print('Generating Directory...',NILASdir)
else:
    
    NILASdir = path
    print("Directory Present:",NILASdir)



# #Global Declarations:--------------------------------------------------------------------------------------------------------------

 
CHROMend = {'1':307041717,'2':244442276,'3':235667834,'4':246994605,'5':223902240,'6':174033170,
            '7':182381542,'8':181122637, '9':159769782, '10':150982314}

B73len =  2106338117 

Scheme_dict = {'g11':'CML10 x 2369','g21':'CML258 x 2369','g31':'CML277 x 2369','g41':'CML341 x 2369','g51':'CML373 x 2369',
               'g61':'Tzi8 x 2369','g71':'Tzi9 x 2369',
               'g12':'CML10 x LH123Ht', 'g22':'CML258 x LH123Ht', 'g32':'CML277 x LH123Ht', 'g42':'CML341 x LH123Ht', 'g52':'CML373 x LH123Ht',
               'g62':'Tzi8 x LH123Ht', 'g72':'Tzi9 x LH123Ht'}

Geno_dict = {'A':'DP','G':'RP', '0':'NA','C':'Het', 0:'NA'}

ZmPR_list = ['Sample','Scheme', 'ZmPR', 'CHROM', 'Hap_Geno', 'Haplotype_Block', 'Hap_No',
             'Marker_No', '5_BP_Density', '3_BP_Density', '5_Breakpoint', '3_Breakpoint', 'PhyLen', 'ZmPR+', 'Foreground_Perc',
             'FG_5_BP_Density', 'FG_3_BP_Density', 'BG_5_BP_Density', 'BG_3_BP_Density', 'Background_Perc']

NILAS_Dict = {}
Coordinates = ['CHROM', 'POS']
#--------------------------------------------------------------------------------------------------------------------

# Import GT Files (manual CoO filter)---------------------------------------------------------------------------------------------
print ("Importing GT file----- YOU HAVE A MANUAL COORDINATE FILTER ON...")
gtdata = pd.read_csv(dir + '/' + NILAS_gt, sep="\t", dtype =object)

NonCon = pd.read_csv(dir + '/' +  'multimode.txt', sep="\t", dtype =object)               #<--------- Manual Coordinate Filter
NonConlist = NonCon['POS']
gtdata = gtdata[~gtdata.POS.isin(NonConlist)]

Consensus = pd.read_csv(dir + '/' + 'NewCon.txt', sep="\t", dtype =object)      #<--------- DP RP consensus genotypes w/ MultiMode Sites Filtered 
Consensus = Consensus[~Consensus.POS.isin(NonConlist)]         

# Select Founder Lines________________________________________________________________________________________________________________
print("Selecting Founder Lines from Consensus File.... ")

Parents = ['CHROM','POS',sys.argv[1],sys.argv[2]]
Consensus = Consensus.filter((Parents))

# Select NILAS Lines________________________________________________________________________________________________________________
print("Selecting NILAS Lines from GT File.... ")

Select = gtdata[gtdata.columns[gtdata.columns.to_series().str.contains(sys.argv[3])]]   #<--------- Selects NILAS lines of breeding group
NILAS  = Select.rename(columns = lambda x : str(x)[:-7])                       #<--------- Removes ".single" from samples
NILASlist = (NILAS.columns[NILAS.columns.to_series().str.contains(sys.argv[3])])
SampleNo = len(NILASlist)
print(SampleNo)

# Concatenate Consensus and NILAS DFs________________________________________________________________________________________________________________
print("Concatenate Consensus and NILAS DFs.....")

gtdataNILAS = pd.concat([Consensus, NILAS], axis=1)

# Filter DP/RP Sites-------------------------------------------
print('Filtering DP/RP Sites: NonVariant/Missing Data....')

# Compile Het Sites Coordinate and Export:
NILAShet = gtdataNILAS[(gtdataNILAS.iloc[:,2] == '0/1') | (gtdataNILAS.iloc[:,3] == '0/1')]
NILAShet = NILAShet.iloc[:,0:4]
NILAShet.to_csv(dir+ '/' + sys.argv[1]+'x'+sys.argv[2]+r'_DP_RP_HET.txt', sep='\t')  #-------------> Parental Heterozygous sites{Export}

# Filter Missing and Het Sites from Consensus Parents: 
gtdataNILAS = gtdataNILAS[(gtdataNILAS.iloc[:,2] != gtdataNILAS.iloc[:,3]) &
    (gtdataNILAS.iloc[:,2] != './.') & (gtdataNILAS.iloc[:,3] != './.') &
    (gtdataNILAS.iloc[:,2] != '0/1') & (gtdataNILAS.iloc[:,3] != '0/1')]

#gtdataNILAS.to_csv(sys.argv[1]+'x'+sys.argv[2]+'_GENO9.txt', sep='\t')       #-------------> ORG Genotype Outfile{Export}
# Slice = gtdataNILAS.drop(gtdataNILAS.columns[6:], axis=1)

# Conditional Block:________________________________________________________________________________________________________________
print("Encoding NILAS.....")
# 
for sample in NILASlist:
    gtdataNILAS.loc[(gtdataNILAS.iloc[:,2] != gtdataNILAS[(sample)]) & ((gtdataNILAS[(sample)]) != gtdataNILAS.iloc[:,3]) &
    (gtdataNILAS[(sample)] != './.') & (gtdataNILAS[(sample)] != '0/1'), [(sample)]] = 'T' #Ambiguous 
    gtdataNILAS.loc[(gtdataNILAS.iloc[:,2].astype(str) == gtdataNILAS[(sample)]) & (gtdataNILAS.iloc[:,3].astype(str) != gtdataNILAS[(sample)]), [(sample)]] = 'G' #Recurrent
    gtdataNILAS.loc[(gtdataNILAS.iloc[:,3].astype(str) == gtdataNILAS[(sample)]) & (gtdataNILAS.iloc[:,2].astype(str) != gtdataNILAS[(sample)]),[(sample)]] = 'A' #Donor
    gtdataNILAS.loc[gtdataNILAS[(sample)] == './.', [(sample)]] = '0' #No Call
    gtdataNILAS.loc[gtdataNILAS[(sample)] == '0/1', [(sample)]] = 'C' #Het

# Slice2 = gtdataNILAS.drop(gtdataNILAS.columns[6:], axis=1)
# Slice2.to_csv(sys.argv[1]+'x'+sys.argv[2]+'_sliceEncode5.1.txt', sep='\t')
# gtdataNILAS.to_csv(os.path.join(path,sys.argv[1]+'x'+sys.argv[2]+r'_Encoded_Geno.txt'), sep='\t') #-------------> Encoded Genotype Outfile{Export}

# Task 6:Impute RP2 Loci:________________________________________________________________________________________________________________
print("Imputing RP2 Loci...")
sample_list = (gtdataNILAS.columns[gtdataNILAS.columns.to_series().str.contains('NILAS')])
Num = len(sample_list)

gtdataNILAS['Count_Donor'] = (gtdataNILAS.iloc[:,1:Num].astype(str) == 'A').sum(axis=1) #Donor
gtdataNILAS['Count_Amb'] = (gtdataNILAS.iloc[:,1:Num].astype(str) == 'T').sum(axis=1) #Ambiguous 
gtdataNILAS['Count_Het'] = (gtdataNILAS.iloc[:,1:Num].astype(str) == 'C').sum(axis=1) #Het

gtdataNILAS['RP2_Loci'] = ((gtdataNILAS[['Count_Donor','Count_Amb','Count_Het']]).sum(axis=1).div(Num))*100

for sample in sample_list:

    gtdataNILAS[(sample)] = np.where((gtdataNILAS['RP2_Loci'] > 20) & (gtdataNILAS[(sample)] != 'G'), '1', gtdataNILAS[(sample)])
    
#gtdataNILAS.to_csv(os.path.join(path,sys.argv[1]+'x'+sys.argv[2]+r'_RP2.txt'), sep='\t') #-------------> Residual Heterozygous Outfile{Export}

# ZmPR Definition:________________________________________________________________________________________________________________
ZmPR1_QTL = gtdataNILAS[(gtdataNILAS['CHROM'] ==  '1') &  (gtdataNILAS['POS'].astype(int).between(57191440, 185363076, inclusive=True))]
ZmPR2_QTL = gtdataNILAS[(gtdataNILAS['CHROM'] ==  '8') &  (gtdataNILAS['POS'].astype(int).between(97682806, 167995313, inclusive=True))]
ZmPR3_QTL = gtdataNILAS[(gtdataNILAS['CHROM'] ==  '9') &  (gtdataNILAS['POS'].astype(int).between(7384635, 119857067, inclusive=True))]
ZmPR4_QTL = gtdataNILAS[(gtdataNILAS['CHROM'] ==  '10') &  (gtdataNILAS['POS'].astype(int).between(10126632, 136329803, inclusive=True))]

ZmPR1_POS = ZmPR1_QTL['POS'].tolist()
ZmPR2_POS = ZmPR2_QTL['POS'].tolist()
ZmPR3_POS = ZmPR3_QTL['POS'].tolist()
ZmPR4_POS = ZmPR4_QTL['POS'].tolist()

ZmPR_POS_list = ZmPR1_POS + ZmPR2_POS + ZmPR3_POS + ZmPR4_POS

# Intersect Group BED file with ISC BED:________________________________________________________________________________________________________________

CoBed = gtdataNILAS.filter((Coordinates)).reset_index(drop=True)
CoBed = CoBed[~CoBed['CHROM'].str.contains('B73V4_ctg')]

CoBed['stop'] = CoBed.iloc[:,1].astype(int) + 1

CoBed.to_csv(path +sys.argv[1]+'x'+sys.argv[2]+r'.bed', sep='\t', header= False, index= False)

GroupBed = pybedtools.BedTool(path + '/'+sys.argv[1]+'x'+sys.argv[2]+r'.bed')
B73isc = pybedtools.BedTool(dir + '/' + 'isdigB73v4_pid96_convert.bed')

U_isc = pybedtools.BedTool.intersect(GroupBed,B73isc, wb=True)

U_isc2 = U_isc.saveas(path +'//' +sys.argv[1]+'x'+sys.argv[2]+'.bed')

isc = pd.read_csv(path+sys.argv[1]+'x'+sys.argv[2]+'.bed', sep="\t", dtype =object,names = ["CHR", "Start", "Stop", "Coverage","Contig_Start","Contig_Stop","isc"])

# Imputation & Annotation:________________________________________________________________________________________________________________
print("Imputing NILAS Haplotypes.....")

gtdataNILAS = gtdataNILAS[~gtdataNILAS['CHROM'].str.contains('B73V4_ctg')]

for sample in NILASlist:
   
    gtdataNILAS[(sample)+'_Haplotype_Block'] = ((gtdataNILAS.CHROM != gtdataNILAS.CHROM.shift(1)) |
    (gtdataNILAS[(sample)] != ((gtdataNILAS[(sample)].astype(str)).shift(1)))).astype(int).cumsum()
    
    gtdataNILAS[(sample)+'_Start'] = gtdataNILAS.groupby([(sample)+'_Haplotype_Block'])['POS'].transform('first')
    gtdataNILAS[(sample)+'_Stop'] = gtdataNILAS.groupby([(sample)+'_Haplotype_Block'])['POS'].transform('last')
    
    gtdataNILAS[(sample)+'_Five_Flank'] = gtdataNILAS.groupby([(sample)+'_Stop',(sample)+'_Haplotype_Block'])['POS'].transform('last').shift(1)
    gtdataNILAS[(sample)+'_Three_Flank'] = gtdataNILAS.groupby([(sample)+'_Start',(sample)+'_Haplotype_Block'])['POS'].transform('first').shift(-1)
    
    gtdataNILAS.loc[(gtdataNILAS[(sample)+'_Start'].astype(str) == gtdataNILAS[(sample)+'_Three_Flank'].astype(str)), [(sample)+'_Three_Flank']] = 'NaN'
    gtdataNILAS[(sample)+'_Three_Flank'] = gtdataNILAS[(sample)+'_Three_Flank'].replace(to_replace=['NaN'], regex=True, method='bfill')
    
    gtdataNILAS.loc[(gtdataNILAS[(sample)+'_Stop'].astype(str) == gtdataNILAS[(sample)+'_Five_Flank'].astype(str)), [(sample)+'_Five_Flank']] = 'NaN'
    gtdataNILAS.loc[(gtdataNILAS[(sample)+'_Five_Flank'].astype(str) == 'NaN'), [(sample)+'_Three_Flank']] = 'NaN'
    
    gtdataNILAS = gtdataNILAS.fillna(0)
    gtdataNILAS = gtdataNILAS.replace(to_replace=['NaN'], value=0)
    # 
    gtdataNILAS[(sample)+'_Three_Flank'] = np.where(gtdataNILAS[(sample)+'_Three_Flank'].astype(int) < gtdataNILAS[(sample)+'_Stop'].astype(int),
                                                     (gtdataNILAS['CHROM'].map(CHROMend).astype(int)),gtdataNILAS[(sample)+'_Three_Flank'])   #Designates chromosome ends 
    
    
    #__________________________ Zero Distance Scenerio__________________________ 
    gtdataNILAS[(sample)+'_Five_Distance'] = gtdataNILAS[(sample)+'_Stop'].astype(int) - gtdataNILAS[(sample)+'_Five_Flank'].astype(int)
    gtdataNILAS[(sample)+'_Three_Distance'] = gtdataNILAS[(sample)+'_Three_Flank'].astype(int) - gtdataNILAS[(sample)+'_Stop'].astype(int)

    gtdataNILAS.loc[(gtdataNILAS[(sample)+'_Five_Flank'] == 0), [(sample)+'_Five_Distance']] = 0
    gtdataNILAS.loc[(gtdataNILAS[(sample)+'_Three_Flank'] == 0), [(sample)+'_Three_Distance']] = 0
    
    gtdataNILAS[(sample)+'_Five_Distance'] = gtdataNILAS[(sample)+'_Five_Distance'].replace(to_replace= 0, regex=True, method='ffill')
    gtdataNILAS[(sample)+'_Three_Distance'] = gtdataNILAS[(sample)+'_Three_Distance'].replace(to_replace= 0, regex=True, method='ffill')
     
    gtdataNILAS[(sample)+'_Five_Flank_Encode'] = gtdataNILAS.groupby([(sample)+'_Haplotype_Block'])[(sample)].transform('last').shift(1)
    gtdataNILAS[(sample)+'_Five_Flank_Encode'] = gtdataNILAS[(sample)+'_Five_Flank_Encode'].replace(to_replace= '0', regex=True, method='ffill')
    
    gtdataNILAS[(sample)+'_Three_Flank_Encode'] = gtdataNILAS.groupby([(sample)+'_Haplotype_Block'])[(sample)].transform('first').shift(-1)
    gtdataNILAS[(sample)+'_Three_Flank_Encode'] = gtdataNILAS[(sample)+'_Three_Flank_Encode'].replace(to_replace= '0', regex=True, method='bfill')
    
    gtdataNILAS[(sample)+'_Imputed'] = gtdataNILAS[(sample)]
    
    #__________________________Impute Missing Data by Flanking Marker Homology__________________________ 
    
    gtdataNILAS[(sample)+'_Imputed'] = np.where((gtdataNILAS[(sample)+'_Imputed'] == '0') &
     (gtdataNILAS[(sample)+'_Five_Flank_Encode'] == gtdataNILAS[(sample)+'_Three_Flank_Encode']),
     gtdataNILAS[(sample)+'_Imputed'].replace(to_replace= '0', regex=True, method='ffill'), gtdataNILAS[(sample)+'_Imputed'])
    
    gtdataNILAS[(sample)+'_Haplotype_Block'] = ((gtdataNILAS.CHROM != gtdataNILAS.CHROM.shift(1)) |
        (gtdataNILAS[(sample)+'_Imputed'] != ((gtdataNILAS[(sample)+'_Imputed'].astype(str)).shift(1)))).astype(int).cumsum()
    
    gtdataNILAS[(sample)+'_Hap_Marker_No'] = gtdataNILAS.groupby(['CHROM', (sample)+'_Haplotype_Block'])['CHROM'].transform('size')    
 
#__________________________ Impute by Contig Congruence __________________________ 

isc = isc.reset_index(drop=True)
gtdataNILAS = gtdataNILAS.reset_index(drop=True)

gtdataNILAS = pd.concat([isc['isc'],gtdataNILAS], axis=1)

for sample in NILASlist:
    
    gtdataNILAS[(sample)+'_isc_no'] = gtdataNILAS.groupby((sample)+'_Haplotype_Block')['isc'].transform('nunique')
    
    gtdataNILAS[(sample)+'_isc_no'] = np.where((gtdataNILAS[(sample)+'_isc_no'] < 2) &
     (gtdataNILAS[(sample)+'_Imputed'] == 'A') & ((gtdataNILAS['POS']).isin(ZmPR_POS_list)), 5, gtdataNILAS[(sample)+'_isc_no'])  #<-- Increases sensitivity for Donor gentoypes in ZmPRs
        
    gtdataNILAS[(sample)+'_isc_no'] = np.where((gtdataNILAS[(sample)+'_isc_no'] <= 3) &
     (gtdataNILAS[(sample)+'_Imputed'] != 'A') & (gtdataNILAS[(sample)+'_Imputed'] != 'A2') &
     ((gtdataNILAS['POS']).isin(ZmPR_POS_list) &
        ((gtdataNILAS[(sample)+'_Stop'].astype(int)) - (gtdataNILAS[(sample)+'_Start'].astype(int))) < 10000000), 1, gtdataNILAS[(sample)+'_isc_no']) #<-- imputation of non-donor markers for near distance contigs
    
    gtdataNILAS[(sample)+'_Imputed'] = np.where((gtdataNILAS[(sample)+'_isc_no'] < 2), 'NaN', gtdataNILAS[(sample)+'_Imputed'])     #Where haplotype block is <2 contigs, replace with NaN, otherwise use orginal genotype
     
    gtdataNILAS.loc[((gtdataNILAS[(sample)+'_isc_no'] < 2) &
        (gtdataNILAS[(sample)+'_Five_Distance'] < gtdataNILAS[(sample)+'_Three_Distance'])),
        [(sample)+'_Imputed']] = gtdataNILAS[(sample)+'_Imputed'].replace(to_replace='NaN', regex=True, method='ffill')
    
    gtdataNILAS.loc[((gtdataNILAS[(sample)+'_isc_no'] < 2) &
        (gtdataNILAS[(sample)+'_Five_Distance'] > gtdataNILAS[(sample)+'_Three_Distance'])),
        [(sample)+'_Imputed']] = gtdataNILAS[(sample)+'_Imputed'].replace(to_replace='NaN', regex=True, method='bfill')
    
    gtdataNILAS.loc[((gtdataNILAS[(sample)+'_isc_no'] < 2) &
         (gtdataNILAS[(sample)+'_Five_Distance'].abs() == 0)), [(sample)+'_Imputed']] = gtdataNILAS[(sample)+'_Imputed'].replace(to_replace='NaN', regex=True, method='bfill')
    
    gtdataNILAS[(sample)+'_Imp_HapBlock'] = ((gtdataNILAS.CHROM != gtdataNILAS.CHROM.shift(1)) |
    (gtdataNILAS[(sample)+'_Imputed'] != ((gtdataNILAS[(sample)+'_Imputed'].astype(str)).shift(1)))).astype(int).cumsum()
       

Slice3 = gtdataNILAS[gtdataNILAS.columns[gtdataNILAS.columns.to_series().str.contains('NILASq2g41i08s2p1r1')]]
Coordinates2 = gtdataNILAS[['CHROM', 'POS', 'isc']]
Slice3 = pd.concat([Coordinates2,Slice3], axis=1)
# Slice3.to_csv(sys.argv[1]+'x'+sys.argv[2]+'_sliceINQ_X.txt', sep='\t')    #-------------> Single Genotype  Outfile{Export}
###_________________________________________________________________________________________________________________________________________________________
#Task 8: Calculating 2nd Alt Genotype Frequency...
print("Calculating 2nd Alt Genotype Frequency...")

###ZmPR1_________________________________________________________________________________________________
ZmPR1 = gtdataNILAS[(gtdataNILAS['CHROM'] ==  '1') &  (gtdataNILAS['POS'].astype(int).between(57191440, 185363076, inclusive=True))]
ZmPR1CoO = ZmPR1.filter((Coordinates))
Compliment = (ZmPR1.columns[~ZmPR1.columns.to_series().str.contains('C')])
ZmPR1_C = (ZmPR1.columns[ZmPR1.columns.to_series().str.contains('C')])
ZmPR1_COMP = ZmPR1.filter(ZmPR1_C)
ZmPR1 = ZmPR1.filter(Compliment)
ZmPR1 = ZmPR1[ZmPR1.columns[ZmPR1.columns.to_series().str.contains('1g')]]
ZmPR1 = ZmPR1[ZmPR1.columns[ZmPR1.columns.to_series().str.contains('Imputed')]]
ZmPR1 = pd.concat([ZmPR1CoO, ZmPR1], axis=1)

Zm1sample_list = (ZmPR1.columns[ZmPR1.columns.to_series().str.contains('1g')])

ZmPR1_COMP = ZmPR1_COMP[ZmPR1_COMP.columns[ZmPR1_COMP.columns.to_series().str.contains('1g')]]
ZmPR1_COMP = ZmPR1_COMP[ZmPR1_COMP.columns[ZmPR1_COMP.columns.to_series().str.contains('Imputed')]]
ZmPR1_COMP = pd.concat([ZmPR1CoO, ZmPR1_COMP], axis=1)

Zm1sample_complist = (ZmPR1_COMP.columns[ZmPR1_COMP.columns.to_series().str.contains('1g')])



Zm1Num = len(Zm1sample_list)
ZmPR1['Count_RP'] = (ZmPR1.iloc[:,2:(Zm1Num +2)].astype(str) == 'G').sum(axis=1) #Recurrent
ZmPR1['Count_Amb'] = (ZmPR1.iloc[:,2:(Zm1Num +2)].astype(str) == 'T').sum(axis=1) #Ambiguous 
ZmPR1['Count_Het'] = (ZmPR1.iloc[:,2:(Zm1Num +2)].astype(str) == 'C').sum(axis=1) #Het

ZmPR1['DP2_Loci_1'] = ((ZmPR1[['Count_RP','Count_Amb','Count_Het']]).sum(axis=1).div(Zm1Num))*100

###ZmPR2_________________________________________________________________________________________________
ZmPR2 = gtdataNILAS[(gtdataNILAS['CHROM'] ==  '8') &  (gtdataNILAS['POS'].astype(int).between(97682806, 167995313, inclusive=True))]
ZmPR2CoO = ZmPR2.filter((Coordinates))
Compliment = (ZmPR2.columns[~ZmPR2.columns.to_series().str.contains('C')])
ZmPR2_C = (ZmPR2.columns[ZmPR2.columns.to_series().str.contains('C')])
ZmPR2_Comp = ZmPR2.filter(ZmPR2_C)
ZmPR2 = ZmPR2.filter(Compliment)
ZmPR2 = ZmPR2[ZmPR2.columns[ZmPR2.columns.to_series().str.contains('2g')]]
ZmPR2 = ZmPR2[ZmPR2.columns[ZmPR2.columns.to_series().str.contains('Imputed')]]
ZmPR2 = pd.concat([ZmPR2CoO, ZmPR2], axis=1)

Zm2sample_list = (ZmPR2.columns[ZmPR2.columns.to_series().str.contains('2g')])

ZmPR2_Comp = ZmPR2_Comp[ZmPR2_Comp.columns[ZmPR2_Comp.columns.to_series().str.contains('2g')]]
ZmPR2_Comp = ZmPR2_Comp[ZmPR2_Comp.columns[ZmPR2_Comp.columns.to_series().str.contains('Imputed')]]
ZmPR2_Comp = pd.concat([ZmPR2CoO, ZmPR2_Comp], axis=1)

Zm2sample_complist = (ZmPR2_Comp.columns[ZmPR2_Comp.columns.to_series().str.contains('2g')])
print(Zm2sample_complist)

Zm2Num = len(Zm2sample_list)
ZmPR2['Count_RP'] = (ZmPR2.iloc[:,2:(Zm2Num +2)].astype(str) == 'G').sum(axis=1) #Recurrent
ZmPR2['Count_Amb'] = (ZmPR2.iloc[:,2:(Zm2Num +2)].astype(str) == 'T').sum(axis=1) #Ambiguous 
ZmPR2['Count_Het'] = (ZmPR2.iloc[:,2:(Zm2Num +2)].astype(str) == 'C').sum(axis=1) #Het

ZmPR2['DP2_Loci_2'] = ((ZmPR2[['Count_RP','Count_Amb','Count_Het']]).sum(axis=1).div(Zm2Num))*100

###ZmPR3_________________________________________________________________________________________________
ZmPR3 = gtdataNILAS[(gtdataNILAS['CHROM'] ==  '9') &  (gtdataNILAS['POS'].astype(int).between(7384635, 119857067, inclusive=True))]
ZmPR3CoO = ZmPR3.filter((Coordinates))
Compliment = (ZmPR3.columns[~ZmPR3.columns.to_series().str.contains('C')])
ZmPR3_C = (ZmPR3.columns[ZmPR3.columns.to_series().str.contains('C')])
ZmPR3_Comp = ZmPR3.filter(ZmPR3_C)
ZmPR3 = ZmPR3.filter(Compliment)
ZmPR3 = ZmPR3[ZmPR3.columns[ZmPR3.columns.to_series().str.contains('3g')]]
ZmPR3 = ZmPR3[ZmPR3.columns[ZmPR3.columns.to_series().str.contains('Imputed')]]
ZmPR3 = pd.concat([ZmPR3CoO, ZmPR3], axis=1)

Zm3sample_list = (ZmPR3.columns[ZmPR3.columns.to_series().str.contains('3g')])


ZmPR3_Comp = ZmPR3_Comp[ZmPR3_Comp.columns[ZmPR3_Comp.columns.to_series().str.contains('3g')]]
ZmPR3_Comp = ZmPR3_Comp[ZmPR3_Comp.columns[ZmPR3_Comp.columns.to_series().str.contains('Imputed')]]
ZmPR3_Comp = pd.concat([ZmPR3CoO, ZmPR3_Comp], axis=1)

Zm3sample_complist = (ZmPR3_Comp.columns[ZmPR3_Comp.columns.to_series().str.contains('3g')])
print(Zm3sample_complist)


Zm3Num = len(Zm3sample_list)
ZmPR3['Count_RP'] = (ZmPR3.iloc[:,2:(Zm3Num +2)].astype(str) == 'G').sum(axis=1) #Recurrent
ZmPR3['Count_Amb'] = (ZmPR3.iloc[:,2:(Zm3Num +2)].astype(str) == 'T').sum(axis=1) #Ambiguous 
ZmPR3['Count_Het'] = (ZmPR3.iloc[:,2:(Zm3Num +2)].astype(str) == 'C').sum(axis=1) #Het

ZmPR3['DP2_Loci_3'] = ((ZmPR3[['Count_RP','Count_Amb','Count_Het']]).sum(axis=1).div(Zm3Num))*100

###ZmPR4_________________________________________________________________________________________________
ZmPR4 = gtdataNILAS[(gtdataNILAS['CHROM'] ==  '10') &  (gtdataNILAS['POS'].astype(int).between(10126632, 136329803, inclusive=True))]
ZmPR4CoO = ZmPR4.filter((Coordinates))
Compliment = (ZmPR4.columns[~ZmPR4.columns.to_series().str.contains('C')])
ZmPR4_C = (ZmPR4.columns[ZmPR4.columns.to_series().str.contains('C')])
ZmPR4_Comp = ZmPR4.filter(ZmPR4_C)
ZmPR4 = ZmPR4.filter(Compliment)
ZmPR4 = ZmPR4[ZmPR4.columns[ZmPR4.columns.to_series().str.contains('4g')]]
ZmPR4 = ZmPR4[ZmPR4.columns[ZmPR4.columns.to_series().str.contains('Imputed')]]
ZmPR4 = pd.concat([ZmPR4CoO, ZmPR4], axis=1)

Zm4sample_list = (ZmPR4.columns[ZmPR4.columns.to_series().str.contains('4g')])

ZmPR4_Comp = ZmPR4_Comp[ZmPR4_Comp.columns[ZmPR4_Comp.columns.to_series().str.contains('4g')]]
ZmPR4_Comp = ZmPR4_Comp[ZmPR4_Comp.columns[ZmPR4_Comp.columns.to_series().str.contains('Imputed')]]
ZmPR4_Comp = pd.concat([ZmPR4CoO, ZmPR4_Comp], axis=1)

Zm4sample_complist = (ZmPR4_Comp.columns[ZmPR4_Comp.columns.to_series().str.contains('4g')])
print(Zm4sample_complist)

Zm4Num = len(Zm4sample_list)
ZmPR4['Count_RP'] = (ZmPR4.iloc[:,2:(Zm4Num +2)].astype(str) == 'G').sum(axis=1) #Recurrent
ZmPR4['Count_Amb'] = (ZmPR4.iloc[:,2:(Zm4Num +2)].astype(str) == 'T').sum(axis=1) #Ambiguous 
ZmPR4['Count_Het'] = (ZmPR4.iloc[:,2:(Zm4Num +2)].astype(str) == 'C').sum(axis=1) #Het

ZmPR4['DP2_Loci_4'] = ((ZmPR4[['Count_RP','Count_Amb','Count_Het']]).sum(axis=1).div(Zm4Num))*100

sample_list = (gtdataNILAS.columns[gtdataNILAS.columns.to_series().str.contains('NILAS')])
Num = len(sample_list)

#Second Alternate Allele BLOCK________________________________________________________________________________________
CoO = gtdataNILAS.filter((Coordinates))
CoO = CoO[~CoO['CHROM'].str.contains('B73V4_ctg')]
NILASimpute = gtdataNILAS[gtdataNILAS.columns[gtdataNILAS.columns.to_series().str.contains('Imputed')]]
HapBlock = gtdataNILAS[gtdataNILAS.columns[gtdataNILAS.columns.to_series().str.contains('_Haplotype_Block')]]

NILAShap = pd.concat([CoO, NILASimpute, ZmPR1['DP2_Loci_1'],ZmPR2['DP2_Loci_2'],ZmPR3['DP2_Loci_3'],ZmPR4['DP2_Loci_4']], axis=1) #<--------RP2 added
NILAShap = NILAShap[NILAShap.CHROM.notnull()]
NILAShap.columns = NILAShap.columns.str.replace('_Imputed','')

print("Imputing 2nd Alt Genotype...")

for sample in NILASlist:
    NILAShap[(sample)+'_Haplotype_Block'] = ((NILAShap.CHROM != NILAShap.CHROM.shift(1)) |
        (NILAShap[(sample)] != ((NILAShap[(sample)].astype(str)).shift(1)))).astype(int).cumsum()
    NILAShap[(sample)+'_Start'] = NILAShap.groupby([(sample)+'_Haplotype_Block'])['POS'].transform('first')
    NILAShap[(sample)+'_Stop'] = NILAShap.groupby([(sample)+'_Haplotype_Block'])['POS'].transform('last')
        
    NILAShap[(sample)+'_Five_Flank'] = NILAShap.groupby([(sample)+'_Stop',(sample)+'_Haplotype_Block'])['POS'].transform('last').shift(1)
    NILAShap[(sample)+'_Three_Flank'] = NILAShap.groupby([(sample)+'_Start',(sample)+'_Haplotype_Block'])['POS'].transform('first').shift(-1)
    
    NILAShap.loc[(NILAShap[(sample)+'_Start'].astype(str) == NILAShap[(sample)+'_Three_Flank'].astype(str)), [(sample)+'_Three_Flank']] = 'NaN'
    NILAShap[(sample)+'_Three_Flank'] = NILAShap[(sample)+'_Three_Flank'].replace(to_replace=['NaN'], regex=True, method='bfill')
    
    NILAShap.loc[(NILAShap[(sample)+'_Stop'].astype(str) == NILAShap[(sample)+'_Five_Flank'].astype(str)), [(sample)+'_Five_Flank']] = 'NaN'
    NILAShap.loc[(NILAShap[(sample)+'_Five_Flank'].astype(str) == 'NaN'), [(sample)+'_Three_Flank']] = 'NaN'
    
    NILAShap = NILAShap.fillna(0)
    NILAShap = NILAShap.replace(to_replace=['NaN'], value=0)

    NILAShap[(sample)+'_Five_Flank_Encode'] = NILAShap.groupby([(sample)+'_Haplotype_Block'])[(sample)].transform('last').shift(1)
    NILAShap[(sample)+'_Five_Flank_Encode'] = NILAShap[(sample)+'_Five_Flank_Encode'].replace(to_replace= '0', regex=True, method='ffill')
    
    NILAShap[(sample)+'_Three_Flank_Encode'] = NILAShap.groupby([(sample)+'_Haplotype_Block'])[(sample)].transform('first').shift(-1)
    NILAShap[(sample)+'_Three_Flank_Encode'] = NILAShap[(sample)+'_Three_Flank_Encode'].replace(to_replace= '0', regex=True, method='bfill')
        
    NILAShap[(sample)+'_Hap_5_Flank'] = NILAShap.groupby(['CHROM', (sample)+'_Haplotype_Block'])[(sample)+'_Five_Flank_Encode'].transform('first') 
    NILAShap[(sample)+'_Hap_3_Flank'] = NILAShap.groupby(['CHROM', (sample)+'_Haplotype_Block'])[(sample)+'_Three_Flank_Encode'].transform('last')
    
    NILAShap[(sample)] = np.where((NILAShap['DP2_Loci_1'] > 95) & (NILAShap[(sample)] != 'C') &
         (NILAShap[(sample)+'_Hap_5_Flank'] == 'A') & (NILAShap[(sample)+'_Hap_3_Flank'] == 'A') &
         (NILAShap.groupby([(sample)+'_Haplotype_Block'])['DP2_Loci_1'].transform('first') ==
          (NILAShap.groupby([(sample)+'_Haplotype_Block'])['DP2_Loci_1'].transform('last'))) &
         (((NILAShap[(sample)+'_Stop'].astype(int)) - (NILAShap[(sample)+'_Start'].astype(int))) < 100000000 ), 'A2',  NILAShap[(sample)])

    NILAShap[(sample)] = np.where((NILAShap['DP2_Loci_2'] > 95) & (NILAShap[(sample)] != 'C') &
         (NILAShap[(sample)+'_Hap_5_Flank'] == 'A') & (NILAShap[(sample)+'_Hap_3_Flank'] == 'A') &
         (NILAShap.groupby([(sample)+'_Haplotype_Block'])['DP2_Loci_2'].transform('first') ==
          (NILAShap.groupby([(sample)+'_Haplotype_Block'])['DP2_Loci_2'].transform('last'))) &
         (((NILAShap[(sample)+'_Stop'].astype(int)) - (NILAShap[(sample)+'_Start'].astype(int))) < 100000000 ), 'A2',  NILAShap[(sample)])
        
    NILAShap[(sample)] = np.where((NILAShap['DP2_Loci_3'] > 95) & (NILAShap[(sample)] != 'C') &
         (NILAShap[(sample)+'_Hap_5_Flank'] == 'A') & (NILAShap[(sample)+'_Hap_3_Flank'] == 'A') &
         (NILAShap.groupby([(sample)+'_Haplotype_Block'])['DP2_Loci_3'].transform('first') ==
          (NILAShap.groupby([(sample)+'_Haplotype_Block'])['DP2_Loci_3'].transform('last'))) &
          (((NILAShap[(sample)+'_Stop'].astype(int)) - (NILAShap[(sample)+'_Start'].astype(int))) < 100000000 ), 'A2',  NILAShap[(sample)])

    NILAShap[(sample)] = np.where((NILAShap['DP2_Loci_4'] > 95) & (NILAShap[(sample)] != 'C') &
         (NILAShap[(sample)+'_Hap_5_Flank'] == 'A') & (NILAShap[(sample)+'_Hap_3_Flank'] == 'A') &
         (NILAShap.groupby([(sample)+'_Haplotype_Block'])['DP2_Loci_4'].transform('first') ==
          (NILAShap.groupby([(sample)+'_Haplotype_Block'])['DP2_Loci_4'].transform('last'))) &
         (((NILAShap[(sample)+'_Stop'].astype(int)) - (NILAShap[(sample)+'_Start'].astype(int))) < 100000000 ), 'A2',  NILAShap[(sample)])
     
    NILAShap[(sample)] = np.where((NILAShap[(sample)] == '1') &
         (NILAShap[(sample)+'_Hap_5_Flank'] == 'A') & (NILAShap[(sample)+'_Hap_3_Flank'] == 'A') &
         (NILAShap.groupby([(sample)+'_Haplotype_Block'])[(sample)+'_Five_Flank_Encode'].transform('first') ==
          (NILAShap.groupby([(sample)+'_Haplotype_Block'])[(sample)+'_Three_Flank_Encode'].transform('last'))), 'A',  NILAShap[(sample)]) ####Condition for re-imputing Donor allele (A) that were called G2
    
    NILAShap[(sample)] = np.where((NILAShap[(sample)] == 'C') &
         (NILAShap[(sample)+'_Hap_5_Flank'] == 'A') & (NILAShap[(sample)+'_Hap_3_Flank'] == 'A') &
         (NILAShap.groupby([(sample)+'_Haplotype_Block'])[(sample)+'_Five_Flank_Encode'].transform('first') ==
          (NILAShap.groupby([(sample)+'_Haplotype_Block'])[(sample)+'_Three_Flank_Encode'].transform('last'))) &
         (((NILAShap[(sample)+'_Stop'].astype(int)) - (NILAShap[(sample)+'_Start'].astype(int))) < 2000000 ), 'A',  NILAShap[(sample)])
    
##Output Options--------------------------------------------------------------------------------------------
#FULL DATA:--------------------------------------------------------------------------------------------
#gtdataNILAS.to_csv(os.path.join(path,sys.argv[1]+'x'+sys.argv[2]+r'_fulldata9.txt'), sep='\t', header= True, index= True)       #------------->Full Outfile {Export}

CoO_Hap = NILAShap.filter((Coordinates))

Zm1sample_list = Zm1sample_list.tolist()
Zm1sample_list = [x[:-8] for x in Zm1sample_list]
Zm2sample_list = Zm2sample_list.tolist()
Zm2sample_list = [x[:-8] for x in Zm2sample_list]
Zm3sample_list = Zm3sample_list.tolist()
Zm3sample_list = [x[:-8] for x in Zm3sample_list]
Zm4sample_list = Zm4sample_list.tolist()
Zm4sample_list = [x[:-8] for x in Zm4sample_list]

#______ FIG BLOCK 1 ____________________#

NILASZmPR1int  = NILAShap[Zm1sample_list]
NILASintro1 = pd.concat([CoO_Hap, NILASZmPR1int], axis=1)
NILAS1Fig= NILASintro1[(NILASintro1['CHROM'] ==  '1') &  (NILASintro1['POS'].astype(int).between(57191440, 185363076, inclusive=True))]
NILAS1Fig= NILAS1Fig.drop(NILAS1Fig.columns[0:2], axis=1)
NILAS1Fig= NILAS1Fig.transpose()
NILAS1Fig.to_csv(path +'/' + sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR1_Fig.txt', sep='\t')

NILASZmPR2int  = NILAShap[Zm2sample_list]
NILASintro2 = pd.concat([CoO_Hap, NILASZmPR2int], axis=1)
NILAS2Fig= NILASintro2[(NILASintro2['CHROM'] ==  '8') &  (NILASintro2['POS'].astype(int).between(97682806, 167995313, inclusive=True))]
NILAS2Fig= NILAS2Fig.drop(NILAS2Fig.columns[0:2], axis=1)
NILAS2Fig= NILAS2Fig.transpose()
NILAS2Fig.to_csv(path +'/' + sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR2_Fig.txt', sep='\t')

NILASZmPR3int  = NILAShap[Zm3sample_list]
NILASintro3 = pd.concat([CoO_Hap, NILASZmPR3int], axis=1)
NILAS3Fig= NILASintro3[(NILASintro3['CHROM'] ==  '9') &  (NILASintro3['POS'].astype(int).between(7384635, 119857067, inclusive=True))]
NILAS3Fig= NILAS3Fig.drop(NILAS3Fig.columns[0:2], axis=1)
NILAS3Fig= NILAS3Fig.transpose()
NILAS3Fig.to_csv(path +'/' + sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR3_Fig.txt', sep='\t')

NILASZmPR4int  = NILAShap[Zm4sample_list]
NILASintro4 = pd.concat([CoO_Hap, NILASZmPR4int], axis=1)
NILAS4Fig= NILASintro4[(NILASintro4['CHROM'] ==  '10') &  (NILASintro4['POS'].astype(int).between(10126632, 136329803, inclusive=True))]
NILAS4Fig= NILAS4Fig.drop(NILAS4Fig.columns[0:2], axis=1)
NILAS4Fig= NILAS4Fig.transpose()
NILAS4Fig.to_csv(path +'/' + sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR4_Fig.txt', sep='\t')

#___________________________________________________________________________________________________________________#


#_________________________FIG BLOCK 2____________________________________________________________________________#
Zm1sample_complist = Zm1sample_complist.tolist()
Zm1sample_complist = [x[:-8] for x in Zm1sample_complist]
Zm2sample_complist = Zm2sample_complist.tolist()
Zm2sample_complist = [x[:-8] for x in Zm2sample_complist]
Zm3sample_complist = Zm3sample_complist.tolist()
Zm3sample_complist = [x[:-8] for x in Zm3sample_complist]
Zm4sample_complist = Zm4sample_complist.tolist()
Zm4sample_complist = [x[:-8] for x in Zm4sample_complist]


NILASZmPR1comp  = NILAShap[Zm1sample_complist]
NILAS_comp1 = pd.concat([CoO_Hap, NILASZmPR1comp], axis=1)
NILAS1figC= NILAS_comp1[(NILAS_comp1['CHROM'] ==  '1') &  (NILAS_comp1['POS'].astype(int).between(57191440, 185363076, inclusive=True))]
NILAS1figC= NILAS1figC.drop(NILAS1figC.columns[0:2], axis=1)
NILAS1figC= NILAS1figC.transpose()
NILAS1figC.to_csv(path +'/' + sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR1_COMP_fig.txt', sep='\t')

NILASZmPR2comp  = NILAShap[Zm2sample_complist]
NILAS_comp2 = pd.concat([CoO_Hap, NILASZmPR2comp], axis=1)
NILAS2figC= NILAS_comp2[(NILAS_comp2['CHROM'] ==  '8') &  (NILAS_comp2['POS'].astype(int).between(97682806, 167995313, inclusive=True))]
NILAS2figC= NILAS2figC.drop(NILAS2figC.columns[0:2], axis=1)
NILAS2figC= NILAS2figC.transpose()
NILAS2figC.to_csv(path +'/' + sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR2_COMP_fig.txt', sep='\t')

NILASZmPR3comp  = NILAShap[Zm3sample_complist]
NILAS_comp3 = pd.concat([CoO_Hap, NILASZmPR3comp], axis=1)
NILAS3figC= NILAS_comp3[(NILAS_comp3['CHROM'] ==  '9') &  (NILAS_comp3['POS'].astype(int).between(7384635, 119857067, inclusive=True))]
NILAS3figC= NILAS3figC.drop(NILAS3figC.columns[0:2], axis=1)
NILAS3figC= NILAS3figC.transpose()
NILAS3figC.to_csv(path +'/' + sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR3_COMP_fig.txt', sep='\t')

NILASZmPR4comp  = NILAShap[Zm4sample_complist]
NILAS_comp4 = pd.concat([CoO_Hap, NILASZmPR4comp], axis=1)
NILAS4figC= NILAS_comp4[(NILAS_comp4['CHROM'] ==  '10') &  (NILAS_comp4['POS'].astype(int).between(10126632, 136329803, inclusive=True))]
NILAS4figC= NILAS4figC.drop(NILAS4figC.columns[0:2], axis=1)
NILAS4figC= NILAS4figC.transpose()
NILAS4figC.to_csv(path +'/' + sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR4_COMP_fig.txt', sep='\t')

##______________________________________________________________________________________________________________##

NILASencode = gtdataNILAS[gtdataNILAS.columns[~gtdataNILAS.columns.to_series().str.contains('_')]]   
NILASencode = NILASencode[~NILASencode['CHROM'].str.contains('B73V4_ctg')]
NILASencode = NILASencode.drop(NILASencode.columns[0:5], axis=1)
NILASencode = NILASencode.transpose()
NILASencode.to_csv(path +'/' + sys.argv[1]+'x'+sys.argv[2]+r'_Encoded_Genotypes.txt', sep='\t')       #------------->Genotype Outfile {Export}

HapMap = NILAShap.filter((Coordinates))
HapMap.to_csv(path + '/' +sys.argv[1]+'x'+sys.argv[2]+r'_Coordinates.txt', sep='\t', header= False, index= True) #------------->MAP Outfile {Export}
 
NILAS_Imputed = NILAShap.iloc[:,0:SampleNo + 2]

NILAS_Imputed = NILAS_Imputed[NILAS_Imputed.columns[NILAS_Imputed.columns.to_series().str.contains('NILAS')]]
NILASgt = NILAS_Imputed.transpose()
NILASgt.to_csv(path +'/' + sys.argv[1]+'x'+sys.argv[2]+r'_Imputed_Genotypes.txt', sep='\t')       #------------->Impute Outfile{ Export}

#Redefine  Haploytpe Blocks--------------------------------------------------------------------------------------------

for sample in NILASlist:
    NILAShap[(sample)] = NILAShap[(sample)].map({'1':'G', 'A2':'A','A':'A','C':'C','G':'G','0':'0'})

    NILAShap[(sample)+'_Haplotype_Block'] = ((NILAShap.CHROM != NILAShap.CHROM.shift(1)) |
        (NILAShap[(sample)] != ((NILAShap[(sample)].astype(str)).shift(1)))).astype(int).cumsum()
    NILAShap[(sample)+'_Start'] = NILAShap.groupby([(sample)+'_Haplotype_Block'])['POS'].transform('first')
    NILAShap[(sample)+'_Stop'] = NILAShap.groupby([(sample)+'_Haplotype_Block'])['POS'].transform('last')
   
    NILAShap[(sample)+'_Five_Flank'] = NILAShap.groupby(['CHROM',(sample)+'_Stop',(sample)+'_Haplotype_Block'])['POS'].transform('last').shift(1)
    NILAShap[(sample)+'_Three_Flank'] = NILAShap.groupby(['CHROM',(sample)+'_Start',(sample)+'_Haplotype_Block'])['POS'].transform('first').shift(-1)
    
    NILAShap.loc[(NILAShap[(sample)+'_Start'].astype(str) == NILAShap[(sample)+'_Three_Flank'].astype(str)), [(sample)+'_Three_Flank']] = 'NaN'
    NILAShap[(sample)+'_Three_Flank'] = NILAShap[(sample)+'_Three_Flank'].replace(to_replace=['NaN'], regex=True, method='bfill')
    #where start = three flank -> three flank = blank
    #back fill blanks with next value 
    
    NILAShap.loc[(NILAShap[(sample)+'_Stop'].astype(str) == NILAShap[(sample)+'_Five_Flank'].astype(str)), [(sample)+'_Five_Flank']] = 'NaN'
    NILAShap.loc[(NILAShap[(sample)+'_Five_Flank'].astype(str) == 'NaN'), [(sample)+'_Three_Flank']] = 'NaN'
    #where stop = five flank -> five flank = blank
    
    NILAShap = NILAShap.fillna(0)
    NILAShap = NILAShap.replace(to_replace=['NaN'], value=0)
    
#Phyiscal Translation Block____________________________________________________________________________________________________________________________________________________________
print("Translating HapBlock to Physical Coordinates...")
 
for sample in NILASlist:    
    NILAShap[(sample)+'_Hap_No'] = NILAShap.groupby(['CHROM'])[(sample)+'_Haplotype_Block'].transform('nunique') #<-- Number of haplotype blocks per chromosome
    
    NILAShap[(sample)+'_Five_Flank'] = np.where((NILAShap[(sample)+'_Five_Flank'].astype(int) > NILAShap[(sample)+'_Start'].astype(int)),0, NILAShap[(sample)+'_Five_Flank'])
       
    NILAShap[(sample)+'_Three_Flank'] = np.where((NILAShap[(sample)+'_Five_Flank'].astype(int) > NILAShap[(sample)+'_Three_Flank'].astype(int)),
        (NILAShap['CHROM'].map(CHROMend)), NILAShap[(sample)+'_Three_Flank'])
    
    NILAShap[(sample)+'_Three_Flank'] = np.where((NILAShap[(sample)+'_Hap_No'] == 1),(NILAShap['CHROM'].map(CHROMend)), NILAShap[(sample)+'_Three_Flank'])
    
    NILAShap[(sample)+'_Five_Flank_Dis'] = NILAShap[(sample)+'_Start'].astype(int) - NILAShap[(sample)+'_Five_Flank'].astype(int)
    
    NILAShap[(sample)+'_Five_Flank_Dis'] = np.where((NILAShap[(sample)+'_Five_Flank'] == 0),
    (NILAShap[(sample)+'_Start'].astype(int) - NILAShap[(sample)+'_Five_Flank'].astype(int)).abs(), NILAShap[(sample)+'_Five_Flank_Dis'] )
    
    NILAShap[(sample)+'_Three_Flank_Dis'] = (NILAShap[(sample)+'_Three_Flank'].astype(int) - NILAShap[(sample)+'_Stop'].astype(int))  ####<==============
    
    NILAShap[(sample)+'_5_BP_Density'] = NILAShap[(sample)+'_Five_Flank_Dis'].div(2)
    NILAShap[(sample)+'_3_BP_Density'] = NILAShap[(sample)+'_Three_Flank_Dis'].div(2)
    
    NILAShap[(sample)+'_5_Breakpoint'] =  NILAShap[(sample)+'_Start'].astype(int) -  NILAShap[(sample)+'_5_BP_Density']
    
    NILAShap[(sample)+'_3_Breakpoint'] = NILAShap[(sample)+'_Stop'].astype(int) + NILAShap[(sample)+'_3_BP_Density']

    
    NILAShap[(sample)+'_5_Breakpoint'] = np.where((NILAShap[(sample)+'_Five_Flank'] == 0) ,0, NILAShap[(sample)+'_5_Breakpoint'])   #<------- One hap block per chromo logic
    NILAShap[(sample)+'_5_Breakpoint'] = np.where((NILAShap[(sample)+'_Hap_No'] == 1) ,0, NILAShap[(sample)+'_5_Breakpoint'])
    NILAShap[(sample)+'_3_Breakpoint'] = np.where((NILAShap[(sample)+'_Hap_No'] == 1) ,0, NILAShap[(sample)+'_3_Breakpoint'])
    NILAShap[(sample)+'_3_Breakpoint'] = np.where((NILAShap[(sample)+'_Three_Flank'] == (NILAShap['CHROM'].map(CHROMend))),(NILAShap['CHROM'].map(CHROMend)), NILAShap[(sample)+'_3_Breakpoint'])
    
    NILAShap[(sample)+'_3_Breakpoint'] = np.where((NILAShap[(sample)+'_3_BP_Density'] <= 0),(NILAShap['CHROM'].map(CHROMend)), NILAShap[(sample)+'_3_Breakpoint']) ####<============== NEW LOGIC

    
    NILAShap[(sample)+'_PhyLen'] = NILAShap[(sample)+'_3_Breakpoint'].astype(int) - NILAShap[(sample)+'_5_Breakpoint'].astype(int)
    
    NILAShap[(sample)+'_Marker_No']= NILAShap.groupby(['CHROM',(sample)+'_Haplotype_Block'])['CHROM'].transform('size')

#Global Summary Block____________________________________________________________________________________________________________________________________________________________

NILAShap.set_index('CHROM', inplace=True)
NILAShap = NILAShap.drop(['POS','DP2_Loci_1','DP2_Loci_2','DP2_Loci_3','DP2_Loci_4'], axis=1)
NILAShap = NILAShap[NILAShap.columns[~NILAShap.columns.to_series().str.contains('Flank')]]

# NILAShap = NILAShap.drop(NILAShap.columns[0], axis=1)

# print(NILAShap['NILASq1g11i01s2Cp1r1'])


for i, g in NILAShap.groupby(by=lambda x: x.split('_')[0], axis=1):
    NILAS_Dict[(i)] = g.groupby(i+'_Haplotype_Block').head(1)
    NILAS_Dict[(i)].columns = NILAS_Dict[(i)].columns.str.replace(i + '_','')
    NILAS_Dict[(i)].columns = NILAS_Dict[(i)].columns.str.replace(i,'Hap_Geno')

NILAS_Values = pd.concat(NILAS_Dict.values(),ignore_index=False, keys = NILAS_Dict.keys())
NILAS_Values['Sample'] = NILAS_Values.index.get_level_values(0)
NILAS_Values['CHROM'] = NILAS_Values.index.get_level_values(1)

###_______ZMPR1______________________________________________________________________________________________###

###-------- FOREGROUND---------#####
ZmPR1 =NILAS_Values[NILAS_Values['Sample'].str.contains('1g')].copy()


ZmPR1['Scheme'] = ZmPR1['Sample'].str.extract('(g..)', expand=False)
ZmPR1['Scheme'] = ZmPR1['Scheme'].map(Scheme_dict)

ZmPR1['ZmPR'] = ZmPR1['Sample'].str.extract('q(.)', expand=False)

ZmPR1['Hap_Geno'] = ZmPR1['Hap_Geno'].map(Geno_dict)

ZmPR1['5_Breakpoint'] = ZmPR1['5_Breakpoint'].astype(float)
ZmPR1['3_Breakpoint'] = ZmPR1['3_Breakpoint'].astype(float)

ZmPR1['ZmPR+'] = np.where(((ZmPR1['5_Breakpoint'].astype(float)).between(57191440, 185363076, inclusive=True) & (ZmPR1['CHROM'] ==  '1') |
     (ZmPR1['3_Breakpoint'].astype(float).between(57191440, 185363076, inclusive=True)) & (ZmPR1['CHROM'] ==  '1')),1,0)

ZmPR1['Foreground_Start'] = np.where((ZmPR1['ZmPR+'] == 1) & (ZmPR1['5_Breakpoint'].astype(int) < 57191440 ), 57191440, ZmPR1['5_Breakpoint'])
ZmPR1['Foreground_Start'] = np.where((ZmPR1['ZmPR+'] != 1), 0, ZmPR1['Foreground_Start'])

ZmPR1['Foreground_Stop'] = np.where((ZmPR1['ZmPR+'] == 1) & (ZmPR1['3_Breakpoint'].astype(int) > 185363076 ), 185363076, ZmPR1['3_Breakpoint'])
ZmPR1['Foreground_Stop'] = np.where((ZmPR1['ZmPR+'] != 1), 0, ZmPR1['Foreground_Stop'])
 
ZmPR1['Foreground_Perc'] = np.where((ZmPR1['ZmPR+'] == 1),
    ((ZmPR1['Foreground_Stop'].astype(int) - ZmPR1['Foreground_Start'].astype(int)).div(B73len)*100),0)

ZmPR1['FG_5_BP_Density']= np.where((ZmPR1['ZmPR+'] == 1),ZmPR1['5_BP_Density'],0)
ZmPR1['FG_3_BP_Density']= np.where((ZmPR1['ZmPR+'] == 1),ZmPR1['3_BP_Density'],0)

###-------- Background---------#####

ZmPR1['BG_5_BP_Density']= np.where((ZmPR1['ZmPR+'] != 1),ZmPR1['5_BP_Density'],0)
ZmPR1['BG_3_BP_Density']= np.where((ZmPR1['ZmPR+'] != 1),ZmPR1['3_BP_Density'],0)

#Starts in background; stops in ZmPR scenerio
ZmPR1['Background_Start'] = np.where(((ZmPR1['ZmPR+'] == 1) & (ZmPR1['5_Breakpoint'].astype(int) < 57191440 ) &
    (ZmPR1['3_Breakpoint'].astype(int) < 185363076 )),(ZmPR1['5_Breakpoint']), (0))
ZmPR1['Background_Stop'] = np.where(((ZmPR1['ZmPR+'] == 1) & (ZmPR1['5_Breakpoint'].astype(int) < 57191440 ) &
    (ZmPR1['3_Breakpoint'].astype(int) < 185363076 )),(57191440), (0))

#Starts in ZmPR; stops in background scenerio
ZmPR1['Background_Start'] = np.where(((ZmPR1['ZmPR+'] == 1) & (ZmPR1['5_Breakpoint'].astype(int) < 185363076 ) &
    (ZmPR1['3_Breakpoint'].astype(int) > 185363076 )),(185363076), (ZmPR1['Background_Start'])) 
ZmPR1['Background_Stop'] = np.where(((ZmPR1['ZmPR+'] == 1) & (ZmPR1['5_Breakpoint'].astype(float)).between(57191440, 185363076, inclusive=True) &
    (ZmPR1['3_Breakpoint'].astype(int) > 185363076 )),(ZmPR1['3_Breakpoint']), (ZmPR1['Background_Stop']))

#Starts and stops in ZmPR scenerio
ZmPR1['Background_Start'] = np.where(((ZmPR1['ZmPR+'] == 1) & (ZmPR1['5_Breakpoint'].astype(int) > 57191440 ) &
    (ZmPR1['3_Breakpoint'].astype(int) < 185363076 )),(0), (ZmPR1['Background_Start']))
ZmPR1['Background_Stop'] = np.where(((ZmPR1['ZmPR+'] == 1) & (ZmPR1['5_Breakpoint'].astype(int) > 57191440 ) &
    (ZmPR1['3_Breakpoint'].astype(int) < 185363076 )),(0), (ZmPR1['Background_Stop']))

#Starts and stops in background scenerio
ZmPR1['Background_Start'] = np.where((ZmPR1['ZmPR+'] != 1), (ZmPR1['5_Breakpoint']), (ZmPR1['Background_Start']))
ZmPR1['Background_Stop'] = np.where((ZmPR1['ZmPR+'] != 1), (ZmPR1['3_Breakpoint']), (ZmPR1['Background_Stop']))

ZmPR1['Background_Perc'] = (ZmPR1['Background_Stop'].astype(int) - ZmPR1['Background_Start'].astype(int)).div(B73len)*100

# ZmPR1.to_csv(os.path.join(path,sys.argv[1]+'x'+sys.argv[2]+r'_REVIEW2.txt'), sep='\t', header=True,index=True) #-------------> INQUIRY FILE {Export}


ZmPR1 = ZmPR1[ZmPR1.columns[~ZmPR1.columns.to_series().str.contains('Start')]]
ZmPR1 = ZmPR1[ZmPR1.columns[~ZmPR1.columns.to_series().str.contains('Stop')]]
ZmPR1 = ZmPR1[ZmPR_list]

FG1 = ['Sample','Scheme', 'ZmPR','Hap_Geno','Foreground_Perc']
ZmPR1_FG_Geno = ZmPR1.filter((FG1))
ZmPR1_FG_Geno['Hap_Geno'] = ZmPR1_FG_Geno['Hap_Geno'].map({'DP':'DP_FG','RP':'RP_FG', 'NA':'NA_FG','Het':'Het_FG'})


BG1 = ['Sample','Scheme', 'ZmPR','Hap_Geno','Background_Perc']
ZmPR1_BG_Geno = ZmPR1.filter((BG1))
ZmPR1_BG_Geno['Hap_Geno'] = ZmPR1_BG_Geno['Hap_Geno'].map({'DP':'DP_BG','RP':'RP_BG', 'NA':'NA_BG','Het':'Het_BG'})


ZmPR1_Geno_Sum = pd.concat([ZmPR1_FG_Geno,ZmPR1_BG_Geno], axis=0,sort=True)
# ZmPR1_Geno_Sum = ZmPR1_Geno_Sum[~ZmPR1_Geno_Sum['Sample'].str.contains('C')]
ZmPR1_Geno_Sum = ZmPR1_Geno_Sum.fillna(0)
ZmPR1_Geno_Sum['Perc_Geno'] = ZmPR1_Geno_Sum['Foreground_Perc'] + ZmPR1_Geno_Sum['Background_Perc']

ZmPR1_Geno_Sum = ZmPR1_Geno_Sum.pivot_table(index=['Sample','Scheme','ZmPR'], columns=['Hap_Geno'], values='Perc_Geno', aggfunc='sum')
ZmPR1_Geno_Sum = ZmPR1_Geno_Sum.fillna(0)
# print(ZmPR1)

ZmPR1I = ZmPR1[~ZmPR1['Sample'].str.contains('C')]
ZmPR1I.to_csv(path + '/' + sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR1_SUM.txt', sep='\t', header=True,index=False)#------------->Global {Export}
ZmPR1C = ZmPR1[ZmPR1['Sample'].str.contains('C')]
ZmPR1C.to_csv(path + '/' + sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR1_SUMC.txt', sep='\t', header=True,index=False)#------------->Global {Export}
ZmPR1_Geno_Sum.to_csv(path + '/' +sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR1_Geno_Sum.txt', sep='\t', header=True,index=True) #------------->Genotype Percentages {Export}
##_______ZMPR2______________________________________________________________________________________________###

###-------- FOREGROUND---------#####
ZmPR2 =NILAS_Values[NILAS_Values['Sample'].str.contains('2g')].copy()

ZmPR2['Scheme'] = ZmPR2['Sample'].str.extract('(g..)', expand=False)
ZmPR2['Scheme'] = ZmPR2['Scheme'].map(Scheme_dict)

ZmPR2['ZmPR'] = ZmPR2['Sample'].str.extract('q(.)', expand=False)

ZmPR2['Hap_Geno'] = ZmPR2['Hap_Geno'].map(Geno_dict)

ZmPR2['5_Breakpoint'] = ZmPR2['5_Breakpoint'].astype(float)
ZmPR2['3_Breakpoint'] = ZmPR2['3_Breakpoint'].astype(float)

ZmPR2['ZmPR+'] = np.where(((ZmPR2['5_Breakpoint'].astype(float)).between(97682806, 167995313, inclusive=True) & (ZmPR2['CHROM'] ==  '8') |
     (ZmPR2['3_Breakpoint'].astype(float).between(97682806, 167995313, inclusive=True)) & (ZmPR2['CHROM'] ==  '8')),1,0)

ZmPR2['Foreground_Start'] = np.where((ZmPR2['ZmPR+'] == 1) & (ZmPR2['5_Breakpoint'].astype(int) < 97682806 ), 97682806, ZmPR2['5_Breakpoint'])
ZmPR2['Foreground_Start'] = np.where((ZmPR2['ZmPR+'] != 1), 0, ZmPR2['Foreground_Start'])

ZmPR2['Foreground_Stop'] = np.where((ZmPR2['ZmPR+'] == 1) & (ZmPR2['3_Breakpoint'].astype(int) > 167995313 ), 167995313, ZmPR2['3_Breakpoint'])
ZmPR2['Foreground_Stop'] = np.where((ZmPR2['ZmPR+'] != 1), 0, ZmPR2['Foreground_Stop'])
 
ZmPR2['Foreground_Perc'] = np.where((ZmPR2['ZmPR+'] == 1),
    ((ZmPR2['Foreground_Stop'].astype(int) - ZmPR2['Foreground_Start'].astype(int)).div(B73len)*100),0)

ZmPR2['FG_5_BP_Density']= np.where((ZmPR2['ZmPR+'] == 1),ZmPR2['5_BP_Density'],0)
ZmPR2['FG_3_BP_Density']= np.where((ZmPR2['ZmPR+'] == 1),ZmPR2['3_BP_Density'],0)

###-------- Background---------#####

ZmPR2['BG_5_BP_Density']= np.where((ZmPR2['ZmPR+'] != 1),ZmPR2['5_BP_Density'],0)
ZmPR2['BG_3_BP_Density']= np.where((ZmPR2['ZmPR+'] != 1),ZmPR2['3_BP_Density'],0)

#Starts in background; stops in ZmPR scenario
ZmPR2['Background_Start'] = np.where(((ZmPR2['ZmPR+'] == 1) & (ZmPR2['5_Breakpoint'].astype(int) < 97682806 ) &
    (ZmPR2['3_Breakpoint'].astype(int) < 167995313 )),(ZmPR2['5_Breakpoint']), (0))
ZmPR2['Background_Stop'] = np.where(((ZmPR2['ZmPR+'] == 1) & (ZmPR2['5_Breakpoint'].astype(int) < 97682806 ) &
    (ZmPR2['3_Breakpoint'].astype(int) < 167995313 )),(97682806), (0))

#Starts in ZmPR; stops in background scenario
ZmPR2['Background_Start'] = np.where(((ZmPR2['ZmPR+'] == 1) & (ZmPR2['5_Breakpoint'].astype(int) < 167995313 ) &
    (ZmPR2['3_Breakpoint'].astype(int) > 167995313 )),(167995313), (ZmPR2['Background_Start'])) 
ZmPR2['Background_Stop'] = np.where(((ZmPR2['ZmPR+'] == 1) & (ZmPR2['5_Breakpoint'].astype(float)).between(97682806, 167995313, inclusive=True) &
    (ZmPR2['3_Breakpoint'].astype(int) > 167995313 )),(ZmPR2['3_Breakpoint']), (ZmPR2['Background_Stop']))

#Starts and stops in ZmPR scenario
ZmPR2['Background_Start'] = np.where(((ZmPR2['ZmPR+'] == 1) & (ZmPR2['5_Breakpoint'].astype(int) > 97682806 ) &
    (ZmPR2['3_Breakpoint'].astype(int) < 167995313 )),(0), (ZmPR2['Background_Start']))
ZmPR2['Background_Stop'] = np.where(((ZmPR2['ZmPR+'] == 1) & (ZmPR2['5_Breakpoint'].astype(int) > 97682806 ) &
    (ZmPR2['3_Breakpoint'].astype(int) < 167995313 )),(0), (ZmPR2['Background_Stop']))

#Starts and stops in background scenario
ZmPR2['Background_Start'] = np.where((ZmPR2['ZmPR+'] != 1), (ZmPR2['5_Breakpoint']), (ZmPR2['Background_Start']))
ZmPR2['Background_Stop'] = np.where((ZmPR2['ZmPR+'] != 1), (ZmPR2['3_Breakpoint']), (ZmPR2['Background_Stop']))

ZmPR2['Background_Perc'] = (ZmPR2['Background_Stop'].astype(int) - ZmPR2['Background_Start'].astype(int)).div(B73len)*100

ZmPR2 = ZmPR2[ZmPR2.columns[~ZmPR2.columns.to_series().str.contains('Start')]]
ZmPR2 = ZmPR2[ZmPR2.columns[~ZmPR2.columns.to_series().str.contains('Stop')]]
ZmPR2 = ZmPR2[ZmPR_list]

FG2 = ['Sample','Scheme', 'ZmPR','Hap_Geno','Foreground_Perc']
ZmPR2_FG_Geno = ZmPR2.filter((FG2))
ZmPR2_FG_Geno['Hap_Geno'] = ZmPR2_FG_Geno['Hap_Geno'].map({'DP':'DP_FG','RP':'RP_FG', 'NA':'NA_FG','Het':'Het_FG'})


BG2 = ['Sample','Scheme', 'ZmPR','Hap_Geno','Background_Perc']
ZmPR2_BG_Geno = ZmPR2.filter((BG2))
ZmPR2_BG_Geno['Hap_Geno'] = ZmPR2_BG_Geno['Hap_Geno'].map({'DP':'DP_BG','RP':'RP_BG', 'NA':'NA_BG','Het':'Het_BG'})


ZmPR2_Geno_Sum = pd.concat([ZmPR2_FG_Geno,ZmPR2_BG_Geno], axis=0)
# ZmPR2_Geno_Sum = ZmPR2_Geno_Sum[~ZmPR2_Geno_Sum['Sample'].str.contains('C')]
ZmPR2_Geno_Sum = ZmPR2_Geno_Sum.fillna(0)
ZmPR2_Geno_Sum['Perc_Geno'] = ZmPR2_Geno_Sum['Foreground_Perc'] + ZmPR2_Geno_Sum['Background_Perc']

ZmPR2_Geno_Sum = ZmPR2_Geno_Sum.pivot_table(index=['Sample','Scheme','ZmPR'], columns=['Hap_Geno'], values='Perc_Geno', aggfunc='sum')
ZmPR2_Geno_Sum = ZmPR2_Geno_Sum.fillna(0)

ZmPR2I = ZmPR2[~ZmPR2['Sample'].str.contains('C')]
ZmPR2I.to_csv(path + '/' +sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR2_SUM.txt', sep='\t', header=True,index=False)#------------->Global {Export}
ZmPR2C = ZmPR2[ZmPR2['Sample'].str.contains('C')]
ZmPR2C.to_csv(path + '/' +sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR2_SUMC.txt', sep='\t', header=True,index=False)#------------->Global {Export}
ZmPR2_Geno_Sum.to_csv(path + '/' +sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR2_Geno_Sum.txt', sep='\t', header=True,index=True) #------------->Genotype Percentages {Export} {Export}

###_______ZMPR3______________________________________________________________________________________________###

###-------- FOREGROUND---------#####
ZmPR3 =NILAS_Values[NILAS_Values['Sample'].str.contains('3g')].copy()

ZmPR3['Scheme'] = ZmPR3['Sample'].str.extract('(g..)', expand=False)
ZmPR3['Scheme'] = ZmPR3['Scheme'].map(Scheme_dict)

ZmPR3['ZmPR'] = ZmPR3['Sample'].str.extract('q(.)', expand=False)

ZmPR3['Hap_Geno'] = ZmPR3['Hap_Geno'].map(Geno_dict)

ZmPR3['5_Breakpoint'] = ZmPR3['5_Breakpoint'].astype(float)
ZmPR3['3_Breakpoint'] = ZmPR3['3_Breakpoint'].astype(float)

ZmPR3['ZmPR+'] = np.where(((ZmPR3['5_Breakpoint'].astype(float)).between(7384635, 119857067, inclusive=True) & (ZmPR3['CHROM'] ==  '9') |
     (ZmPR3['3_Breakpoint'].astype(float).between(7384635, 119857067, inclusive=True)) & (ZmPR3['CHROM'] ==  '9')),1,0)

ZmPR3['Foreground_Start'] = np.where((ZmPR3['ZmPR+'] == 1) & (ZmPR3['5_Breakpoint'].astype(int) < 7384635 ), 7384635, ZmPR3['5_Breakpoint'])
ZmPR3['Foreground_Start'] = np.where((ZmPR3['ZmPR+'] != 1), 0, ZmPR3['Foreground_Start'])

ZmPR3['Foreground_Stop'] = np.where((ZmPR3['ZmPR+'] == 1) & (ZmPR3['3_Breakpoint'].astype(int) > 119857067 ), 119857067, ZmPR3['3_Breakpoint'])
ZmPR3['Foreground_Stop'] = np.where((ZmPR3['ZmPR+'] != 1), 0, ZmPR3['Foreground_Stop'])
 
ZmPR3['Foreground_Perc'] = np.where((ZmPR3['ZmPR+'] == 1),
    ((ZmPR3['Foreground_Stop'].astype(int) - ZmPR3['Foreground_Start'].astype(int)).div(B73len)*100),0)

ZmPR3['FG_5_BP_Density']= np.where((ZmPR3['ZmPR+'] == 1),ZmPR3['5_BP_Density'],0)
ZmPR3['FG_3_BP_Density']= np.where((ZmPR3['ZmPR+'] == 1),ZmPR3['3_BP_Density'],0)

###-------- Background---------#####

ZmPR3['BG_5_BP_Density']= np.where((ZmPR3['ZmPR+'] != 1),ZmPR3['5_BP_Density'],0)
ZmPR3['BG_3_BP_Density']= np.where((ZmPR3['ZmPR+'] != 1),ZmPR3['3_BP_Density'],0)

#Starts in background; stops in ZmPR scenario
ZmPR3['Background_Start'] = np.where(((ZmPR3['ZmPR+'] == 1) & (ZmPR3['5_Breakpoint'].astype(int) < 7384635 ) &
    (ZmPR3['3_Breakpoint'].astype(int) < 119857067 )),(ZmPR3['5_Breakpoint']), (0))
ZmPR3['Background_Stop'] = np.where(((ZmPR3['ZmPR+'] == 1) & (ZmPR3['5_Breakpoint'].astype(int) < 7384635 ) &
    (ZmPR3['3_Breakpoint'].astype(int) < 119857067 )),(7384635), (0))

#Starts in ZmPR; stops in background scenario
ZmPR3['Background_Start'] = np.where(((ZmPR3['ZmPR+'] == 1) & (ZmPR3['5_Breakpoint'].astype(int) < 119857067 ) &
    (ZmPR3['3_Breakpoint'].astype(int) > 119857067 )),(119857067), (ZmPR3['Background_Start'])) 
ZmPR3['Background_Stop'] = np.where(((ZmPR3['ZmPR+'] == 1) & (ZmPR3['5_Breakpoint'].astype(float)).between(7384635, 119857067, inclusive=True) &
    (ZmPR3['3_Breakpoint'].astype(int) > 119857067 )),(ZmPR3['3_Breakpoint']), (ZmPR3['Background_Stop']))

#Starts and stops in ZmPR scenario
ZmPR3['Background_Start'] = np.where(((ZmPR3['ZmPR+'] == 1) & (ZmPR3['5_Breakpoint'].astype(int) > 7384635 ) &
    (ZmPR3['3_Breakpoint'].astype(int) < 119857067 )),(0), (ZmPR3['Background_Start']))
ZmPR3['Background_Stop'] = np.where(((ZmPR3['ZmPR+'] == 1) & (ZmPR3['5_Breakpoint'].astype(int) > 7384635 ) &
    (ZmPR3['3_Breakpoint'].astype(int) < 119857067 )),(0), (ZmPR3['Background_Stop']))

#Starts and stops in background scenario
ZmPR3['Background_Start'] = np.where((ZmPR3['ZmPR+'] != 1), (ZmPR3['5_Breakpoint']), (ZmPR3['Background_Start']))
ZmPR3['Background_Stop'] = np.where((ZmPR3['ZmPR+'] != 1), (ZmPR3['3_Breakpoint']), (ZmPR3['Background_Stop']))

ZmPR3['Background_Perc'] = (ZmPR3['Background_Stop'].astype(int) - ZmPR3['Background_Start'].astype(int)).div(B73len)*100

ZmPR3 = ZmPR3[ZmPR3.columns[~ZmPR3.columns.to_series().str.contains('Start')]]
ZmPR3 = ZmPR3[ZmPR3.columns[~ZmPR3.columns.to_series().str.contains('Stop')]]
ZmPR3 = ZmPR3[ZmPR_list]

FG3 = ['Sample','Scheme', 'ZmPR','Hap_Geno','Foreground_Perc']
ZmPR3_FG_Geno = ZmPR3.filter((FG3))
ZmPR3_FG_Geno['Hap_Geno'] = ZmPR3_FG_Geno['Hap_Geno'].map({'DP':'DP_FG','RP':'RP_FG', 'NA':'NA_FG','Het':'Het_FG'})


BG3 = ['Sample','Scheme', 'ZmPR','Hap_Geno','Background_Perc']
ZmPR3_BG_Geno = ZmPR3.filter((BG3))
ZmPR3_BG_Geno['Hap_Geno'] = ZmPR3_BG_Geno['Hap_Geno'].map({'DP':'DP_BG','RP':'RP_BG', 'NA':'NA_BG','Het':'Het_BG'})


ZmPR3_Geno_Sum = pd.concat([ZmPR3_FG_Geno,ZmPR3_BG_Geno], axis=0)
# ZmPR3_Geno_Sum = ZmPR3_Geno_Sum[~ZmPR3_Geno_Sum['Sample'].str.contains('C')]
ZmPR3_Geno_Sum = ZmPR3_Geno_Sum.fillna(0)
ZmPR3_Geno_Sum['Perc_Geno'] = ZmPR3_Geno_Sum['Foreground_Perc'] + ZmPR3_Geno_Sum['Background_Perc']

ZmPR3_Geno_Sum = ZmPR3_Geno_Sum.pivot_table(index=['Sample','Scheme','ZmPR'], columns=['Hap_Geno'], values='Perc_Geno', aggfunc='sum')
ZmPR3_Geno_Sum = ZmPR3_Geno_Sum.fillna(0)

ZmPR3I = ZmPR3[~ZmPR3['Sample'].str.contains('C')]
ZmPR3I.to_csv(path + '/' +sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR3_SUM.txt', sep='\t', header=True,index=False)#------------->Global {Export}
ZmPR3C = ZmPR3[ZmPR3['Sample'].str.contains('C')]
ZmPR3C.to_csv(path + '/' +sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR3_SUMC.txt', sep='\t', header=True,index=False)#------------->Global {Export}
ZmPR3_Geno_Sum.to_csv(path + '/' +sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR3_Geno_Sum.txt', sep='\t', header=True,index=True) #------------->Genotype Percentages {Export}

###_______ZMPR4______________________________________________________________________________________________###

###-------- FOREGROUND---------#####
ZmPR4 =NILAS_Values[NILAS_Values['Sample'].str.contains('4g')].copy()

ZmPR4['Scheme'] = ZmPR4['Sample'].str.extract('(g..)', expand=False)
ZmPR4['Scheme'] = ZmPR4['Scheme'].map(Scheme_dict)

ZmPR4['ZmPR'] = ZmPR4['Sample'].str.extract('q(.)', expand=False)

ZmPR4['Hap_Geno'] = ZmPR4['Hap_Geno'].map(Geno_dict)

ZmPR4['5_Breakpoint'] = ZmPR4['5_Breakpoint'].astype(float)
ZmPR4['3_Breakpoint'] = ZmPR4['3_Breakpoint'].astype(float)

ZmPR4['ZmPR+'] = np.where(((ZmPR4['5_Breakpoint'].astype(float)).between(10126632, 136329803, inclusive=True) & (ZmPR4['CHROM'] ==  '10') |
     (ZmPR4['3_Breakpoint'].astype(float).between(10126632, 136329803, inclusive=True)) & (ZmPR4['CHROM'] ==  '10')),1,0)

ZmPR4['Foreground_Start'] = np.where((ZmPR4['ZmPR+'] == 1) & (ZmPR4['5_Breakpoint'].astype(int) < 10126632 ), 10126632, ZmPR4['5_Breakpoint'])
ZmPR4['Foreground_Start'] = np.where((ZmPR4['ZmPR+'] != 1), 0, ZmPR4['Foreground_Start'])

ZmPR4['Foreground_Stop'] = np.where((ZmPR4['ZmPR+'] == 1) & (ZmPR4['3_Breakpoint'].astype(int) > 136329803 ), 136329803, ZmPR4['3_Breakpoint'])
ZmPR4['Foreground_Stop'] = np.where((ZmPR4['ZmPR+'] != 1), 0, ZmPR4['Foreground_Stop'])
 
ZmPR4['Foreground_Perc'] = np.where((ZmPR4['ZmPR+'] == 1),
    ((ZmPR4['Foreground_Stop'].astype(int) - ZmPR4['Foreground_Start'].astype(int)).div(B73len)*100),0)

ZmPR4['FG_5_BP_Density']= np.where((ZmPR4['ZmPR+'] == 1),ZmPR4['5_BP_Density'],0)
ZmPR4['FG_3_BP_Density']= np.where((ZmPR4['ZmPR+'] == 1),ZmPR4['3_BP_Density'],0)

###-------- Background---------#####

ZmPR4['BG_5_BP_Density']= np.where((ZmPR4['ZmPR+'] != 1),ZmPR4['5_BP_Density'],0)
ZmPR4['BG_3_BP_Density']= np.where((ZmPR4['ZmPR+'] != 1),ZmPR4['3_BP_Density'],0)

#Starts in background; stops in ZmPR scenario
ZmPR4['Background_Start'] = np.where(((ZmPR4['ZmPR+'] == 1) & (ZmPR4['5_Breakpoint'].astype(int) < 10126632 ) &
    (ZmPR4['3_Breakpoint'].astype(int) < 136329803 )),(ZmPR4['5_Breakpoint']), (0))
ZmPR4['Background_Stop'] = np.where(((ZmPR4['ZmPR+'] == 1) & (ZmPR4['5_Breakpoint'].astype(int) < 10126632 ) &
    (ZmPR4['3_Breakpoint'].astype(int) < 136329803 )),(10126632), (0))

#Starts in ZmPR; stops in background scenario
ZmPR4['Background_Start'] = np.where(((ZmPR4['ZmPR+'] == 1) & (ZmPR4['5_Breakpoint'].astype(int) < 136329803 ) &
    (ZmPR4['3_Breakpoint'].astype(int) > 136329803 )),(136329803), (ZmPR4['Background_Start'])) 
ZmPR4['Background_Stop'] = np.where(((ZmPR4['ZmPR+'] == 1) & (ZmPR4['5_Breakpoint'].astype(float)).between(10126632, 136329803, inclusive=True) &
    (ZmPR4['3_Breakpoint'].astype(int) > 136329803 )),(ZmPR4['3_Breakpoint']), (ZmPR4['Background_Stop']))

#Starts and stops in ZmPR scenario
ZmPR4['Background_Start'] = np.where(((ZmPR4['ZmPR+'] == 1) & (ZmPR4['5_Breakpoint'].astype(int) > 10126632 ) &
    (ZmPR4['3_Breakpoint'].astype(int) < 136329803 )),(0), (ZmPR4['Background_Start']))
ZmPR4['Background_Stop'] = np.where(((ZmPR4['ZmPR+'] == 1) & (ZmPR4['5_Breakpoint'].astype(int) > 10126632 ) &
    (ZmPR4['3_Breakpoint'].astype(int) < 136329803 )),(0), (ZmPR4['Background_Stop']))

#Starts and stops in background scenario
ZmPR4['Background_Start'] = np.where((ZmPR4['ZmPR+'] != 1), (ZmPR4['5_Breakpoint']), (ZmPR4['Background_Start']))
ZmPR4['Background_Stop'] = np.where((ZmPR4['ZmPR+'] != 1), (ZmPR4['3_Breakpoint']), (ZmPR4['Background_Stop']))

ZmPR4['Background_Perc'] = (ZmPR4['Background_Stop'].astype(int) - ZmPR4['Background_Start'].astype(int)).div(B73len)*100

ZmPR4 = ZmPR4[ZmPR4.columns[~ZmPR4.columns.to_series().str.contains('Start')]]
ZmPR4 = ZmPR4[ZmPR4.columns[~ZmPR4.columns.to_series().str.contains('Stop')]]
ZmPR4 = ZmPR4[ZmPR_list]

FG4 = ['Sample','Scheme', 'ZmPR','Hap_Geno','Foreground_Perc']
ZmPR4_FG_Geno = ZmPR4.filter((FG4))
ZmPR4_FG_Geno['Hap_Geno'] = ZmPR4_FG_Geno['Hap_Geno'].map({'DP':'DP_FG','RP':'RP_FG', 'NA':'NA_FG','Het':'Het_FG'})


BG4 = ['Sample','Scheme', 'ZmPR','Hap_Geno','Background_Perc']
ZmPR4_BG_Geno = ZmPR4.filter((BG4))
ZmPR4_BG_Geno['Hap_Geno'] = ZmPR4_BG_Geno['Hap_Geno'].map({'DP':'DP_BG','RP':'RP_BG', 'NA':'NA_BG','Het':'Het_BG'})


ZmPR4_Geno_Sum = pd.concat([ZmPR4_FG_Geno,ZmPR4_BG_Geno], axis=0)
# ZmPR4_Geno_Sum = ZmPR4_Geno_Sum[~ZmPR4_Geno_Sum['Sample'].str.contains('C')]
ZmPR4_Geno_Sum = ZmPR4_Geno_Sum.fillna(0)
ZmPR4_Geno_Sum['Perc_Geno'] = ZmPR4_Geno_Sum['Foreground_Perc'] + ZmPR4_Geno_Sum['Background_Perc']

ZmPR4_Geno_Sum = ZmPR4_Geno_Sum.pivot_table(index=['Sample','Scheme','ZmPR'], columns=['Hap_Geno'], values='Perc_Geno', aggfunc='sum')
ZmPR4_Geno_Sum = ZmPR4_Geno_Sum.fillna(0)

ZmPR4I = ZmPR4[~ZmPR4['Sample'].str.contains('C')]
ZmPR4I.to_csv(path + '/' +sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR4_SUM.txt', sep='\t', header=True,index=False)#------------->Global {Export}
ZmPR4C = ZmPR4[ZmPR4['Sample'].str.contains('C')]
ZmPR4C.to_csv(path + '/' +sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR4_SUMC.txt', sep='\t', header=True,index=False)#------------->Global {Export}
ZmPR4_Geno_Sum.to_csv(path + '/' +sys.argv[1]+'x'+sys.argv[2]+r'_ZmPR4_Geno_Sum.txt', sep='\t', header=True,index=True) #------------->Genotype Percentages {Export}


print ('!!!!----Script Complete----!!!!')
