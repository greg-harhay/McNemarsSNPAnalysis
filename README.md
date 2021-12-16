# McNemar's SNP Analysis
Matlab scripts for performing McNemars Test on SNPs of the Illumina Bovine HD770 Chip. A summary of the McNemar's test can be found [here](https://en.wikipedia.org/wiki/McNemar%27s_test).  These scripts, written in 100% MatLab, are the computational engine that generated the results presented in the publication titled "Association of ARRDC3 and NFIA variants with bovine congestive heart failure in feedlot cattle" with link HERE
 
##  PLINK PED File Grooming
Each row represents a single animal's genotype in the PLINK PED file. The PED file must consist of genotypes from matched pairs of animals, with a single case animal matched with a single control animal and designated as a matched pair with a pair identifier.  The PED file must have the pair identifier, a number, in first column. The case and control phenotype (PHENO) are coded in column 6, with the  control coded with 1 (PHENO = 1) and the case with a 2 (PHENO = 2). For each match pair, the control animal shall precede the case animal, top to bottom, in the PED file.    

Use UNIX/Linux **sort** to unscramble PED. The command below takes an unsorted PED of matched case and control animals and sorts first by number in column 1 (pair identifier) and the by column 6 (control ID preceding the case ID)

`sort -k1,1n -k6,6n  BCHF102pairsHD770Filtered.ped > BCHF102pairsHD770FilteredSorted.ped`

To make it easier to check that the sorting worked, use Unix/Linux **awk** to extract the first 6 columns of the PED to create a FAM file of animal metadata.

`awk '{print $1,$2,$3,$4,$5,$6}'  BCHF102pairsHD770FilteredSorted.ped >  BCHF102pairsHD770FilteredSorted.fam`

The FAM file is included in this repository.

### Mike Heaton What are other PED filtering requirements ?
 
## Input Requirements
* PED file groomed as described above
* MAP file of genome coordinates for each SNP in PED

## Running McNemar's Analysis: Chromosome 28 example
McNemars\_Chip\_Analysis.mlx - Matlab Live Script - currently configured to run SNP resident only on chr 28 to demonstrate running the script and output. Can easily change analysis input files within McNemars\_Chip\_Analysis.mlx by providing file names for the PED, MAP and basename for the  plinkPED, plinkMAP, and basename file handles within this script. This script can be run on a local workstation, assuming that the Matlab function files (*.m) in this repository are in the users path at runtime. Please make sure that the PED and MPA files are in the same directory as the McNemars\_Chip\_Analysis.mlx Live Script, if not, please change path to the files for the file handles. Alternatively, this analysis can be run in the [McNemarsSNP Analysis CodeOcean (CO) Compute Capsule](https://codeocean.com/capsule/0870729/tree) using McNemars\_Chip\_Analysis.m **function**, a version of the Live Script modified to run in CO. *(Note this link won't work for you unless you given permission by Greg or this Compute Capsule is made public ... it is currently private until the paper is published)*
###  Evaluating McNemar's contingency table


                                      Control
                             +           |         -
                     +                   |
                             Qa          |         Qb
           Case   --------------------------------------------
                                         |
                             Qc          |         Qd
                                         |
                     -                   |
                                         |



### How to score
* Qa if allele in both case and control, Qa counter increments by 1
* Qb if case has allele, but not control, Qb counter increments by 1
* Qc if case does not have allele, but control does, Qc counter increments by 1
* Qd if neither control or case have allele, then Qd increments by 1

**Only Qb and Qc are informative for risk or protection**


## Outputs
* CSV file of McNemar's test scores, occupancy of McNeamar's contingency table quandrants,  chi-square, chi-square continuity correction , exact p-values, & mid p-values
* pdf representation of MatLab Livescript (.mlx.pdf) 


### CSV output file: Explanation of column header abbreviations

HiFreqAllele = high frequency allele

For ExactlyOneAllele case with high frequency allele

* ExactOneHF_a = exactly one allele (heterozygotes, hets) only) in quadrant Qa  
* ExactOneHF_b = exactly one allele (hets only) in quadrant Qb
* ExactOneHF_c = exactly one allele (hets only) in quadrant Qc
* ExactOneHF_d = exactly one allele (hets only) in quadrant Qd
* ExactOneHF_bPc = Qb + Qc
* ExactOneHF_aPbPcPd = Qa + Qb + Qc + Qd
* ExactOneHF_ bPC_ D_ aPbPcPd = (Qa + Qc) / (Qa + Qb + Qc + Qd)
* ExactOneHF_bDc = Qb / Qc
* ExactOneHF_ChiSquared = ((Qb - Qc)^2)/(Qb + Qc)
* ExactOneHF\_ChiSquared_CC =  (abs(Qb-Qc)-1)^2/(Qb + Qc)
* ExactOneHF\_ p\_exact =  exact p-value
* neg\_log10\_ExactOneHF\_p_exact =  -log10(exact p-value)
* ExactOneHF\_p\_mid = mid p-value
* neg\_log10\_ExactOneHF\_p\_mid =  -log10(mid p-value)

The abbreviations schema is the same for the ExactlyOneOrTwoAllele (hets or homozygotes) and ExactlyTwoAllele (homozygotes) cases using both the high and low frequency alleles.


##  ![](https://unlicense.org/pd-icon.png) License
This code is released into the public domain under the [UnLicence](https://unlicense.org) 

	This is free and unencumbered software released into the public domain.
	
	Anyone is free to copy, modify, publish, use, compile, sell, or
	distribute this software, either in source code form or as a compiled
	binary, for any purpose, commercial or non-commercial, and by any
	means.
	
	In jurisdictions that recognize copyright laws, the author or authors
	of this software dedicate any and all copyright interest in the
	software to the public domain. We make this dedication for the benefit
	of the public at large and to the detriment of our heirs and
	successors. We intend this dedication to be an overt act of
	relinquishment in perpetuity of all present and future rights to this
	software under copyright law.
	
	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
	MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
	IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
	OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
	ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
	OTHER DEALINGS IN THE SOFTWARE.
	
	For more information, please refer to <http://unlicense.org/>
 
