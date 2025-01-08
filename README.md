# McNemar's SNP Analysis 
Matlab scripts for performing McNemars Test on SNPs of the Illumina Bovine HD770 Chip. A summary of the McNemar's test can be found [here](https://en.wikipedia.org/wiki/McNemar%27s_test).  These scripts, written in 100% MatLab, are the computational engine that generated the results presented in the publication titled "Association of *ARRDC3* and *NFIA* variants with bovine congestive heart failure in feedlot cattle" with link HERE
## Code, Data, Results 
Our CodeOcean Compute capsule can be found [here](https://codeocean.com/capsule/4541362/tree/v1).

##  PLINK PED File Grooming
Each [PLINK](https://zzz.bwh.harvard.edu/plink/download.shtml#download) PED file row represents a single animal's six metadata columns, plus their diploid genotype columns, in the order of the SNPs present in the PLINK MAP file.   The .ped file should be modified to include the pair identifier (a number) in the first column, and the rows sorted in descending pair order from 1 to 102 with the unaffected animal (i.e., “control”) listed first in the pair.  The phenotypes (PHENO) are coded in the sixth column, with the control coded with 1 (PHENO = 1) and the case with a 2 (PHENO = 2). 

Use UNIX/Linux **sort** to unscramble the PED file if necessary. The command below takes an unsorted PED file of matched case and control animals and sorts first by number in column 1 (pair identifier) and the by column 6 (control ID preceding the case ID) consistent with the requirements stated in the previous paragraph.

`sort -k1,1n -k6,6n  BCHF102pairsHD770Filtered.ped > BCHF102pairsHD770FilteredSorted.ped`

To make it easier to check that the sorting worked, use Unix/Linux **awk** to extract the first 6 columns of the PED to create a FAM file of animal metadata.

`awk '{print $1,$2,$3,$4,$5,$6}'  BCHF102pairsHD770FilteredSorted.ped >  BCHF102pairsHD770FilteredSorted.fam`

The FAM file is included in this repository.

For example, the first 5 animal pairs in correct sorted order are

	1    NE01_610_65613    0    0    2    1
	1    NE01_610_65547    0    0    2    2
	2    NE01_6DRS_18287    0    0    1    1
	2    NE01_6DRS_18803    0    0    1    2
	3    NE01_6DRS_18228    0    0    1    1
	3    NE01_6DRS_18211    0    0    1    2
	4    NE01_6DRS_18667    0    0    1    1
	4    NE01_6DRS_18688    0    0    1    2
	5    NE01_6DRH_65178    0    0    2    1
	5    NE01_6DRH_65555    0    0    2    2

 
## Input Requirements
* PED file groomed as described above
* MAP file of genome coordinates for each SNP in PED

These should be placed in the `/data` folder. Please change the filenames associated with the plinkMAP, plinkPED, and  basename in [McNemars\_Chip\_Analysis.m ](/code/McNemars_Chip_Analysis.m)file accordingly as they are currently set as defined below.

for chromosome 28 only 

	plinkMAP = 'BCHF102pairsHD770FFFSortFiltered_extract_chr28.map';
	plinkPED = 'BCHF102pairsHD770FFFSortFiltered_extract_chr28.ped';
	basename = 'BCHF102pairsHD770FFFSortFiltered_extract_chr28';
	
for entire dataset

	plinkMAP = 'BCHF102pairsHD770FFFSortFiltered.map';
	plinkPED = 'BCHF102pairsHD770FFFSortFiltered.ped';
	basename = 'BCHF102pairsHD770FFFSortFiltered';


# CodeOcean Matlab scripts
Code Ocean does not interactively run using Matlab .mlx Live Scripts. Therefor, these scripts were converted to .m non-interactive scripts to run on Code Ocean. 

## Running McNemar's Analysis 


The compute capsule is run by clicking on the **Reproducible Run** button on the upper right hand side of the Code Ocean window. The action start the **run** script in the code directory in the left side of the window.  The default run script is to run SNP on chromosome 28

### Test on SNPs from chromosome 28 only (default)

Because there is no hash sign '#' starting the lng, the script will run the McNemars\_Chip\_Analysis\_chr28.m script

	
	# Master script to run analyses on control and genomic SSR 
	# 

	matlab -nodisplay -nodesktop -r "run McNemars_Chip_Analysis_chr28.m"
	
	# matlab -nodisplay -nodesktop -r "run McNemars_Chip_Analysis_GitHub_FullChip.m"

### Run full analysis on all SNPs (edit run script)
	
It is easy to change input files used by the analysis by changing the input file names for the PED, MAP and basename for the  plinkPED, plinkMAP, and basename variables within the script. For example, a second Matlab script is available to run the analysis on all the SNP. To run the fill analysis,  edit the run script to run McNemars\_Chip\_Analysis\_GitHub\_FullChip.m to remove the '#' in front of

	# matlab -nodisplay -nodesktop -r "run McNemars_Chip_Analysis_GitHub_FullChip.m"
	
and place one in the first column position 

	matlab -nodisplay -nodesktop -r "run McNemars_Chip_Analysis_chr28.m"

resulting in a **run** file that looks like


	# Master script to run GWAS using McNemar's test 
	# 

	# matlab -nodisplay -nodesktop -r "run McNemars_Chip_Analysis_chr28.m"
	
	matlab -nodisplay -nodesktop -r "run McNemars_Chip_Analysis_GitHub_FullChip.m"
	
Running the entire dataset usually takes a few hours. The compute intensive step is calling diplotypes for every animal and SNP.

# Run Matlab Live Scripts on workstation

Run the analogous McNemars\_Chip\_Analysis\_chr28.mlx chromosome 28 only test data set or the full dataset with McNemars\_Chip\_Analysis\_GitHub\_FullChip.mlx. ***The Parallel Computing and Statistics & Machine Learning Toolboxes are required***. Download the 	BCHF102pairsHD770FFFSortFiltered.map & BCHF102pairsHD770FFFSortFiltered.ped full chip data from the Code Ocean data directory as these files are too big to store on GitHub. A pdf is provided for each Live Scripts after they were run. 

# Outputs
A CSV file of McNemar's test scores, occupancy of McNeamar's contingency table quandrants,  chi-square, chi-square continuity correction , exact p-values, & mid p-values in the  `/results` directory for Code Ocean.  MatLab Live Scripts will write the CSV into the directory the Live Script were run from as well as a binary MAT file that captures all variables generated.
 

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


### CSV output file: Explanation of column header abbreviations

HiFreqAllele = high frequency biallelic SNP allele in the group of 204 cases and controls (i.e. major allele)

For ExactlyOneAllele case with high frequency (HF) allele

* ExactOneHF_a = exactly one allele (heterozygotes, hets only) in quadrant Qa  
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


## Environment and Dependencies

Click **Environment** on the left to find a computational environment to accommodate your software (languages, frameworks) or hardware (GPU) requirements. You can then further customize the environment by installing additional packages. The changes you make will be applied the next time your capsule runs. See our help articles on [the computational environment](https://help.codeocean.com/getting-started/the-computational-environment/configuring-your-computational-environment-an-overview) for more information.

### Environment Caching

The next time you run your capsule after making changes to any part of the environment, a custom environment will be built and cached for future runs.

When you publish a capsule, its computational environment will be preserved with it, thereby ensuring computational reproducibility.

## ![](https://unlicense.org/pd-icon.png) License 

This code is released into the public domain under the [UnLicense](https://unlicense.org) 
 

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
