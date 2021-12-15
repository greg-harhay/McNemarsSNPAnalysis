# McNemars SNP Analysis
Matlab scripts for performing McNemars Test on SNPs of the Illumina Bovine HD770 Chip. A summary of the McNemar's test can be found [here](https://en.wikipedia.org/wiki/McNemar%27s_test).  These scripts are written in 100% MatLab are the computational engine generating the results discussed in the publication titled "Association of ARRDC3 and NFIA variants with bovine congestive heart failure in feedlot cattle" with link HERE
 
##  PLINK File Grooming
Each row represents a single animal's genotype in the PLINK file. The PLINK file must consist of genotypes from matched pairs of animals, with a single case animal matched with a single control animal and designated as a matched pair with a pair identifier.  The PED file must have the pair identifier, a number, in first column. The case and control phenotype (PHENO) are coded in column 6, with the  control coded with 1 (PHENO = 1) and the case with a 2 (PHENO = 2). For each match pair, the control animal shall precede the case animal, top to bottom, in the PED file.    

Use UNIX/Linux **sort** to unscramble PED. The command below takes an unsorted PED of matched case and control animals and sorts first by number in column 1 (pair identifier) and the by column 6 (control ID preceding the case ID)

`sort -k1,1n -k6,6n  BCHF102pairsHD770Filtered.ped > BCHF102pairsHD770FilteredSorted.ped`

To make it easier to check that the sorting worked, use Unix/Linux **awk** to extract the first 6 columns of the PED to create a FAM file of animal metadata.

`awk '{print $1,$2,$3,$4,$5,$6}'  BCHF102pairsHD770FilteredSorted.ped >  BCHF102pairsHD770FilteredSorted.fam`

The FAM file is included in this repository.

### Mike Heaton What are other PED filtering requirements ?
 
## Input Requirements
* PED file groomed as described above
* MAP file of genome coordinates for each SNP in PED

## Outputs
* CSV file of McNemar's test scores, occupancy of McNeamar's contingency table quandrants,  chi-square, chi-square continuity correction , exact p-values, & mid p-values
* Matlab Live Script (.mlx)
* pdf representation of MatLab Livescript (.mlx.pdf) 

###  Evaluation McNemar's contingency table


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


### CSV output file: Explaination of column header abbreviations

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
* ExactOneHF_ ChiSquared_CC =  (abs(Qb-Qc)-1)^2/(Qb + Qc)
* ExactOneHF _ p _exact =  exact p-value
* neg_ log10_ ExactOneHF_ p_ exact =  -log10(exact p-value)
* ExactOneHF_ p_ mid = mid p-value
* neg_ log10_ ExactOneHF_ p_ mid =  -log10(mid p-value)

The abbreviations schema is the same for the ExactlyOneOrTwoAllele (hets or homozygotes) and ExactlyTwoAllele (homozygotes) cases using both the high and low frequency alleles.


## License
 
