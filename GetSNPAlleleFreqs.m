function [AlleleHighFreq,AlleleLowFreq] = GetSNPAlleleFreqs(diploT_array_sort,SNP_Num)
%UNTITLED2 SNP_Num is the column number of diploT_array_sort
%   Detailed explanation goes here
% for each snp figure our which alleles are present and hets, figure our
% dominant allelele and put that at the top of the 2 x 2 contingency table
% each SNP diplotype in own column, need to figure out all allleles and
% most and least frequent ones

% Define Homozygous allele for diplotype

HomoAllele = char( 'A','C','G','T');

alleleFreqs = tabulate(diploT_array_sort(:,SNP_Num)); % used MATLAB baked in function for reporting

alleleFreqsSort = sortrows(alleleFreqs,3,'descend');

% may need to include the  statement below if there are 'N' no call genotypes
% if  char(alleleFreqsSortTable.alleleFreqsSort1) == 'N'


% domininant allele will be closer to the top of this list, if not the top
% and therefor will have the lowest position number

alleleFreqsSortTable = cell2table(alleleFreqsSort);

% which single alleles present in the SNP, a character vector of alleles
% present

HomoAllelePresent = char(intersect(alleleFreqsSortTable.alleleFreqsSort1,HomoAllele));

num_homo = length(HomoAllelePresent);

homo_allele_rel_pos = zeros(1,num_homo);

for i = 1:num_homo
    pos_array = char(alleleFreqsSortTable.alleleFreqsSort1) == HomoAllelePresent(i);
    homo_allele_rel_pos(i) = find(pos_array);
end

% Dominant and minor alleles

AlleleHighFreq = char(alleleFreqsSortTable.alleleFreqsSort1(min(homo_allele_rel_pos)));

if ( char(alleleFreqsSortTable.alleleFreqsSort1(max(homo_allele_rel_pos))) == 'A' || 'C' ||'G'||'T')
    AlleleLowFreq = char(alleleFreqsSortTable.alleleFreqsSort1(max(homo_allele_rel_pos)));
else % could be an N or 0 at the bottom of the list, the lowest frequency occurance
    AlleleLowFreq = char(alleleFreqsSortTable.alleleFreqsSort1(max(homo_allele_rel_pos)-1));
end


end

