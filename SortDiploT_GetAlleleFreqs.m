function [diploT_array_sort, paired_metadata_sort] = SortDiploT(diploT_array,Pairedmetadataarray)
%UNTITLED Summary of this function goes here
%   [
% Automated Process for setting up McNemars test for each snp
% Create array of diplotypes in ascending pair order with control first and
% case second

% read in pair information  from  Paired_Metadata_array.xls

%get unique values in list of pairs  

pair_ids = unique(Pairedmetadataarray.Pairs);
diploT_array_sort = char();
paired_metadata_sort = table();

% Sort and make new sorted array

k = 0; % new iterator

for i = 1:numel(pair_ids)
    for j = 1:2 
        m = pair_ids(i);
        k = k+1;
        new_idx = Pairedmetadataarray.Pairs == m & Pairedmetadataarray.Phenotype == j;       
        l = int32(find(new_idx));   
        diploT_array_sort(k,:) = diploT_array(l,:);
        paired_metadata_sort(k,:) = Pairedmetadataarray(l,:);
    end
end


% for each snp figure our which alleles are present and hets, figure our
% dominant allelele and put that at the top of the 2 x 2 contingency table
% each SNP diplotype in own column, need to figure out all allleles and
% most and least frequent ones

% Define Homozygous allele for diplotype

HomoAllele = char( 'A','C','G','T');

alleleFreqs = tabulate(diploT_array_sort(:,104));

alleleFreqsSort = sortrows(alleleFreqs,3,'descend'); 

% domininant allele will be closer to the top of this list, if not the top
% and therefor will have the lowest position number 

alleleFreqsSortTable = cell2table(alleleFreqsSort);

 % which single alleles present in the SNP, a character vector of alleles
 % present
 
 HomoAllelePresent = char(intersect(alleleFreqsSortTable.alleleFreqsSort1,HomoAllele));
 
 num_homo = length(HomoAllelePresent);
 
 % if num_homo = 1, skip SNP, not informative
 % figure out which allele is the most populated
 
 
 if num_homo > 1     % don't bother with completely homozygotic SNP
     
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

end

