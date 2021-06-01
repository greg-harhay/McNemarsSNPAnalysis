function [AlleleHighFreq,AlleleLowFreq] = GetSNPAlleleFreqs(diploT_array_sort,SNP_Num)
%UNTITLED2 SNP_Num is the column number of diploT_array_sort
%   Detailed explanation goes here
% for each snp figure our which alleles are present and hets, figure our
% dominant allelele and put that at the top of the 2 x 2 contingency table
% each SNP diplotype in own column, need to figure out all allleles and
% most and least frequent ones

% Define Homozygous allele for diplotype

HomoAllele = char( 'A','C','G','T');

diploFreqs = tabulate(diploT_array_sort(:,SNP_Num)); % used MATLAB baked in function for reporting

diploFreqsSort = sortrows(diploFreqs,3,'descend');

% may need to include the  statement below if there are 'N' no call genotypes
% if  char(diploFreqsSortTable.diploFreqsSort1) == 'N'


% domininant allele will be closer to the top of this list, if not the top
% and therefor will have the lowest position number

diploFreqsSortTable = cell2table(diploFreqsSort);

% which single alleles present in the SNP, a character vector of alleles
% present

HomoAllelePresent = char(intersect(diploFreqsSortTable.diploFreqsSort1,HomoAllele));

num_homo = length(HomoAllelePresent);

homo_allele_rel_pos = zeros(1,num_homo);

for i = 1:num_homo
    pos_array = char(diploFreqsSortTable.diploFreqsSort1) == HomoAllelePresent(i);
    homo_allele_rel_pos(i) = find(pos_array);
end

% Dominant and minor alleles

AlleleHighFreq = char(diploFreqsSortTable.diploFreqsSort1(min(homo_allele_rel_pos)));


% If num_homo = 1, then need to compute minor allelle from AlleleHighFreq
% and diplotype -- assume this conditions only occurs when number of
% diplotypes at SNP is 2, so that if the het allele is at postion 1, homo is
% at 2 or that het allele is at pos 2 if homo is at 1

if num_homo == 1
    if homo_allele_rel_pos == 1
        het_diplo_pos = 2 ;
    else  % homo_allele_rel_pos must be = 2
        het_diplo_pos = 1 ;
    end
    
    if AlleleHighFreq == 'A'
        if char(table2cell(diploFreqsSortTable(het_diplo_pos,1))) == 'R'
            AlleleLowFreq = 'G';
        elseif char(table2cell(diploFreqsSortTable(het_diplo_pos,1))) == 'W'
            AlleleLowFreq = 'T';
        elseif char(table2cell(diploFreqsSortTable(het_diplo_pos,1))) == 'M'
            AlleleLowFreq = 'C';
        end
    end
    
    if AlleleHighFreq == 'C'
        
        if char(table2cell(diploFreqsSortTable(het_diplo_pos,1))) == 'Y'
            AlleleLowFreq = 'T';
        elseif char(table2cell(diploFreqsSortTable(het_diplo_pos,1))) == 'S'
            AlleleLowFreq = 'G';
        elseif char(table2cell(diploFreqsSortTable(het_diplo_pos,1))) == 'M'
            AlleleLowFreq = 'A';
        end      
    end
    
    if AlleleHighFreq == 'G'
        if char(table2cell(diploFreqsSortTable(het_diplo_pos,1))) ==  'R'
            AlleleLowFreq = 'A';
        elseif char(table2cell(diploFreqsSortTable(het_diplo_pos,1))) == 'S'
            AlleleLowFreq = 'C';
        elseif char(table2cell(diploFreqsSortTable(het_diplo_pos,1))) == 'K'
            AlleleLowFreq = 'T';
        end
    end
    
    if AlleleHighFreq == 'T'
        if char(table2cell(diploFreqsSortTable(het_diplo_pos,1))) ==  'Y'
            AlleleLowFreq = 'C';
        elseif char(table2cell(diploFreqsSortTable(het_diplo_pos,1))) == 'W'
            AlleleLowFreq = 'A';
        elseif char(table2cell(diploFreqsSortTable(het_diplo_pos,1))) == 'K'
            AlleleLowFreq = 'G';
        end
    end
 
    % else num_homo > 1, just look the lower frequency homozygous allele
    
elseif ( char(diploFreqsSortTable.diploFreqsSort1(max(homo_allele_rel_pos))) == 'A' || 'C' ||'G'||'T')
    AlleleLowFreq = char(diploFreqsSortTable.diploFreqsSort1(max(homo_allele_rel_pos)));
    % else % could be an N or 0 at the bottom of the list, the lowest frequency occurance
    %     AlleleLowFreq = char(diploFreqsSortTable.alleleFreqsSort1(max(homo_allele_rel_pos)-1));
end


end

