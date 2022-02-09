function [AlleleHighFreq,AlleleLowFreq] = GetSNPAlleleFreqs(diploT_array_sort,SNP_Num)
% 
% [AlleleHighFreq,AlleleLowFreq] =  ... 
% GetSNPAlleleFreqs(diploT_array_sort,SNP_Num)
% 
% The inputs are a diploT_array_sort (sorted diplotype array) and SNP_Num
% (the SNP number or column number from diplotype array). The outputs are
% the AlleleHighFreq (the high frequency allele) and  AlleleLowFreq (low
% frequency allele) at each SNP.

% Define Homozygous allele for diplotype

HomoAllele = char( 'A','C','G','T');

diploFreqs = tabulate(diploT_array_sort(:,SNP_Num)); % used MATLAB baked in function for reporting

diploFreqsSort = sortrows(diploFreqs,3,'descend');

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

elseif ( char(diploFreqsSortTable.diploFreqsSort1(max(homo_allele_rel_pos))) == 'A' || 'C' ||'G'||'T')
    AlleleLowFreq = char(diploFreqsSortTable.diploFreqsSort1(max(homo_allele_rel_pos)));
end


end

