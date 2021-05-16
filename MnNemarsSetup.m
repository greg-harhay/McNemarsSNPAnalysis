% Automated Process for setting up McNemars test for each snp
% Create array of diplotypes in ascending pair order with control first and
% case second

% read in pair infor from  Paired_Metadata_array.xls

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
% dominiant one

% Define Homozygous allele for diplotype

HomoAllele = char( 'A','C','G','T');

alleleFreqs = tabulate(diploT_array_sort(:,k/2));

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


%% Score McNemars

% Probably should set this up as a function
% controls in row i, cases in i+1 

% get number of rows in diploT_array_sort

diploT_size = size(diploT_array_sort);

% case control iterator will incremnt by 2 each time as we are dealing with
% 2 rows for each pair value

idx = 1:2:diploT_size(1)-1;

% outside loop will be for each column or SNP

% parfor l = 1:diploT_size(2)

parfor i = 1:numel(idx) 
    j = idx(i);
    control_no_effect(i,1) = diploT_array_sort(j,10);
    case_affected(i,1) = diploT_array_sort(j+1,10);  
end



%                                      Control
%                             +           |         -
%                     +                   |
%                             Qa          |         Qb
%           Case   --------------------------------------------
%                                         |
%                             Qc          |         Qd
%                                         |
%                     -                   |         
%                                         |
%
%
%
%
% Define that high frequency allele to located in Qa if present in both
% control and case. This guarantees that Qb/Qc is > 1
%
% Next need to see of high frequency single allele is present in each case
% and control diplotype  - 
%   Qa if in both case and control, Qa counter increments by 1
%   Qb if case has allele, and not control then Qb counter increments by 1
%   Qc if case does not have allele, but control does, then Qc counter increments by 1
%   Qd if neither control or case have high frequency allele, then Qd
%      increments by 1
%

%test = cat(2,control_no_effect,case_affected);

Qa = 0; Qb = 0; Qc=0; Qd = 0;


parfor i = 1:numel(control_no_effect)
    control_test = singleAllelePresent(AlleleHighFreq,control_no_effect(i));
    case_test = singleAllelePresent(AlleleHighFreq,case_affected(i));
    
    
    if ( control_test == 1 && case_test == 1)
        Qa = Qa +1;
    elseif ( case_test == 1 && control_test == 0)
        Qb = Qb +1;
    elseif ( case_test == 0 && control_test == 1)
        Qc = Qc +1;
    elseif ( case_test == 0 && control_test == 0)
        Qd = Qd +1;
    end
end

% The occupancy of the quandrants have been computed 
% 
% Compute Statistics
%
%

% McNemar's chi_sqr   chi squared

chi_sqr = ((Qb -Qc)^2)/(Qb + Qc);

% McNemars chi_sqr with continuity correction chi_sqr_cc

chi_sqr_cc = (abs(Qb-Qc)-1)^2/(Qb + Qc);

% McNemar's exact p-value

n = Qb + Qc;

if Qc <= Qb  
    p_exact = 2*binocdf(Qc,n,0.5);
else
    p_exact = 2*binocdf(Qb,n,0.5);
end



