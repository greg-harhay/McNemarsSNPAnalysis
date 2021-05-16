function [diploT_array_sort, paired_metadata_sort] = SortDiploT(diploT_array,Pairedmetadataarray)
%UNTITLED Summary of this function goes here
%   
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

end

