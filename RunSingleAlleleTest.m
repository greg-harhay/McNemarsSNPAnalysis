

% Assume starting with unsorted diploT array, need to sort using metadata
% in Paired_metadata_array

% Whole Array Operations
% First, read in paired_metadata array as Pairedmetadataarray

[DiploT_sort, Paired_metadata_sort] = SortDiploT(diploT_array,Pairedmetadataarray);

% Command inside a loop to process every SNP in diplotype array of SNP

numsnp = size(DiploT_sort);


tic

parfor i = 1:numsnp(2)

[AlleleHighFreq(i),AlleleLowFreq(i)] = GetSNPAlleleFreqs(DiploT_sort,i);

[ChiSQR(i),ChiSQR_CC(i),p_Exact(i),a(i),b(i),c(i),d(i)] = McNemarsScoreSNP(AlleleHighFreq(i), ...
    DiploT_sort,i);

end

toc

