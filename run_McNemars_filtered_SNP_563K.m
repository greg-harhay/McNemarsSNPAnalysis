

% Assume starting with unsorted diploT array, need to sort using metadata
% in Paired_metadata_array

% Whole Array Operations
% First, read in paired_metadata array as Pairedmetadataarray

[DiploT_sort, Paired_metadata_sort] = SortDiploT(diploT_array,Pairedmetadataarray);

numsnp = size(DiploT_sort);

% Commands inside a loop to process every SNP in diplotype array of SNP

tic

parfor i = 1:numsnp(2)

[AlleleHighFreq(i),AlleleLowFreq(i)] = GetSNPAlleleFreqs(DiploT_sort,i);

[singleA_ChiSQR(i),singleA_ChiSQR_CC(i),singleA_p_Exact(i),singleA_a(i),singleA_b(i),...
    singleA_c(i),singleA_d(i)] = McNemarsScoreSNPSingleA(AlleleHighFreq(i), DiploT_sort,i);

[doubleA_ChiSQR(i),doubleA_ChiSQR_CC(i),doubleA_p_Exact(i),doubleA_a(i),doubleA_b(i),...
   doubleA_c(i),doubleA_d(i)] = McNemarsScoreSNPDoubleA(AlleleHighFreq(i), DiploT_sort,i);

end

toc

McNeamersTable = table(singleA_a',singleA_b',singleA_c',singleA_d',... 
   singleA_ChiSQR',singleA_ChiSQR_CC',singleA_p_Exact', ...
   doubleA_a',doubleA_b',doubleA_c',doubleA_d',... 
   doubleA_ChiSQR',doubleA_ChiSQR_CC',doubleA_p_Exact', ...
   'VariableNames',{'singleA_a','singleA_b','singleA_c',...
   'singleA_d','singleA_ChiSquared','singleA_ChiSquared_CC','singleA_p_exact',...
   'doubleA_a','doubleA_b','doubleA_c',...
   'doubleA_d','doubleA_ChiSquared','doubleA_ChiSquared_CC','doubleA_p_exact'});

McNemarsWithSNPCoordinates = horzcat(BCHF102pairsRecode2Filtered,McNeamersTable);

writetable(McNemarsWithSNPCoordinates,'McNemarsWithSNPCoordintes.csv')
