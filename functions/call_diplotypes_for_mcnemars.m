function [diploT_array] = call_diplotypes_for_mcnemars(PED)
% [diploT_array] = call_diplotypes_for_mcnemars(PED)
% Input is a PED file where it is assumed that columns 1 - 6 contain sample
% (animal) metadata including animal ID and phenotype. However, for
% downstream computational efficiency, the PED rows should first be sorted 
% by animal pair, with the control animal preceeding the case animal from
% top to bottom. This function only operates on data starting in column 7 
% to the last column in the PED. This function assumes biallelic alleles 
% occuring in pairs on a given row (animal or sample). This function 
% returns an array of diplotypes for every row (animal, sample) in the PED.  

%% start clock with tic
tic

fid = fopen(PED);
line_ex = fgetl(fid);

% There  are 6 columns with metadata at the begining
% Some number, animal id, ,some number ,some number, sex, and pair id?
% diplotype pairs start at column 7

% convert diplotypes to single character representation 
%
%       C = C,C
%       G = G,G
% rules:
%       A = A,A
%       T = T,T
%       R = A,G or G,A
%       Y = C,T or T,C
%       S = C,G or G,C
%       W = A,T or T,A
%       K = G,T or T,G
%       M = A,C or C,A

% run through the pairs, one pair at a time

col_elements = strsplit(line_ex);
total_num_cols = length(col_elements);
genotype_cols(1,:) = col_elements(1,7:total_num_cols);
total_num_geno = length(genotype_cols);

diploT_pair  = char();
diploT_array= char(); 
metadata_array = char();

col_idx = int32(1:2:total_num_geno-1);

% a container is an associative arrray or hash

keySet = {'AA','CC','GG','TT','AG','GA','CT','TC','CG','GC','AT','TA','GT','TG','AC','CA', ...
    '00','0A','A0','0C','C0','0T','T0','0G','G0'};
valueSet = {'A','C','G','T','R','R','Y','Y','S','S','W','W','K','K','M','M', ...
    '0','0','0','0','0','0','0','0','0'};
diplo = containers.Map(keySet,valueSet);

while ischar (line_ex)
    
    % need to do this for every new line
    col_elements = strsplit(line_ex);
    metadata(1,:) =  col_elements(1,1:6);
    metadata_array = cat(1,  metadata_array,  metadata);
    genotype_cols(1,:) = col_elements(1,7:total_num_cols);
       
    parfor i = 1:numel(col_idx)
        j = col_idx(i);
        pair = strcat(genotype_cols(j),genotype_cols(j+1));
        diploT_pair(i) = diplo(string(pair)); % use hash for lookup
    end
    
    diploT_array = cat(1, diploT_array, diploT_pair);
    line_ex = fgetl(fid);
end

toc
%% Stop clock, time to compute

