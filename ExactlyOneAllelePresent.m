function y = ExactlyOneAllelePresent(allele,diplotype)
% ExactlyOneAllelePresent(allele,diplotype)
% Takes as input an allele (A,C,G, or T) and an animal diplotype at a 
% given SNP position and returns 1 if the diplotype only has a single copy 
% of the input allele, 0 if diplotype has two copies of the allele or 
% no copy, or 100 if there is no diplotype call for the animal at the 
% given SNP position.

% allele either A, C, G, or T
% diplotype is
%       R = A,G or G,A
%       Y = C,T or T,C
%       S = C,G or G,C
%       W = A,T or T,A
%       K = G,T or T,G
%       M = A,C or C,A

if (diplotype == 'N' || diplotype == '0')
    y = 100; % look for 100 in calling function as signal to throw out pair
    return
elseif (allele =='A' && ismember(diplotype, ['R','W','M']))
    y= 1;
    return
elseif (allele == 'C' && ismember(diplotype, ['Y','S','M']))
    y= 1;
    return
elseif (allele == 'G' && ismember(diplotype, ['R','S','K']))
    y= 1;
    return
elseif (allele == 'T' && ismember(diplotype, ['Y','W','K']))
    y= 1;
    return
else
    y= 0;
    return
end

end