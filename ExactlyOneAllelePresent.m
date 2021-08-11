% membership for single allele diplotype - returns 1 if single allele
% present in diplotype

function y = ExactlyOneAllelePresent(allele,diplotype)

% allele either A, C, G, or T
% diplotype is
%       R = A,G or G,A
%       Y = C,T or T,C
%       S = C,G or G,C
%       W = A,T or T,A
%       K = G,T or T,G
%       M = A,C or C,A

if (allele =='A' && ismember(diplotype, ['R','W','M']))
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
elseif (diplotype == 'N' || diplotype == '0')
    y = 100; % look for 100 in calling function as signal to throw out pair
    return
else
    y= 0;
    return
end

end