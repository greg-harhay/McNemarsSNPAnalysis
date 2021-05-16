% membership for single allele diplotype - returns 1 if single allele
% present in diplotype

function y = singleAllelePresent(allele,diplotype)

% allele either A, C, G, or T
% diplotype is
%       R = A,G or G,A
%       Y = C,T or T,C
%       S = C,G or G,C
%       W = A,T or T,A
%       K = G,T or T,G
%       M = A,C or C,A

if (allele == 'A' && diplotype == 'A' || allele == 'A' && diplotype == 'R' || ...
        allele == 'A' && diplotype == 'W'|| allele == 'A' &&  diplotype == 'M')  
    y = 1;
    return
elseif (allele == 'C' && diplotype == 'C' || allele == 'C' && diplotype == 'Y' || ...
        allele == 'C' && diplotype == 'S'|| allele == 'C' && diplotype == 'M')   
    y= 1;
    return
elseif (allele == 'G' && diplotype == 'G' || allele == 'G' && diplotype == 'R' || ... 
        allele == 'G' && diplotype == 'S'|| allele == 'G' && diplotype == 'K')    
    y= 1;   
    return
elseif (allele == 'T' && diplotype == 'T' || allele == 'T' && diplotype == 'Y' || ...
        allele == 'T' && diplotype == 'W'|| allele == 'T' && diplotype == 'K')
    y= 1;
    return 
elseif (allele == 'N' || allele == '0')
    y = 100; % look for 100 in calling function as signal to throw out pair
    return
else
    y= 0;
    return
end

end