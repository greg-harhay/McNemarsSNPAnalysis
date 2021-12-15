function y = ExactlyTwoAllelePresent(allele,diplotype)
% ExactlyTwoAllelePresent(allele,diplotype)
% Takes as input an allele (A,C,G, or T) and an animal diplotype at a 
% given SNP position and returns 1 only if the input diplotype has two 
% copies of the input allele, 0 if diplotype has a single copy of the 
% allele or no allele copies present, or 100 if there is no diplotype call
% for the animal at the given SNP position

% allele either A, C, G, or T
% diplotype is
%% rules:
%       A = A,A
%       C = C,C
%       G = G,G
%       T = T,T


if (diplotype == 'N' || diplotype == '0')
    y = 100; % look for 100 in calling function as signal to throw out pair
    return
elseif (allele == 'A' && diplotype == 'A' )  
    y = 1;
    return
elseif (allele == 'C' && diplotype == 'C' )   
    y= 1;
    return
elseif (allele == 'G' && diplotype == 'G')    
    y= 1;   
    return
elseif (allele == 'T' && diplotype == 'T')
    y= 1;
    return
else
    y= 0;
    return
end

end