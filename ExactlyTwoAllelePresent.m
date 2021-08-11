% membership for single allele diplotype - returns 1 if single allele
% present in diplotype

function y = ExactlyTwoAllelePresent(allele,diplotype)

% allele either A, C, G, or T
% diplotype is
%% rules:
%       A = A,A
%       C = C,C
%       G = G,G
%       T = T,T

if (allele == 'A' && diplotype == 'A' )  
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
elseif (diplotype == 'N' || diplotype == '0')
    y = 100; % look for 100 in calling function as signal to throw out pair
    return
else
    y= 0;
    return
end

end