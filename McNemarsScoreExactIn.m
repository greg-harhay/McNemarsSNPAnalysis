function [chi_sqr,chi_sqr_cc,p_exact,p_mid,Qa,Qb,Qc,Qd] = ...
    McNemarsScoreExactIn(Allele,ExactIn,diploT_array_sort,SNP_Num)
%
% [chi_sqr,chi_sqr_cc,p_exact,p_mid,Qa,Qb,Qc,Qd] =  ...
% McNemarsScoreExactIn(Allele,ExactIn,diploT_array_sort,SNP_Num)
%
% Inputs are allele (A,C,G, or T), ExactIn (a flag = 1, 1.5, for 2) 
% for the three different cases considered, diployT_array_sort is an array 
% of sorted diplotypes where the animals have been presorted by animal pair
% first, then case-control with control (phenotype = 1) animal always 
% preceeding the case (phenotype = 2) having been provided in the input  
% PED, and finally SNP_Num, the SNP number from PED. The outputs are 
% statistics associated with the contingencey table chi-square, chi-square 
% with continuity correction, p_exact (exact p value, from binomial 
% distro), p_mid (mid p value from binomial distro) and
% occupancy of the four coordinates Qa,Qb,Qc,Qd of the contingency table 
%
% Score McNemars Test for three different cases. 
%
% ExactlyOneAllelePresent: ExactIn = 1 
% Where an animal in the pair has exactly one matching allele (hets) in
% its diplotype.
%
% ExactlyOneorTwoAllelePresent: ExactIn = 1.5 
% Where an animal in the pair has one or two matching alleles in 
% its diplotype. 
%
% ExactlyTwoAllelePresent: ExactIn = 2  
% Where an animal in the pair has exactly two matching alleles (homos) in
% its diplotype.
%
% McNemars scheme for evaluation contingency table evaluation
%
%
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
%   Qa if allele in both case and control, Qa counter increments by 1
%   Qb if case has allele, and not control then Qb counter increments by 1
%   Qc if case does not have allele, but control does, then Qc counter increments by 1
%   Qd if neither control or case have allele, then Qd increments by 1
%


% get number of rows in diploT_array_sort

diploT_size = size(diploT_array_sort);

% Need input diplotype array pre-sorted by pair number (first) and controls in 
% row i, cases in i+1.
% Case control iterator will increment by 2 each time as we are dealing with
% 2 rows for each pair value

% diploT_size(1) gives total number of rows.

idx = 1:2:diploT_size(1)-1;

% gives idx = 1,3,5,7 ....

% Get diplotypes for each animal in pair

% Diplotypes arrays
% Control preceeds Case in sorted diplotype array, the array
% was sorted as a consequence of being derived from a sorted PED file
% provide by the user.

for i = 1:numel(idx)
    j = idx(i);
    control_no_effect(i,1) = diploT_array_sort(j,SNP_Num);
    case_affected(i,1) = diploT_array_sort(j+1,SNP_Num);
end


%test = cat(2,control_no_effect,case_affected);

Qa = 0; Qb = 0; Qc=0; Qd = 0;

if ExactIn == 1  % animals have single allele only

    for i = 1:numel(control_no_effect)

        control_test = ExactlyOneAllelePresent(Allele,control_no_effect(i));
        case_test =  ExactlyOneAllelePresent(Allele,case_affected(i));

        if (case_test == 100 || control_test == 100) % N or 0, no score for animal(s)
            continue  % don't score pair when there is no score for the alleles for one or both animals
        end

        % look or single copy of allele in

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
end

if ExactIn == 1.5
    for i = 1:numel(control_no_effect)
        control_test = ExactlyOneOrTwoAllelePresent(Allele,control_no_effect(i));
        case_test = ExactlyOneOrTwoAllelePresent(Allele,case_affected(i));
        if (case_test == 100 || control_test == 100) % N or 0, no score for animal(s)
            continue  % don't score pair when there is no score for the alleles for one or both animals
        end
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
end

if ExactIn == 2
    for i = 1:numel(control_no_effect)
        control_test = ExactlyTwoAllelePresent(Allele,control_no_effect(i));
        case_test = ExactlyTwoAllelePresent(Allele,case_affected(i));
        if (case_test == 100 || control_test == 100) % N or 0, no score for animal(s)
            continue  % don't score pair when there is no score for the alleles for one or both animals
        end
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
    p_mid = 2*binocdf(Qc,n,0.5) - binopdf(Qc,n,0.5);
else
    p_exact = 2*binocdf(Qb,n,0.5);
    p_mid = 2*binocdf(Qb,n,0.5) - binopdf(Qb,n,0.5);
end

end

