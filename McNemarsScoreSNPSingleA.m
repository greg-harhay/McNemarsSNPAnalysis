function [chi_sqr,chi_sqr_cc,p_exact,Qa,Qb,Qc,Qd] = McNemarsScoreSNPSingleA(AlleleHighFreq, ...
    diploT_array_sort,SNP_Num)
%% UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Score McNemars
% need input diplotype array pre-sorted by pair number (first) and
% phenotype control roll i, case in i+1
% controls in row i, cases in i+1

% get number of rows in diploT_array_sort

diploT_size = size(diploT_array_sort);

% case control iterator will incremnt by 2 each time as we are dealing with
% 2 rows for each pair value

% diploT_size(1) gives total number of rows.

idx = 1:2:diploT_size(1)-1;

% gives idx = 1,3,5,7 ....

% outside loop will be for each column or SNP

% parfor l = 1:diploT_size(2)

for i = 1:numel(idx)
    j = idx(i);
    control_no_effect(i,1) = diploT_array_sort(j,SNP_Num);
    case_affected(i,1) = diploT_array_sort(j+1,SNP_Num);
end

%%%% McNemars scheme for evaluation
% if case has allele, and control doesn't, pair score in Qb

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
% Define that high frequency allele to located in Qa if present in both
% control and case. This guarantees that Qb/Qc is > 1
%
% Next need to see of high frequency single allele is present in each case
% and control diplotype  -
%   Qa if in both case and control, Qa counter increments by 1
%   Qb if case has allele, and not control then Qb counter increments by 1
%   Qc if case does not have allele, but control does, then Qc counter increments by 1
%   Qd if neither control or case have high frequency allele, then Qd
%      increments by 1
%

%test = cat(2,control_no_effect,case_affected);

Qa = 0; Qb = 0; Qc=0; Qd = 0;


for i = 1:numel(control_no_effect)
    control_test = singleAllelePresent(AlleleHighFreq,control_no_effect(i));
    case_test = singleAllelePresent(AlleleHighFreq,case_affected(i));
    if (case_test == 100 || control_test == 100) %  0, no score for animal(s)
        continue  % don't score pair when there is no score for the alleles for one or both animals
    elseif ( control_test == 1 && case_test == 1)
        Qa = Qa +1;
    elseif ( case_test == 1 && control_test == 0)
        Qb = Qb +1;
    elseif ( case_test == 0 && control_test == 1)
        Qc = Qc +1;
    elseif ( case_test == 0 && control_test == 0)
        Qd = Qd +1;
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
else
    p_exact = 2*binocdf(Qb,n,0.5);
end

end

