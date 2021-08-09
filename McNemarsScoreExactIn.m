function [chi_sqr,chi_sqr_cc,p_exact,p_mid,Qa,Qb,Qc,Qd] = ...
    McNemarsScoreExactIn(Allele,ExactIn,diploT_array_sort,SNP_Num)
%% McNemarsScoreExactIn
% Score McNemars for matched in three different cases. Where animals in
% pairs have exactly one matching allele (hets). This is the ExactlyOne case
% has has ExactIn = 1. In the case where there can be one or two matching
% allles in each animal, the case is ExactIn = 1.5. Where there must be two
% matching alleles to the input Allele, ExactIn = 2 and is referered to
% ExactlyTwo.

% need input diplotype array pre-sorted by pair number (first) and
% controls in row i, cases in i+1

% get number of rows in diploT_array_sort


diploT_size = size(diploT_array_sort);

% case control iterator will incremnt by 2 each time as we are dealing with
% 2 rows for each pair value

% diploT_size(1) gives total number of rows.

idx = 1:2:diploT_size(1)-1;

% gives idx = 1,3,5,7 ....

% Get diplotypes for each animal in pair

% Diplotypes arrays
%   control_no_effect
%   case_affected

for i = 1:numel(idx)
    j = idx(i);
    control_no_effect(i,1) = diploT_array_sort(j,SNP_Num);
    case_affected(i,1) = diploT_array_sort(j+1,SNP_Num);
end

%%%% McNemars scheme for evaluation for ExactOneOrTwo and ExactTwo
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
% For ExactlyOneorTwo or Exactly Two alleles for each animal

% Scoring scheme for exactly one allele (ExactlyOne) - ExactIn = 1
%  is different as one has to keep track of diplotypes
% not containing the query allele or haves homozygotes of non matching
% allele in one or both animals. Discussed more carefully in
% Next need to see  allele is present in each case
% and control diplotype  -
%   Qa if in both case and control, Qa counter increments by 1
%   Qb if case has allele, and not control then Qb counter increments by 1
%   Qc if case does not have allele, but control does, then Qc counter increments by 1
%   Qd if neither control or case have high frequency allele, then Qd
%      increments by 1
%

%test = cat(2,control_no_effect,case_affected);

Qa = 0; Qb = 0; Qc=0; Qd = 0;

if ExactIn == 1  % animals have single allele only
    
    for i = 1:numel(control_no_effect)
        
        control_test = ExactlyOneAllelePresent(Allele,control_no_effect(i));
        case_test =  ExactlyOneAllelePresent(Allele,case_affected(i));
        
        if (case_test == 100 || control_test == 100) % N or 0, no score for animal(s)
            continue  % don't score pair when there is no score for the alleles for one or both animals
        end
        
        % if looking for A genotype and get animal with homozygous
        % G, then the G is scored in Qb, Qc, or Qd. Figure out first if Allele has
        % matching diplotype A genotype with matching A diplotype. Skip
        % past this animal as it has two copies.
        % It het diplotype with genotype A, need to figure out other
        % diplotype to know where to score in. Two hets go in Qa.
        % But if have A allele, and a het with A, then the other animal
        %  G diplotype to score in Qb or Qc. If two G diploypes, then score
        %  in Qd
        
        % Compare genotype to diplotype below -
        
        % homozygotes matching allele
        
        %         if Allele == control_no_effect(i,1) || Allele == case_affected(i,1)
        %
        %             continue; % have two copies of allele in case or control or both
        %             % throw out animal from scoring, jump to next animal
        %             % in for loop
        %         end
        
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
end

if ExactIn == 2
    for i = 1:numel(control_no_effect)
        control_test = ExactlyTwoAllelePresent(Allele,control_no_effect(i));
        case_test = ExactlyTwoAllelePresent(Allele,case_affected(i));
        if (case_test == 100 || control_test == 100) % N or 0, no score for animal(s)
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

