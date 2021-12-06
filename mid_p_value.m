function [mid_p_value] = mid_p_value(Qb,Qc,prob)

% 

n = Qb + Qc;

if Qc <= Qb
    mid_p_value = 2*binocdf(Qc,n,prob) - binopdf(Qc,n,prob);
else
    mid_p_value = 2*binocdf(Qb,n,prob) - binopdf(Qb,n,prob);
end

end