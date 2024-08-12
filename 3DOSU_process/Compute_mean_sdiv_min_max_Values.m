function [meanDat, sdivDat, minDat, maxDat, COVDat] = Compute_mean_sdiv_min_max_Values(vals)
% the difference between 
%1) Compute_mean_sdiv_min_max_Matrices(vecOfDat)
% and 
%2) this
% in 1, vecOfDat{i} is a matrix for sample i
% in 2, vals(i,:) is sample i and columns are dataTypes
[sz, m] = size(vals);
if (sz == 0)
    meanDat = 0;
    sdivDat = 0;
    minDat = 0;
    maxDat = 0;
end

meanDat = zeros(m, 1);
sdivDat = zeros(m, 1);
minDat = inf * ones(m, 1);
maxDat = -inf * ones(m, 1);

sm = meanDat;
for I = 1:sz
    mat = vals(I,:);
    for i = 1:m
        v = mat(i);
        if (v < minDat(i))
            minDat(i) = v;
        end
        if (v > maxDat(i))
            maxDat(i) = v;
        end
        sm(i) = sm(i) + v;
    end
end
meanDat = 1.0 / sz * sm;

if (sz > 1)
    for I = 1:sz
        mat = vals(I,:);
        for i = 1:m
            v = mat(i) - meanDat(i);
            sdivDat(i) = sdivDat(i) + v * v;
        end
    end
    sdivDat = sqrt(1.0 / (sz - 1.0) * sdivDat);
end
COVDat = sdivDat ./ meanDat;
