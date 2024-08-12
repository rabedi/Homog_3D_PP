function [meanDat, sdivDat, minDat, maxDat, COVDat] = Compute_mean_sdiv_min_max_Matrices(vecOfDat)
sz = length(vecOfDat);
if (sz == 0)
    meanDat = 0;
    sdivDat = 0;
    minDat = 0;
    maxDat = 0;
end
[m, n] = size(vecOfDat{1});
meanDat = zeros(m, n);
sdivDat = zeros(m, n);
minDat = inf * ones(m, n);
maxDat = -inf * ones(m, n);

sm = meanDat;
for I = 1:sz
    mat = vecOfDat{I};
    for i = 1:m
        for j = 1:n
            v = mat(i, j);
            if (v < minDat(i, j))
                minDat(i, j) = v;
            end
            if (v > maxDat(i, j))
                maxDat(i, j) = v;
            end
            sm(i, j) = sm(i, j) + v;
        end
    end
end
meanDat = 1.0 / sz * sm;

if (sz > 1)
    for I = 1:sz
        mat = vecOfDat{I};
        for i = 1:m
            for j = 1:n
                v = mat(i, j) - meanDat(i, j);
                sdivDat(i, j) = sdivDat(i, j) + v * v;
            end
        end
    end
    sdivDat = sqrt(1.0 / (sz - 1.0) * sdivDat);
end
COVDat = sdivDat ./ meanDat;
