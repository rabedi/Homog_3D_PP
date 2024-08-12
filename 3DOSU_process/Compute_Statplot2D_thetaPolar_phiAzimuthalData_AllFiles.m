function Compute_Statplot2D_thetaPolar_phiAzimuthalData_AllFiles(...
rootFoldePlusFile, polarFromPole, indexMax, indexMin, validIndices)

if (nargin < 1)
    rootFoldePlusFile = '../3DOSU_data_out/fiber_3060/fiber_3060_';
end

if (nargin < 2)
polarFromPole = 1;
end
if (nargin < 3)
    indexMax = 1000;
end
if (nargin < 4)
    indexMin = 0;
end
if (nargin < 5)
    validIndices = indexMin:indexMax;
end
extension = 'mt';

fidconf = fopen('_header_mt.txt', 'r');
[buf, valid] = fscanf(fidconf, '%s', 1);
cntr = 0;
while (valid)
    cntr = cntr + 1;
    fileNameWOSerials{cntr} = [rootFoldePlusFile, buf, '_'];
    [buf, valid] = fscanf(fidconf, '%s', 1);
end
fclose(fidconf);



sz = length(fileNameWOSerials);
for i = 1:sz
%[meanDat, sdivDat, minDat, maxDat, COVDat] = 
    Compute_Statplot2D_thetaPolar_phiAzimuthalData(...
    fileNameWOSerials{i}, extension, polarFromPole, indexMax, indexMin, validIndices);
end
