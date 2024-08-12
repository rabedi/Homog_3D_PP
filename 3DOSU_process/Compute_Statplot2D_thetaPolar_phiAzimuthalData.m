function [meanDat, sdivDat, minDat, maxDat, COVDat] = Compute_Statplot2D_thetaPolar_phiAzimuthalData(...
    fileNameWOSerial, extension, polarFromPole, indexMax, indexMin, validIndices)
if (nargin < 1)
%    fileNameWOSerial = '../3DOSU_data_out/fiber_3060/fiber_3060_Cn_normal_';
%    fileNameWOSerial = '../3DOSU_data_out/fiber_z/fiber_z_sn_normal_';
%    fileNameWOSerial = '../3DOSU_data_out/fiber_z/fiber_z_Cn_normal_';
    fileNameWOSerial = '../3DOSU_data_out/fiber_3060/fiber_3060_sn_normal_';
end
if (nargin < 2)
    extension = 'mt';
end
if (nargin < 3)
polarFromPole = 1;
end
if (nargin < 4)
    indexMax = 1000;
end
if (nargin < 5)
    indexMin = 0;
end
if (nargin < 6)
    validIndices = indexMin:indexMax;
end
meanDat = [];
sdivDat = [];
minDat = [];
maxDat = [];
COVDat = [];

szz = length(validIndices);
cntr = 0;
vecOfDat = cell(0);
for I = 1:szz
   ser =  validIndices(I);
   si = num2str(ser);
   fileName = [fileNameWOSerial, si, '.', extension];
   fid = fopen(fileName, 'r');
   if (fid < 0)
       continue;
   end
   fclose(fid);
   cntr = cntr + 1;
   vecOfDat{cntr} = load(fileName);
end
[omat{1}, omat{2}, omat{3}, omat{4}, omat{5}] = Compute_mean_sdiv_min_max_Matrices(vecOfDat);
meanDat = omat{1};
sdivDat = omat{2};
minDat = omat{3};
maxDat = omat{4};
COVDat = omat{5};

added = {'mean', 'sdiv', 'min', 'max', 'COV'};
addpath('../3DCsn');
for J = 1:5
    fnnExt = [fileNameWOSerial, added{J}];
    fn = [fnnExt, '.', extension];
    dt = omat{J};
    save(fn, 'dt', '-ascii');
    [X, Y, val, successful] = plot2D_thetaPolar_phiAzimuthalData(fnnExt, extension, polarFromPole);
end

