function ComputeWrite_mean_sdiv_min_max_Values(fileNameWOExt, extension, printText)
% opens the input file (each row a new sample) and writes it's stats
if (nargin < 1)
    fileNameWOExt = '../3DOSU_data_out/fiber_3060/fiber_3060_sn_normal';
end
if (nargin < 2)
    extension = 'sum';
end
if (nargin < 3)
    printText = 1;
end
fileName = [fileNameWOExt, '.', extension];
fid = fopen(fileName, 'r');
if (fid < 0)
    return;
end
fclose(fid);
vals = load(fileName, '-ascii');
[meanDat, sdivDat, minDat, maxDat, COVDat] = Compute_mean_sdiv_min_max_Values(vals);
fileNameOut = [fileNameWOExt, '_', extension, '.stat'];
fido = fopen(fileNameOut, 'w');
[m, n] = size(meanDat);
if (printText)
    fprintf(fido, 'mean\t');
end
for i = 1:m
    fprintf(fido, '%f\t', meanDat(i));
end
fprintf(fido, '\n');

if (printText)
    fprintf(fido, 'sdiv\t');
end
for i = 1:m
    fprintf(fido, '%f\t', sdivDat(i));
end
fprintf(fido, '\n');

if (printText)
    fprintf(fido, 'min\t');
end
for i = 1:m
    fprintf(fido, '%f\t', minDat(i));
end
fprintf(fido, '\n');

if (printText)
    fprintf(fido, 'max\t');
end
for i = 1:m
    fprintf(fido, '%f\t', maxDat(i));
end
fprintf(fido, '\n');

if (printText)
    fprintf(fido, 'COVDat\t');
end
for i = 1:m
    fprintf(fido, '%f\t', COVDat(i));
end
fprintf(fido, '\n');
fclose(fido);    