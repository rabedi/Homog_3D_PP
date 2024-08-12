function ComputeWrite_mean_sdiv_min_max_Values_AllFiles_sameExtension(rootPlusStartFileName, extension, printText, plotPDFs, pdf_theta_phiLookFromTop)
if nargin < 1
    rootPlusStartFileName = '../3DOSU_data_out/fiber_3060/fiber_3060_';
end
if nargin < 2
    extension = 'sum';
%    extension = 'bsum';
%    extension = 'csum';
    extension = 'vcAllC6';
%    extension = 'vcAllC21';
end
if nargin < 3
    printText = 1;
end
if nargin < 4
    plotPDFs = 1;
end
if nargin < 5
    pdf_theta_phiLookFromTop = 1;
end

fidHeader = fopen(['_header_', extension, '.txt'], 'r');
if (fidHeader < 0)
    return;
end
headerText = fgetl(fidHeader);
buf = fscanf(fidHeader, '%s', 1);
numAngleFields = fscanf(fidHeader, '%d', 1);
buf = fscanf(fidHeader, '%s', 1);
for i = 1:numAngleFields
    startAngleCols(i) = fscanf(fidHeader, '%d', 1);
end
buf = fscanf(fidHeader, '%s', 1);
for i = 1:numAngleFields
    addedAngleNames{i} = fscanf(fidHeader, '%s', 1);
end
fileName4StatWOExt = [rootPlusStartFileName, extension, '.stat'];

cntr = 0;
[buf, valid] = fscanf(fidHeader, '%s', 1);
while ((valid) && (strcmp(buf, 'END') == 0))
    cntr = cntr + 1;
    fileNameWOExtMainPart{cntr} =  buf;
    [buf, valid] = fscanf(fidHeader, '%s', 1);
end
fclose(fidHeader);
sz = length(fileNameWOExtMainPart);
angleInDeg = 1;
for fi = 1:sz
    fileNameWOExt = [rootPlusStartFileName, fileNameWOExtMainPart{fi}];
    ComputeWrite_mean_sdiv_min_max_Values(fileNameWOExt, extension, printText);

    if (numAngleFields > 0)
    	ExtractPlotSphericalAngleData(fileNameWOExt, extension, fileName4StatWOExt, angleInDeg, startAngleCols, addedAngleNames, plotPDFs, pdf_theta_phiLookFromTop);
    end
end