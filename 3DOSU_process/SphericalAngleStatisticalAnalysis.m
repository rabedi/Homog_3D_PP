function [meanAngles, r, rBased_sdiv, marginal_pdf_polar, marginal_pdf_azimuthal] = SphericalAngleStatisticalAnalysis(polarAzimuthalAngles, angleInDeg, computePDFs)
if nargin < 2
    angleInDeg = 1;
end
if nargin < 3
    computePDFs = 1;
end


% meanAngles and polarAzimuthalAngles are indexed as 1 -> polar (polar
% direction -> theta = 0); 2 -> phi, azimuthal

% marginal_pdf_polar, marginal_pdf_azimuthal, 1st component is x
% coordinate, 2 y coordinate
if (angleInDeg)
        polarAzimuthalAngles = pi / 180.0 * polarAzimuthalAngles;
end
[m, col] = size(polarAzimuthalAngles);
% col should be 2
sm = zeros(3, 1);
for i = 1:m
    polarTheta = polarAzimuthalAngles(i, 1);
    azimuthalPhi = polarAzimuthalAngles(i, 2);
    ct = cos(polarTheta); st = sin(polarTheta);
    ca = cos(azimuthalPhi); sa = sin(azimuthalPhi);
    x = st * ca;
    y = st * sa;
    z = ct;
    sm(1) = sm(1) + x;
    sm(2) = sm(2) + y;
    sm(3) = sm(3) + z;
end

rVec = 1.0 / m * sm;
r = norm(rVec);
if (r < 1e-15)
    meanAngles = [0, 0];
    rBased_sdiv = inf;
else
    dir = rVec / r;
    thetaMean = acos(dir(3));
    azimuthalMean = atan2(dir(2), dir(1));
    if (thetaMean > 0.5 * pi)
        thetaMean = pi - thetaMean;
        azimuthalMean = azimuthalMean + pi;
    end
    if (azimuthalMean < 0)
        azimuthalMean = 2 * pi + azimuthalMean;
    end
    meanAngles = [thetaMean, azimuthalMean];
    if (angleInDeg)
        meanAngles = 180 / pi * meanAngles;
    end
    rBased_sdiv = -log(r);
end

if (computePDFs == 0)
    marginal_pdf_polar = cell(0);
    marginal_pdf_azimuthal = cell(0);
    return;
end

addpath('../../Circular Stats');
addpath('../../Circular Stats/CircStat2012a');


if (angleInDeg)
    polarAzimuthalAngles = 180 / pi * polarAzimuthalAngles;
end
azimuthalAngles = polarAzimuthalAngles(:, 2);
numPts = 180;
vfWeights = [];
isCDF = 0;
angleMin = 0;
isPeriodic = 0;
%option = 2;
option = 1;
[pdf_x, pdf_y] = getStat_PDFCDF_xy_angularData(azimuthalAngles, option, angleInDeg, numPts, vfWeights, isCDF, angleMin, isPeriodic);
marginal_pdf_azimuthal{1} = pdf_x;
marginal_pdf_azimuthal{2} = pdf_y;

%%%%%%%%%%%%
%lims = [-1e-6, pi/2 + 1e-6];
lims = [0, pi/2];
if (angleInDeg)
    lims = lims * 180 / pi;
end
xvls = 0:lims(2)/30:lims(2);
polarAngles = polarAzimuthalAngles(:, 1);
[pdf_y, pdf_x] = ksdensity(polarAngles, xvls);
%[pdf_y, pdf_x] = ksdensity(polarAngles);
%option = 1;
%[pdf_x, pdf_y] = getStat_PDFCDF_xy_angularData(polarAngles, option, angleInDeg, numPts, vfWeights, isCDF, angleMin, 0);
%[pdf_y, pdf_x] = ksdensity(polarAngles, 'Support', lims, 'BoundaryCorrection', 'reflection');
%[pdf_y, pdf_x] = ksdensity(polarAngles, 'Support', lims, 'BoundaryCorrection', 'reflection', 'Bandwidth', 0.001);
%[pdf_y, pdf_x] = ksdensity(polarAngles, 'Support', lims, 'Bandwidth', 1.01);
%[pdf_x, pdf_y] = getStat_PDFCDF_xy_angularData(polarAngles, option, angleInDeg, numPts, vfWeights, isCDF, angleMin, 0);
marginal_pdf_polar{1} = pdf_x;
marginal_pdf_polar{2} = pdf_y;
%plot(marginal_pdf_polar{1},marginal_pdf_polar{2},'--r','LineWidth',1.5)


