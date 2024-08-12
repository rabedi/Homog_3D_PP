% some C
CIn = rand(6, 6);
CIn = CIn * transpose(CIn);

% some sns
num_polarThetaSegments = 30;
num_azimuthalPhiSegments = 120;
polarTheta = 70 * pi / 180;
c_p = cos(polarTheta);
s_p = sin(polarTheta);
azimuthalPhiRad = 60 * pi / 180;
c_a = cos(azimuthalPhiRad);
s_a = sin(azimuthalPhiRad);
eBase = [s_p * c_a, s_p * s_a, c_p];

num_polarThetaPoints = num_polarThetaSegments + 1;
num_azimuthalPhiPoints = num_azimuthalPhiSegments;
for i = 1:num_polarThetaPoints
    polarThetas(i) = (i - 1) / num_polarThetaSegments * pi;
end
for j = 1:num_azimuthalPhiPoints
    azimuthalPhiRads(j) = (j - 1) * (2 * pi) / num_azimuthalPhiSegments;
end
for i = 1:num_polarThetaPoints
    polarTheta = polarThetas(i);
    c_p = cos(polarTheta);
    s_p = sin(polarTheta);
    for j = 1:num_azimuthalPhiPoints
        azimuthalPhiRad = azimuthalPhiRads(j);
        c_a = cos(azimuthalPhiRad);
        s_a = sin(azimuthalPhiRad);
        e = [s_p * c_a, s_p * s_a, c_p];
        value = 2 + abs(e * eBase');
        sns(i, j) = value;
    end
end

Csn = C3D_sns;
AnisoIndexEq_max2minMinus1 = 0;
%Csn = Csn.Compute_C3D_sns(CIn, sns, AnisoIndexEq_max2minMinus1);
saveRootName = 'fiber_y';
polarFromPole = 1;
b_plot = 1;
b_print =1;
plotAddedName = 'SVE0'; 
printScalarVals = 1;
Csn = Csn.ComputePlotPrint_C3D_sns(CIn, sns, AnisoIndexEq_max2minMinus1, saveRootName, plotAddedName, printScalarVals, polarFromPole, b_plot, b_print);


a = 12;

