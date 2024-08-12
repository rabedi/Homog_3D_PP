sd3d = scalarDataAt3DOrientation;
num_polarThetaSegmentsIn = 30;
num_azimuthalPhiSegmentsIn = 120;
sd3d = sd3d.Size_ScalarDataAt3DOrientation(num_polarThetaSegmentsIn, num_azimuthalPhiSegmentsIn);

polarTheta = 70 * pi / 180;
c_p = cos(polarTheta);
s_p = sin(polarTheta);
azimuthalPhiRad = 60 * pi / 180;
c_a = cos(azimuthalPhiRad);
s_a = sin(azimuthalPhiRad);
eBase = [s_p * c_a, s_p * s_a, c_p];

for i = 1:sd3d.num_polarThetaPoints
    polarTheta = sd3d.polarThetasRad(i);
    c_p = cos(polarTheta);
    s_p = sin(polarTheta);
    for j = 1:sd3d.num_azimuthalPhiPoints
        azimuthalPhiRad = sd3d.azimuthalPhisRad(j);
        c_a = cos(azimuthalPhiRad);
        s_a = sin(azimuthalPhiRad);
        e = [s_p * c_a, s_p * s_a, c_p];
        value = 2 + abs(e * eBase');
        sd3d.scalarDat(i, j) = value;
    end
end

sn = scalarDataAt3DOrientation_normal_inPlane;
sn.scalars_normal = sn.scalars_normal.Set_ScalarValues(sd3d.scalarDat);
computeInPlaneVals = 1;
AnisoIndexEq_max2minMinus1 = 0;
sn = sn.ComputeAllFromNormalValues(AnisoIndexEq_max2minMinus1, computeInPlaneVals);
polarFromPole = 1;
saveRootName = 'SVE1_sn';
sn.plot2DSaveAll(polarFromPole, saveRootName, 'SVE1');
baseName = 'SVE_sn';
printScalarVals = 1;
sn.printStat_scalarDataAt3DOrientation_normal_inPlane(baseName, 'SVE1', printScalarVals);
a = 12;
