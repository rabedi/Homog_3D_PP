polarThetaDeg = 30;
azimuthalPhiDeg = 65;
alphaDeg = 10;
r3D = Rotation3D;
r3D = r3D.SetAnglesDeg(polarThetaDeg, azimuthalPhiDeg, alphaDeg);
Q = r3D.Q;
QQT = Q * transpose(Q);
detQ = det(Q);

%Ctest = rand(6, 6) .* rand(6, 6) + rand(6, 6) + 1;
Ctest = rand(6, 6);
Ctest = Ctest * transpose(Ctest);
c3 = C3D;
c3 = c3.SetC(Ctest, 0);


% C voigt and ths tests
Cdir = C3D;
Cdir = Cdir.SetC(Ctest, 1);
Cduplicate = Cdir.SetC4th2CVoigt();
delCs = Cdir.C - Cduplicate.C;
delCb = Cduplicate.C - Ctest;

% rotations 
isRad = 0;
[cRotatedDirect, rot3DDirect]  = c3.RotateDirectVoigt_wAngles(polarThetaDeg, azimuthalPhiDeg, alphaDeg, isRad);
[cRotated4th, rot3D4th]  = c3.RotateUsing4thOrderC_wAngles(polarThetaDeg, azimuthalPhiDeg, alphaDeg, isRad);
delCrotated = cRotatedDirect.C - cRotated4th.C


L2normFromVoigt = c3.getL2NormFromVoigt()
L2normFrom4th = c3.getL2NormFrom4th()
L2normFromVoigtRotated = cRotatedDirect.getL2NormFromVoigt()
L2normFrom4thRotated = cRotatedDirect.getL2NormFrom4th()

b_transverIsotropicNormalizedDiff = 1;
[rotatedC, rotatedCIso, outOfPlaneStiffness, outOfPlane2InplaneNormalStiffness, transverIsotropicNormalizedDiff, r3Dout, diffNorm, CNorm] = c3.ComputeCRotated_CRotatedIso_etc(polarThetaDeg, azimuthalPhiDeg, isRad, b_transverIsotropicNormalizedDiff)
outOfPlane2InplaneNormalStiffness
