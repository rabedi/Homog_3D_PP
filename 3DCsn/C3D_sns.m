classdef C3D_sns
% values inside normal are given, in-plane, etc, are computed,    
    properties
        C = C3D;
        sns = scalarDataAt3DOrientation_normal_inPlane;
        Cns = scalarDataAt3DOrientation_normal_inPlane;
        ACns = scalarDataAt3DOrientation;
        
        % ACn stuff
        ACn_polarThetaIndex;
        ACn_azimuthalPhiIndex;
        ACn_polarThetaRad;
        ACn_azimuthalPhiRad;
        ACn_polarThetaDeg;
        ACn_azimuthalPhiDeg;
        ACn_rotatedC;
        ACn_rotatedCIso;
        ACn_outOfPlaneStiffness;
        ACn_inPlaneStiffness;
        ACn_outOfPlane2InplaneNormalStiffness;
        ACn_ACn;
    end
    methods
        function objout = ComputePlotPrint_C3D_sns(obj, C3DIn, snValsIn, AnisoIndexEq_max2minMinus1, saveRootName, plotAddedName, printScalarVals, polarFromPole, b_plot, b_print)
            objout = obj.Compute_C3D_sns(C3DIn, snValsIn, AnisoIndexEq_max2minMinus1);
            if (b_plot)
                objout.plot2DSaveCsn(polarFromPole, saveRootName, plotAddedName);
            end
            if (b_print)
                objout.printStatCsn(saveRootName, plotAddedName, printScalarVals);
            end
        end
        % AnisoIndexEq_max2minMinus1 == 1, anisotropy = max/min - 1,
        % == 0, anisotropy = sdiv / mean (both for in-plane values)
        % snVals indexed as (polar_i, azimuthal_j)
        function objout = Compute_C3D_sns(obj, C3DIn, snValsIn, AnisoIndexEq_max2minMinus1)
            objout = obj;
            objout.C = objout.C.SetC(C3DIn, 0);
            computeInPlaneVals = 1;
            objout.sns = objout.sns.SetValuesComputeAllFromNormalValues(snValsIn, AnisoIndexEq_max2minMinus1, computeInPlaneVals);
            [num_polarThetaPoints, num_azimuthalPhiPoints] = size(snValsIn);
            CnVals = zeros(num_polarThetaPoints, num_azimuthalPhiPoints);
            ACnVals = zeros(num_polarThetaPoints, num_azimuthalPhiPoints);
            polarThetas = zeros(num_polarThetaPoints, 1);
            azimuthalPhis = zeros(num_azimuthalPhiPoints, 1);
            delPolar = 0.5 * pi / (num_polarThetaPoints - 1);
            delAzimuthal = 2.0 * pi / num_azimuthalPhiPoints;
            for i = 1:num_polarThetaPoints
                polarThetas(i) = (i - 1) * delPolar;
            end
            for j = 1:num_azimuthalPhiPoints
                azimuthalPhis(j) = (j - 1) * delAzimuthal;
            end
            b_transverIsotropicNormalizedDiff = 1;
            minACn = inf;
            i_minACn = -1; j_minACn = -1;
            for i = 1:num_polarThetaPoints
                polarTheta = polarThetas(i);
                for j = 1:num_azimuthalPhiPoints
                    azimuthalPhi = azimuthalPhis(j);

                    [rotatedC, rotatedCIso, outOfPlaneStiffness, ...
                        outOfPlane2InplaneNormalStiffness, ...
                        transverIsotropicNormalizedDiff, ...
                        r3Dout, diffNorm, CNorm] = ...
                        objout.C.ComputeCRotated_CRotatedIso_etc(polarTheta, azimuthalPhi, 1, b_transverIsotropicNormalizedDiff);
                    CnVals(i, j) = outOfPlaneStiffness;
                    ACnVals(i, j) = transverIsotropicNormalizedDiff;
                    if (transverIsotropicNormalizedDiff < minACn)
                        minACn = transverIsotropicNormalizedDiff;
                        i_minACn = i;
                        j_minACn = j;
                    end
                end
            end
            i = i_minACn;
            j = j_minACn;
            polarTheta = polarThetas(i);
            azimuthalPhi = azimuthalPhis(j);
            objout.ACn_polarThetaIndex = i;
            objout.ACn_azimuthalPhiIndex = j;
            objout.ACn_polarThetaRad = polarTheta;
            objout.ACn_azimuthalPhiRad = azimuthalPhi;
            objout.ACn_polarThetaDeg = 180 / pi * polarTheta;
            objout.ACn_azimuthalPhiDeg = 180 / pi * azimuthalPhi;
            [objout.ACn_rotatedC, objout.ACn_rotatedCIso, objout.ACn_outOfPlaneStiffness, ...
                        objout.ACn_outOfPlane2InplaneNormalStiffness, ...
                        objout.ACn_ACn, ...
                        r3Dout, diffNorm, CNorm] = ...
                        objout.C.ComputeCRotated_CRotatedIso_etc(polarTheta, azimuthalPhi, 1, b_transverIsotropicNormalizedDiff);
            objout.ACn_inPlaneStiffness = objout.ACn_rotatedCIso.C(1, 1);

            objout.Cns = objout.sns.SetValuesComputeAllFromNormalValues(CnVals, AnisoIndexEq_max2minMinus1, computeInPlaneVals);
            objout.ACns = objout.ACns.Set_ScalarValuesCompute(ACnVals, 0);
        end            
        function plot2DSaveCsn(obj, polarFromPole, saveRootName, plotAddedName)
            obj.sns.plot2DSaveAll(polarFromPole, [saveRootName, 'sn'], plotAddedName);
            obj.Cns.plot2DSaveAll(polarFromPole, [saveRootName, 'Cn'], plotAddedName);
            obj.ACns.plot2DSave(polarFromPole, [saveRootName, 'ACn'], plotAddedName);
        end
        function printStatCsn(obj, baseName, plotAddedName, printScalarVals)
            obj.sns.printStat_scalarDataAt3DOrientation_normal_inPlane([baseName, 'sn'], plotAddedName, printScalarVals);
            obj.Cns.printStat_scalarDataAt3DOrientation_normal_inPlane([baseName, 'Cn'], plotAddedName, printScalarVals);
            fileName = [baseName, 'ACn.sum'];
            fid = fopen(fileName, 'a');
            fid_plotAddedName = -1;
            if (printScalarVals)
                fileName = [baseName, 'ACn', plotAddedName, '.mt'];
                fid_plotAddedName = fopen(fileName, 'w');
            end
            obj.ACns.printStat(fid, fid_plotAddedName);
            fclose(fid);
            if (fid_plotAddedName > 0)
                fclose(fid_plotAddedName);
            end
            fileName = [baseName, 'ACn_summary.bsum'];
            fid = fopen(fileName, 'a');

            fprintf(fid, '%f\t', obj.ACn_ACn);
            fprintf(fid, '%f\t', obj.ACn_outOfPlaneStiffness);
            fprintf(fid, '%f\t', obj.ACn_inPlaneStiffness);
            fprintf(fid, '%f\t', obj.ACn_outOfPlane2InplaneNormalStiffness);
            fprintf(fid, '%f\t', obj.ACn_polarThetaDeg);
            fprintf(fid, '%f\t', obj.ACn_azimuthalPhiDeg);
            fprintf(fid, '%f\t', obj.ACn_polarThetaRad);
            fprintf(fid, '%f\t', obj.ACn_azimuthalPhiRad);
            fprintf(fid, '%f\t', obj.ACn_polarThetaIndex);
            fprintf(fid, '%f\n', obj.ACn_azimuthalPhiIndex);
            fclose(fid);
            
%            fileName = [baseName, 'ACn_RotatedC', plotAddedName, '.vc'];
%            fid = fopen(fileName, 'w');
            fileName = [baseName, 'ACn_RotatedC', '.vcallC21'];
            fid = fopen(fileName, 'a');
            for i = 1:6
                for j = i:6
                    fprintf(fid, '%f\t', obj.ACn_rotatedC.C(i, j));
                end
            end
            fprintf(fid, '\n');
            fclose(fid);
            
            C11 = obj.ACn_rotatedCIso.C(1, 1);
            C12 = obj.ACn_rotatedCIso.C(1, 2);
            C66 = obj.ACn_rotatedCIso.C(6, 6);
            C13 = obj.ACn_rotatedCIso.C(1, 3);
            C33 = obj.ACn_rotatedCIso.C(3, 3);
            C44 = obj.ACn_rotatedCIso.C(4, 4);
%            fileName = [baseName, 'ACn_Rotated_CIso', plotAddedName, '.txt'];
%            fid = fopen(fileName, 'w');
            fileName = [baseName, 'ACn_Rotated_CIso', '.vcallC6'];
            fid = fopen(fileName, 'a');
            fprintf(fid, '%f\t', C11);
            fprintf(fid, '%f\t', C12);
            fprintf(fid, '%f\t', C66);
            fprintf(fid, '%f\t', C13);
            fprintf(fid, '%f\t', C33);
            fprintf(fid, '%f\n', C44);
            fclose(fid);
        end
    end
end
