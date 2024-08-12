classdef scalarDataAt3DOrientation_normal_inPlane
% values inside normal are given, in-plane, etc, are computed,    
    properties
        % scalars_normal.scalarDat is provided
        scalars_normal = scalarDataAt3DOrientation;
        inPlaneComputed = 0;

        % computed
        % averages in-plane normal to "normal" vector
        scalars_inPlane_average = scalarDataAt3DOrientation;
        % compute anisotropy index from in-plane values
        scalars_inPlane_AnisoIndex = scalarDataAt3DOrientation;
        % computes normal to in-plane average (ratio = normal /
        % average(inPLane))
        scalars_normal2InplaneAverage = scalarDataAt3DOrientation;
        % stores ratio + ratio^-1 - 2
        scalars_normal2InplaneIndex = scalarDataAt3DOrientation;
        
        
        inPlaneVals4TransIsoDir;
        normalVal4TransIsoDir;
        inPlane_average4TransIsoDir;
        inPlane_AnisoIndex4TransIsoDir;
        normal2InplaneAverage4TransIsoDir;
        normal2InplaneIndex4TransIsoDir;
        
        iPolarIndex4TransIsoDir;
        jAzimuthalIndex4TransIsoDir;
        polarAzimuthalAngles4TransIsoDirDeg;
        polarAzimuthalAngles4TransIsoDirRad;
    end
    methods
        % computeInPlaneVals == 1: computes in-plane values
        % AnisoIndexEq_max2minMinus1 == 1, anisotropy = max/min - 1,
        % == 0, anisotropy = sdiv / mean (both for in-plane values)
        
        function objout = SetValuesComputeAllFromNormalValues(obj, scalarValsIn, AnisoIndexEq_max2minMinus1, computeInPlaneVals)
            objout = obj;
            objout.scalars_normal = objout.scalars_normal.Set_ScalarValues(scalarValsIn);
        	objout = objout.ComputeAllFromNormalValues( AnisoIndexEq_max2minMinus1, computeInPlaneVals)
        end
        function objout = ComputeAllFromNormalValues(obj, AnisoIndexEq_max2minMinus1, computeInPlaneVals)
            objout = obj;
            objout.scalars_normal = objout.scalars_normal.ComputeStatistics();
            objout.inPlaneComputed = 0;
            if (~computeInPlaneVals)
                return;
            end
            objout.scalars_normal = objout.scalars_normal.ComputeInPlaneValues();
            % compute, in-plane, min, max, mean, sdiv
            num_polarThetaPoints = objout.scalars_normal.num_polarThetaPoints;
            num_azimuthalPhiPoints = objout.scalars_normal.num_azimuthalPhiPoints;
%            mins = inf * ones(num_polarThetaPoints, num_azimuthalPhiPoints);
%            maxs = -inf * ones(num_polarThetaPoints, num_azimuthalPhiPoints);
            means = zeros(num_polarThetaPoints, num_azimuthalPhiPoints);
%            sdivs = zeros(num_polarThetaPoints, num_azimuthalPhiPoints);
            aInds = zeros(num_polarThetaPoints, num_azimuthalPhiPoints);
            n2ts = zeros(num_polarThetaPoints, num_azimuthalPhiPoints);
            n2tInds = zeros(num_polarThetaPoints, num_azimuthalPhiPoints);
            
            for i = 1:num_polarThetaPoints
                for j = 1:num_azimuthalPhiPoints
                    normalV = objout.scalars_normal.scalarDat(i, j);
                    vc = objout.scalars_normal.inPlaneValues{i, j};
                    mn = min(abs(vc));
                    mx = max(abs(vc));
                    mnv = mean(vc);
                    sdv = std(vc);
                    if (AnisoIndexEq_max2minMinus1)
                        if (abs(mn) > 0)
                            aind = abs(mx)/abs(mn) - 1;
                        else
                            aind = nan;
                        end
                    else
                        if (abs(mnv) > 0)
                            aind = sdv / abs(mnv);
                        else
                            aind = nan;
                        end
                    end
                    if (mnv > 0)
                        n2t = abs(normalV / mnv);
                        n2tInd = n2t + 1/n2t - 2;
                    else
                        n2t = nan;
                        n2tInd = nan;
                    end
%                    mins(i, j) = mn;
%                    maxs(i, j) = mx;
                    means(i, j) = mnv;
%                    sdivs(i, j) = sdv;
                    aInds(i, j) = aind;
                    n2ts(i, j) = n2t;
                    n2tInds(i, j) = n2tInd;
                end
            end
            objout.scalars_inPlane_average = objout.scalars_inPlane_average.Set_ScalarValues(means);
            objout.scalars_inPlane_AnisoIndex = objout.scalars_inPlane_AnisoIndex.Set_ScalarValues(aInds);
            objout.scalars_normal2InplaneAverage = objout.scalars_normal2InplaneAverage.Set_ScalarValues(n2ts);
            objout.scalars_normal2InplaneIndex = objout.scalars_inPlane_average.Set_ScalarValues(n2tInds);

            objout.scalars_inPlane_average = objout.scalars_inPlane_average.ComputeStatistics();
            objout.scalars_inPlane_AnisoIndex = objout.scalars_inPlane_AnisoIndex.ComputeStatistics();
            objout.scalars_normal2InplaneAverage = objout.scalars_normal2InplaneAverage.ComputeStatistics();
            objout.scalars_normal2InplaneIndex = objout.scalars_normal2InplaneIndex.ComputeStatistics();
            
            objout.iPolarIndex4TransIsoDir = objout.scalars_inPlane_AnisoIndex.extermumLocIndex{1}(1);
            objout.jAzimuthalIndex4TransIsoDir = objout.scalars_inPlane_AnisoIndex.extermumLocIndex{1}(2);
            i = objout.iPolarIndex4TransIsoDir;
            j = objout.jAzimuthalIndex4TransIsoDir;
            objout.polarAzimuthalAngles4TransIsoDirDeg(1) = objout.scalars_normal.polarThetasDeg(objout.iPolarIndex4TransIsoDir);
            objout.polarAzimuthalAngles4TransIsoDirDeg(2) = objout.scalars_normal.azimuthalPhisDeg(objout.jAzimuthalIndex4TransIsoDir);
            objout.polarAzimuthalAngles4TransIsoDirRad = objout.polarAzimuthalAngles4TransIsoDirDeg * pi / 180;
            objout.normalVal4TransIsoDir = objout.scalars_normal.scalarDat(i, j);
            objout.inPlaneVals4TransIsoDir = objout.scalars_normal.inPlaneValues{i, j};
            
            objout.inPlane_average4TransIsoDir = objout.scalars_inPlane_average.scalarDat(i, j);
            objout.inPlane_AnisoIndex4TransIsoDir = objout.scalars_inPlane_AnisoIndex.scalarDat(i, j);
            objout.normal2InplaneAverage4TransIsoDir = objout.scalars_normal2InplaneAverage.scalarDat(i, j);
            objout.normal2InplaneIndex4TransIsoDir = objout.scalars_normal2InplaneIndex.scalarDat(i, j);
            
            objout.scalars_normal.inPlaneValues = cell(0, 0);
            objout.inPlaneComputed = 1;
        end
        function plot2DSaveAll(obj, polarFromPole, saveRootName, plotAddedName)
            obj.scalars_normal.plot2DSave(polarFromPole, [saveRootName, '_normal'], plotAddedName);
            if (obj.inPlaneComputed == 1)
                obj.scalars_inPlane_average.plot2DSave(polarFromPole, [saveRootName, '_inPlane_average'], plotAddedName);
                obj.scalars_inPlane_AnisoIndex.plot2DSave(polarFromPole, [saveRootName, '_inPlane_AnisoIndex'], plotAddedName);
                obj.scalars_normal2InplaneAverage.plot2DSave(polarFromPole, [saveRootName, '_normal2InplaneAverage'], plotAddedName);
                obj.scalars_normal2InplaneIndex.plot2DSave(polarFromPole, [saveRootName, '_normal2InplaneIndex'], plotAddedName);
            end
        end
        function printStat_scalarDataAt3DOrientation_normal_inPlane(obj, baseName, plotAddedName, printScalarVals)
            fileName = [baseName, '_normal.sum'];
            fid = fopen(fileName, 'a');
            if (fid < 0)
                fprintf(1, 'Cannot open file %s\n', fileName);
                pause;
            end
            fidScalar = -1;
            if (printScalarVals)
                fileName = [baseName, '_normal', plotAddedName, '.mt'];
                fidScalar = fopen(fileName, 'w');
            end
            obj.scalars_normal.printStat(fid, fidScalar);
            fclose(fid);
            if (fidScalar > 0)
                fclose(fidScalar);
            end
           if (obj.inPlaneComputed == 1)
                fileName = [baseName, '_inPlane_average.sum'];
                fid = fopen(fileName, 'a');
                if (fid < 0)
                    fprintf(1, 'Cannot open file %s\n', fileName);
                    pause;
                end
                fidScalar = -1;
                if (printScalarVals)
                    fileName = [baseName, '_inPlane_average', plotAddedName, '.mt'];
                    fidScalar = fopen(fileName, 'w');
                end
                obj.scalars_inPlane_average.printStat(fid, fidScalar);
                fclose(fid);
                if (fidScalar > 0)
                    fclose(fidScalar);
                end
               
                fileName = [baseName, '_inPlane_AnisoIndex.sum'];
                fid = fopen(fileName, 'a');
                if (fid < 0)
                    fprintf(1, 'Cannot open file %s\n', fileName);
                    pause;
                end
                fidScalar = -1;
                if (printScalarVals)
                    fileName = [baseName, '_inPlane_AnisoIndex', plotAddedName, '.mt'];
                    fidScalar = fopen(fileName, 'w');
                end
                obj.scalars_inPlane_AnisoIndex.printStat(fid, fidScalar);
                fclose(fid);
                if (fidScalar > 0)
                    fclose(fidScalar);
                end
               
                fileName = [baseName, '_normal2InplaneAverage.sum'];
                fid = fopen(fileName, 'a');
                if (fid < 0)
                    fprintf(1, 'Cannot open file %s\n', fileName);
                    pause;
                end
                fidScalar = -1;
                if (printScalarVals)
                    fileName = [baseName, '_normal2InplaneAverage', plotAddedName, '.mt'];
                    fidScalar = fopen(fileName, 'w');
                end
                obj.scalars_normal2InplaneAverage.printStat(fid, fidScalar);
                fclose(fid);
                if (fidScalar > 0)
                    fclose(fidScalar);
                end

                fileName = [baseName, '_normal2InplaneIndex.sum'];
                fid = fopen(fileName, 'a');
                if (fid < 0)
                    fprintf(1, 'Cannot open file %s\n', fileName);
                    pause;
                end
                fidScalar = -1;
                if (printScalarVals)
                    fileName = [baseName, '_normal2InplaneIndex', plotAddedName, '.mt'];
                    fidScalar = fopen(fileName, 'w');
                end
                obj.scalars_normal2InplaneIndex.printStat(fid, fidScalar);
                fclose(fid);
                if (fidScalar > 0)
                    fclose(fidScalar);
                end

                fileName = [baseName, '_inPlane_AnisoIndexSummary.csum'];
                fid = fopen(fileName, 'a');
                if (fid < 0)
                    fprintf(1, 'Cannot open file  %s\n', fileName);
                    pause;
                end
                fprintf(fid, '%f\t', obj.normalVal4TransIsoDir);
                fprintf(fid, '%f\t', obj.inPlane_average4TransIsoDir);
                fprintf(fid, '%f\t', obj.inPlane_AnisoIndex4TransIsoDir);
                fprintf(fid, '%f\t', obj.normal2InplaneAverage4TransIsoDir);
                fprintf(fid, '%f\n', obj.normal2InplaneIndex4TransIsoDir);
                fclose(fid);

                fileName = [baseName, '_inPlane_AnisoIndexInPlaneValues', plotAddedName, '.ipvc'];
                fid = fopen(fileName, 'w');
                if (fid < 0)
                    fprintf(1, 'Cannot open file %s\n', fileName);
                    pause;
                end
                for i = 1:length(obj.inPlaneVals4TransIsoDir)
                    fprintf(fid, '%f\t', obj.inPlaneVals4TransIsoDir(i));
                end
                fprintf(fid, '%\n');
                fclose(fid);
           end
        end
    end
end
