classdef scalarDataAt3DOrientation
    properties
        % inputs
        num_polarThetaSegments;
        num_azimuthalPhiSegments;
        
        % indexed as (polari, azimuthalj)
        % should be filled out after size function
        scalarDat;
        
        % comptued
        num_polarThetaPoints;
        num_azimuthalPhiPoints; % does not repeat 0 -> 360
        
        polarThetasDeg;
        azimuthalPhisDeg;
        polarThetasRad;
        azimuthalPhisRad;
        
        c_polarThetasRad, s_polarThetasRad;
        c_azimuthalPhisRad, s_azimuthalPhisRad;
        
        del_azimuthalPhisRad;
        del_polarThetasRad;
        
        %%%%%%%%%%%%%%% data only used for computing anisotropy indices
        %%%%%%%%%%%%%%% from an original data
        % data useful for computing transverse isotropic indices
        % indexed as {azimuthalAngle, polarIndex}(alphaIndex) ...
        inPlaneValues;
        
        %%%%%%%%%%%%%%%%%%%
        % statistics
        meanV;
        sdiv;
        % 1, min, 2 max
        extremumV;
        % 1 and 2 min and max; 1 -> polar, 2 azimuthal 
        extermumLocIndex;
        extermumLocDeg;
        extermumLocRad;
    end
    methods
        function objout = Set_ScalarValuesCompute(obj, scalarValuesIn, b_computeInPlane)
            objout = obj.Set_ScalarValues(scalarValuesIn);
            if (b_computeInPlane)
                objout = objout.ComputeInPlaneValues();
            end
            objout = objout.ComputeStatistics();
        end
        function objout = Set_ScalarValues(obj, scalarValuesIn)
            [num_polarThetaPoints, num_azimuthalPhiPoints] = size(scalarValuesIn);
            num_polarThetaSegmentsIn = num_polarThetaPoints - 1;
            num_azimuthalPhiSegmentsIn = num_azimuthalPhiPoints;
            objout = obj.Size_ScalarDataAt3DOrientation(num_polarThetaSegmentsIn, num_azimuthalPhiSegmentsIn);
            objout.scalarDat = scalarValuesIn;
        end
        function objout = Size_ScalarDataAt3DOrientation(obj, num_polarThetaSegmentsIn, num_azimuthalPhiSegmentsIn)
            objout = obj;
            objout.num_polarThetaSegments = num_polarThetaSegmentsIn;
            objout.num_azimuthalPhiSegments = num_azimuthalPhiSegmentsIn;
            
            objout.num_polarThetaPoints = objout.num_polarThetaSegments + 1;
            objout.num_azimuthalPhiPoints = objout.num_azimuthalPhiSegments;
            
            objout.scalarDat = zeros(objout.num_polarThetaPoints, objout.num_azimuthalPhiPoints);
            objout.polarThetasDeg = zeros(objout.num_polarThetaPoints, 1);
            del = 90 / objout.num_polarThetaSegments;
            for i = 1:objout.num_polarThetaPoints
                objout.polarThetasDeg(i) = (i - 1) * del;
            end
            objout.azimuthalPhisDeg = zeros(objout.num_azimuthalPhiPoints, 1);
            del = 360 / objout.num_azimuthalPhiSegments;
            for i = 1:objout.num_azimuthalPhiPoints
                objout.azimuthalPhisDeg(i) = (i - 1) * del;
            end
            
            factor = pi() / 180.0;
            objout.polarThetasRad = factor * objout.polarThetasDeg;
            objout.azimuthalPhisRad = factor * objout.azimuthalPhisDeg;
            objout.del_azimuthalPhisRad = objout.azimuthalPhisRad(2);
            objout.del_polarThetasRad = objout.polarThetasRad(2);
            objout.c_polarThetasRad = cos(objout.polarThetasRad);
            objout.s_polarThetasRad = sin(objout.polarThetasRad);
            objout.c_azimuthalPhisRad = cos(objout.azimuthalPhisRad);
            objout.s_azimuthalPhisRad = sin(objout.azimuthalPhisRad);
        end
        function objout = ComputeInPlaneValues(obj)
            objout = obj;
            objout.inPlaneValues = cell(objout.num_polarThetaPoints, objout.num_azimuthalPhiPoints);
            inPlaneValsPol = zeros(objout.num_azimuthalPhiPoints, 1);
            for j = 1:objout.num_azimuthalPhiPoints
                inPlaneValsPol(j) = objout.scalarDat(objout.num_polarThetaPoints, j);
            end
            for j = 1:objout.num_azimuthalPhiPoints
                objout.inPlaneValues{1, j} = inPlaneValsPol;
            end
            for i = 2:objout.num_polarThetaPoints
                %                    polarThetaRad = objout.polarThetasRad(i);
                c_polar = objout.c_polarThetasRad(i);
                s_polar = objout.s_polarThetasRad(i);
                
                for j = 1:objout.num_azimuthalPhiPoints
                    %                        azimuthalPhiRad = objout.azimuthalPhisRad(j);
                    c_az = objout.c_azimuthalPhisRad(j);
                    s_az = objout.s_azimuthalPhisRad(j);
                    
                    inPlaneVals = zeros(objout.num_azimuthalPhiPoints, 1);
                    for k = 1:objout.num_azimuthalPhiPoints
                        c_al = objout.c_azimuthalPhisRad(k);
                        s_al = objout.s_azimuthalPhisRad(k);
                        %                            alphaRad = objout.azimuthalPhisRad(k);
                        
                        eInplane(1) = -(c_al * s_az + s_al * c_az * c_polar);
                        eInplane(2) =   c_al * c_az - s_al * s_az * c_polar;
                        eInplane(3) =   s_al * s_polar;
                        % assume values depend on direction and not +/-
                        % orientation of the vector
                        if (eInplane(3) < 0)
                            eInplane = -eInplane;
                        end
                        polar_Prime = acos(eInplane(3));
                        s_polar_Prime = sin(polar_Prime);
                        if (abs(s_polar_Prime) < 1e-5) % at the polar
                            inPlaneVals(k) = objout.scalarDat(1, 1);
                            continue;
                        end
                        azimuthal_Prime = atan2(eInplane(2), eInplane(1));
                        if (azimuthal_Prime < 0)
                            azimuthal_Prime = azimuthal_Prime + 2 * pi;
                        end
                        tol = 1e-6;
                        tolp = tol * objout.del_polarThetasRad;
                        if (polar_Prime <= tolp)
                            polar_indexLower = 1;
                            polar_indexHigher = 1;
                            weightPolarMax = 1;
                            weightPolarMin = 0;
                        elseif (polar_Prime - pi/2 >= -tolp)
                            polar_indexLower = objout.num_polarThetaPoints;
                            polar_indexHigher = objout.num_polarThetaPoints;
                            weightPolarMax = 1;
                            weightPolarMin = 0;
                        else
                            polar_indexLower = floor(polar_Prime / objout.del_polarThetasRad + 0.5 * tol);
                            weightPolarMax = polar_Prime / objout.del_polarThetasRad - polar_indexLower;
                            weightPolarMin = 1 - weightPolarMax;
                            polar_indexLower = polar_indexLower + 1;
                            polar_indexHigher =  polar_indexLower + 1;
                        end
                        tolaz = objout.del_azimuthalPhisRad * tol;
                        if ((abs(azimuthal_Prime) < tolaz) || (abs(2 * pi - azimuthal_Prime) < tolaz))
                            azimuthal_indexLower = 1;
                            azimuthal_indexHigher = 1;
                            weightAzimuthalMax = 1;
                            weightAzimuthalMin = 0;
                        else
                            azimuthal_indexLower = floor(azimuthal_Prime / objout.del_azimuthalPhisRad + 0.5 * tol);
                            weightAzimuthalMax = azimuthal_Prime / objout.del_azimuthalPhisRad - azimuthal_indexLower;
                            weightAzimuthalMin = 1 - weightAzimuthalMax;
                            azimuthal_indexLower = azimuthal_indexLower + 1;
                            azimuthal_indexHigher = azimuthal_indexLower + 1;
                            if (azimuthal_indexHigher == (objout.num_azimuthalPhiPoints + 1))
                                azimuthal_indexHigher = 1;
                            end
                        end
                        valPolarLow_azimuthalLow = objout.scalarDat(polar_indexLower, azimuthal_indexLower);
                        valPolarLow_azimuthalHigh = objout.scalarDat(polar_indexLower, azimuthal_indexHigher);
                        valPolarHigh_azimuthalLow = objout.scalarDat(polar_indexHigher, azimuthal_indexLower);
                        valPolarHigh_azimuthalHigh = objout.scalarDat(polar_indexHigher, azimuthal_indexHigher);
                        
                        vl = weightPolarMin * weightAzimuthalMin * valPolarLow_azimuthalLow + ...
                            weightPolarMin * weightAzimuthalMax * valPolarLow_azimuthalHigh + ...
                            weightPolarMax * weightAzimuthalMin * valPolarHigh_azimuthalLow + ...
                            weightPolarMax * weightAzimuthalMax * valPolarHigh_azimuthalHigh;
                            
                        inPlaneVals(k) = vl;
                    end
                    objout.inPlaneValues{i, j} = inPlaneVals;
                end
            end
        end
        function [X, Y, val] = plot2D(obj, polarFromPole)
            azimuthal = [obj.azimuthalPhisDeg', 360];
            polar = obj.polarThetasDeg;
            [X, Y] = meshgrid(azimuthal, polar);
            for i = 1:obj.num_azimuthalPhiPoints + 1
                ai = i;
                if (ai == obj.num_azimuthalPhiPoints + 1)
                    ai = 1;
                end
                for j = 1:obj.num_polarThetaPoints
                    pj = j;
                    if (~polarFromPole)
                        pj = obj.num_polarThetaPoints + 1 - j;
                    end
                    val(j, i) = obj.scalarDat(pj, ai);
                end
            end
            contourf(X,Y,val, 'EdgeColor', 'none');
            colorbar;
            fs = 18;
            hx = xlabel('$$\phi$$ (Azimuthal)', 'interpreter', 'latex', 'FontSize', fs);
            set(hx, 'VerticalAlignment','Top');
            hy = ylabel('$$\theta$$ (Polar)', 'interpreter', 'latex', 'FontSize', fs);
            set(hy, 'VerticalAlignment','Bottom');
        end
        
        function [X, Y, val] = plot2DSave(obj, polarFromPole, saveRootName, plotAddedName)
            f = figure(1000);
            [X, Y, val] = plot2D(obj, polarFromPole);
            savefig([saveRootName, plotAddedName, '.fig']);
            print('-dpng', [saveRootName, plotAddedName, '.png']);
            close(f);
        end
        function plot3D(obj)
            dat = zeros(obj.num_azimuthalPhiPoints + 1, obj.num_polarThetaPoints);
            for i = 1:obj.num_azimuthalPhiPoints + 1
                ai = i;
                if (ai == obj.num_azimuthalPhiPoints + 1)
                    ai = 1;
                end
                for j = 1:obj.num_polarThetaPoints
                    % polar is from base not the pole
                    pj = j;
%                    pj = obj.num_polarThetaPoints + 1 - j;
                    dat(j, i) = obj.scalarDat(pj, ai);
                end
            end
%            [x1,y1,z1,c] = sphere3d(dat, 0, 2*pi, 0, pi/2, 1, 1.0, 'off', 'spline', 0.001); 
            [x1,y1,z1,c] = sphere3d(dat, 0, 2*pi, 0, pi/2, 1, 1.0, 'off', 'spline', 0.001); 
            s = surf(x1,y1,z1,c); 
%            s.EdgeColor = 'none';
            colorbar;
            fs = 18;
            hx = xlabel('$$x$$', 'interpreter', 'latex', 'FontSize', fs);
            set(hx, 'VerticalAlignment','Top');
            hy = ylabel('$$y$$', 'interpreter', 'latex', 'FontSize', fs);
            set(hy, 'VerticalAlignment','Bottom');
            hz = zlabel('$$z$$', 'interpreter', 'latex', 'FontSize', fs);
            set(hz, 'VerticalAlignment','Bottom');
        end
        
        function objout = ComputeStatistics(obj)
            objout = obj;
            for xi = 1:2
                objout.extermumLocDeg{xi} = zeros(2, 1);
            end
            sm = 0.0;
            objout.extremumV(1) =  inf;
            objout.extremumV(2) = -inf;
            
            J = objout.del_azimuthalPhisRad * objout.del_polarThetasRad;

            sm = 0.0;
            smArea = 0.0;
            for i = 1:objout.num_polarThetaPoints - 1
                s_polar_i = objout.s_polarThetasRad(i);
                s_polar_ip1 = objout.s_polarThetasRad(i + 1);
                
                for j = 1:objout.num_azimuthalPhiPoints
                    jp1 = j + 1;
                    if (j == objout.num_azimuthalPhiPoints)
                        jp1 = 1;
                    end
                    vij = objout.scalarDat(i, j);
                    vip1j = objout.scalarDat(i + 1, j);
                    vijp1 = objout.scalarDat(i, jp1);
                    vip1jp1 = objout.scalarDat(i + 1, jp1);
                    
                    sm = sm + 0.25 * J * (s_polar_i * (vij + vijp1) + s_polar_ip1 * (vip1j + vip1jp1));
                    smArea = smArea + 0.5 * J * (s_polar_i + s_polar_ip1);
                end
            end
            mnV = sm / smArea;
            objout.meanV = mnV;
            sm = 0.0;
            for i = 1:objout.num_polarThetaPoints - 1
                s_polar_i = objout.s_polarThetasRad(i);
                s_polar_ip1 = objout.s_polarThetasRad(i + 1);
                
                for j = 1:objout.num_azimuthalPhiPoints
                    jp1 = j + 1;
                    if (j == objout.num_azimuthalPhiPoints)
                        jp1 = 1;
                    end
                    vij = objout.scalarDat(i, j) - mnV;
                    vip1j = objout.scalarDat(i + 1, j) - mnV;
                    vijp1 = objout.scalarDat(i, jp1) - mnV;
                    vip1jp1 = objout.scalarDat(i + 1, jp1) - mnV;
                    vij = vij * vij;
                    vip1j = vip1j * vip1j;
                    vijp1 = vijp1 * vijp1;
                    vip1jp1 = vip1jp1 * vip1jp1;
                    
                    sm = sm + 0.25 * J * (s_polar_i * (vij + vijp1) + s_polar_ip1 * (vip1j + vip1jp1));
                end
            end
            vr = sm / smArea;
            objout.sdiv = sqrt(vr);
            
            for i = 1:objout.num_polarThetaPoints
                polar_i = objout.polarThetasRad(i);
                for j = 1:objout.num_azimuthalPhiPoints
                    azimuthalPhi_j = objout.azimuthalPhisRad(j);
                    vij = objout.scalarDat(i, j);
                    
                    if (objout.extremumV(1) > vij)
                        objout.extremumV(1) = vij;
                        objout.extermumLocRad{1}(1) = polar_i; 
                        objout.extermumLocRad{1}(2) = azimuthalPhi_j;
                        objout.extermumLocIndex{1}(1) = i;
                        objout.extermumLocIndex{1}(2) = j;
                    end
                    if (objout.extremumV(2) < vij)
                        objout.extremumV(2) = vij;
                        objout.extermumLocRad{2}(1) = polar_i; 
                        objout.extermumLocRad{2}(2) = azimuthalPhi_j;
                        objout.extermumLocIndex{2}(1) = i;
                        objout.extermumLocIndex{2}(2) = j;
                    end
                end
            end
            for xi = 1:2
                objout.extermumLocDeg{xi} = 180 / pi * objout.extermumLocRad{xi};
            end
        end
        function printHeader(obj, fileName)
            fid = open(fileName, 'w');
            fprintf(fid, 'mean\tsdiv\tCOV\tm\tm_pThetaD\tm_aPhiD\tm_pThetaR\tm_aPhiR\tm_pThetaI\tm_aPhiI\tM\tM_pThetaD\tM_aPhiD\tM_pThetaR\tM_aPhiR\tM_pThetaI\tM_aPhiI\n');
        end      
        function printStat(obj, fid, fid_plotAddedName)
            fprintf(fid, '%f\t', obj.meanV);
            fprintf(fid, '%f\t', obj.sdiv);
            fprintf(fid, '%f\t', obj.sdiv / (abs(obj.meanV) + 1e-20));
            for mM = 1:2
                fprintf(fid, '%f\t', obj.extremumV(mM));
                fprintf(fid, '%f\t', obj.extermumLocDeg{mM}(1));
                fprintf(fid, '%f\t', obj.extermumLocDeg{mM}(2));
                fprintf(fid, '%f\t', obj.extermumLocRad{mM}(1));
                fprintf(fid, '%f\t', obj.extermumLocRad{mM}(2));
                fprintf(fid, '%d\t', obj.extermumLocIndex{mM}(1));
                fprintf(fid, '%d\t', obj.extermumLocIndex{mM}(2));
            end
            fprintf(fid, '\n');
            if (fid_plotAddedName > 0)
                for i = 1:obj.num_polarThetaPoints
                    for j = 1:obj.num_azimuthalPhiPoints
                        fprintf(fid_plotAddedName, '%f\t', obj.scalarDat(i, j));
                    end
                    fprintf(fid_plotAddedName, '\n');
                end
            end
        end
    end
end

