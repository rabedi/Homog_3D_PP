classdef Rotation3D
    properties
        polarThetaRadian = -10001;
        azimuthalPhiRadian;
        alphaRadian;
        Q = zeros(3, 3); % e'1, e'2, e'3
        % https://www.brown.edu/Departments/Engineering/Courses/EN224/anis_general/anis_general.htm
        % 6x6 matrix for rotating voigt C
        T;
    end
    methods
        function objout = SetAngles(obj, polarTheta, azimuthalPhi, alpha, isRad)
            if (isRad)
                objout = SetAnglesRad(obj, polarTheta, azimuthalPhi, alpha);
            else
                objout = SetAnglesDeg(obj, polarTheta, azimuthalPhi, alpha);
            end                
        end
        function objout = SetAnglesRad(obj, polarThetaRad, azimuthalPhiRad, alphaRad)
            obj.polarThetaRadian = polarThetaRad;
            obj.azimuthalPhiRadian = azimuthalPhiRad;
            obj.alphaRadian = alphaRad;
            objout = ComputeQFrom3Angles(obj);         
        end
        function objout = SetAnglesDeg(obj, polarThetaDeg, azimuthalPhiDeg, alphaDeg)
            fact = pi / 180;
            obj.polarThetaRadian = fact * polarThetaDeg;
            obj.azimuthalPhiRadian = fact * azimuthalPhiDeg;
            obj.alphaRadian = fact * alphaDeg;
            objout = ComputeQFrom3Angles(obj);         
        end
        function objout = ComputeQFrom3Angles(obj)
            objout = obj;
            cTheta = cos(obj.polarThetaRadian); sTheta = sin(obj.polarThetaRadian);
            cPhi = cos(obj.azimuthalPhiRadian); sPhi = sin(obj.azimuthalPhiRadian);
            cAlpha = cos(obj.alphaRadian); sAlpha = sin(obj.alphaRadian);
            cThetacPhi = cTheta * cPhi;
            cThetasPhi = cTheta * sPhi;
            objout.Q(1, 1) = -(cAlpha * sPhi + sAlpha * cThetacPhi);
            objout.Q(1, 2) =   cAlpha * cPhi - sAlpha * cThetasPhi;
            objout.Q(1, 3) = sAlpha * sTheta;
            
            objout.Q(2, 1) =   sAlpha * sPhi - cAlpha * cThetacPhi;
            objout.Q(2, 2) = -(sAlpha * cPhi + cAlpha * cThetasPhi);
            objout.Q(2, 3) = cAlpha * sTheta;
            
            objout.Q(3, 1) = sTheta * cPhi;
            objout.Q(3, 2) = sTheta * sPhi;
            objout.Q(3, 3) = cTheta;
        end
        function [stressVoigt, stress, objout] = getVoigtStressCorersponding2Sigma_n_polar_azimuthal_dir(obj, polarTheta, azimuthalPhi, isRad)
            alpha = 0.0;
            if (obj.polarThetaRadian > -10000)
                objout = obj;
            else
                objout = obj.SetAngles(polarTheta, azimuthalPhi, alpha, isRad);
            end
            stressPrime = zeros(3, 3);
            % 3'3' component is 1
            %stressPrime(3, 3) = 1.0;
            %s_ij = Qmi Qnj s'_mn
            stress = zeros(3, 3);
            for i = 1:3
                for j = 1:3
                    stress(i, j) = objout.Q(3, i) * objout.Q(3, j);
                end
            end
            stressVoigt = zeros(6, 1);
            stressVoigt(1) = stress(1, 1);
            stressVoigt(2) = stress(2, 2);
            stressVoigt(3) = stress(3, 3);
            stressVoigt(4) = stress(2, 3);
            stressVoigt(5) = stress(1, 3);
            stressVoigt(6) = stress(1, 2);
        end
    end
end