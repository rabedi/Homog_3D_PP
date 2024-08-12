classdef C3D
    properties
        % 6x6 Voigt notation, [eps11, eps22, eps33, 2*eps23, 2*eps31, 2*eps12]
        % s similar without factors of 2
        C;
        C4thOrder;
    end
    methods
        function objout = SetC(obj, CVoigtIn, compute4thOrder)
            objout = obj;
            objout.C = CVoigtIn;
            if (compute4thOrder)
                objout.C4thOrder = zeros(3, 3, 3, 3);
                objout.C4thOrder(1, 1, 1, 1) = CVoigtIn(1, 1);
                objout.C4thOrder(2, 2, 2, 2) = CVoigtIn(2, 2);
                objout.C4thOrder(3, 3, 3, 3) = CVoigtIn(3, 3);
                tmp = CVoigtIn(1, 2);
                objout.C4thOrder(1, 1, 2, 2) = tmp;
                objout.C4thOrder(2, 2, 1, 1) = tmp;
                tmp = CVoigtIn(2, 3);
                objout.C4thOrder(2, 2, 3, 3) = tmp;
                objout.C4thOrder(3, 3, 2, 2) = tmp;
                tmp = CVoigtIn(1, 3);
                objout.C4thOrder(1, 1, 3, 3) = tmp;
                objout.C4thOrder(3, 3, 1, 1) = tmp;
                tmp = CVoigtIn(1, 4);
                objout.C4thOrder(1, 1, 2, 3) = tmp;
                objout.C4thOrder(1, 1, 3, 2) = tmp;
                objout.C4thOrder(2, 3, 1, 1) = tmp;
                objout.C4thOrder(3, 2, 1, 1) = tmp;
                tmp = CVoigtIn(1, 5);
                objout.C4thOrder(1, 1, 1, 3) = tmp;
                objout.C4thOrder(1, 1, 3, 1) = tmp;
                objout.C4thOrder(1, 3, 1, 1) = tmp;
                objout.C4thOrder(3, 1, 1, 1) = tmp;
                tmp = CVoigtIn(1, 6);
                objout.C4thOrder(1, 1, 1, 2) = tmp;
                objout.C4thOrder(1, 1, 2, 1) = tmp;
                objout.C4thOrder(1, 2, 1, 1) = tmp;
                objout.C4thOrder(2, 1, 1, 1) = tmp;
                tmp = CVoigtIn(2, 4);
                objout.C4thOrder(2, 2, 2, 3) = tmp;
                objout.C4thOrder(2, 2, 3, 2) = tmp;
                objout.C4thOrder(2, 3, 2, 2) = tmp;
                objout.C4thOrder(3, 2, 2, 2) = tmp;
                tmp = CVoigtIn(2, 5);
                objout.C4thOrder(2, 2, 1, 3) = tmp;
                objout.C4thOrder(2, 2, 3, 1) = tmp;
                objout.C4thOrder(1, 3, 2, 2) = tmp;
                objout.C4thOrder(3, 1, 2, 2) = tmp;
                tmp = CVoigtIn(2, 6);
                objout.C4thOrder(2, 2, 1, 2) = tmp;
                objout.C4thOrder(2, 2, 2, 1) = tmp;
                objout.C4thOrder(1, 2, 2, 2) = tmp;
                objout.C4thOrder(2, 1, 2, 2) = tmp;
                tmp = CVoigtIn(3, 4);
                objout.C4thOrder(3, 3, 2, 3) = tmp;
                objout.C4thOrder(3, 3, 3, 2) = tmp;
                objout.C4thOrder(2, 3, 3, 3) = tmp;
                objout.C4thOrder(3, 2, 3, 3) = tmp;
                tmp = CVoigtIn(3, 5);
                objout.C4thOrder(3, 3, 1, 3) = tmp;
                objout.C4thOrder(3, 3, 3, 1) = tmp;
                objout.C4thOrder(1, 3, 3, 3) = tmp;
                objout.C4thOrder(3, 1, 3, 3) = tmp;
                tmp = CVoigtIn(3, 6);
                objout.C4thOrder(3, 3, 1, 2) = tmp;
                objout.C4thOrder(3, 3, 2, 1) = tmp;
                objout.C4thOrder(1, 2, 3, 3) = tmp;
                objout.C4thOrder(2, 1, 3, 3) = tmp;
                tmp = CVoigtIn(4, 4);
                objout.C4thOrder(2, 3, 2, 3) = tmp;
                objout.C4thOrder(2, 3, 3, 2) = tmp;
                objout.C4thOrder(3, 2, 2, 3) = tmp;
                objout.C4thOrder(3, 2, 3, 2) = tmp;
                tmp = CVoigtIn(5, 5);
                objout.C4thOrder(1, 3, 1, 3) = tmp;
                objout.C4thOrder(1, 3, 3, 1) = tmp;
                objout.C4thOrder(3, 1, 1, 3) = tmp;
                objout.C4thOrder(3, 1, 3, 1) = tmp;
                tmp = CVoigtIn(6, 6);
                objout.C4thOrder(1, 2, 1, 2) = tmp;
                objout.C4thOrder(1, 2, 2, 1) = tmp;
                objout.C4thOrder(2, 1, 1, 2) = tmp;
                objout.C4thOrder(2, 1, 2, 1) = tmp;
                tmp = CVoigtIn(4, 5);
                objout.C4thOrder(2, 3, 3, 1) = tmp;
                objout.C4thOrder(2, 3, 1, 3) = tmp;
                objout.C4thOrder(3, 2, 3, 1) = tmp;
                objout.C4thOrder(3, 2, 1, 3) = tmp;
                objout.C4thOrder(3, 1, 2, 3) = tmp;
                objout.C4thOrder(3, 1, 3, 2) = tmp;
                objout.C4thOrder(1, 3, 2, 3) = tmp;
                objout.C4thOrder(1, 3, 3, 2) = tmp;
                tmp = CVoigtIn(4, 6);
                objout.C4thOrder(2, 3, 2, 1) = tmp;
                objout.C4thOrder(2, 3, 1, 2) = tmp;
                objout.C4thOrder(3, 2, 2, 1) = tmp;
                objout.C4thOrder(3, 2, 1, 2) = tmp;
                objout.C4thOrder(2, 1, 2, 3) = tmp;
                objout.C4thOrder(2, 1, 3, 2) = tmp;
                objout.C4thOrder(1, 2, 2, 3) = tmp;
                objout.C4thOrder(1, 2, 3, 2) = tmp;
                tmp = CVoigtIn(5, 6);
                objout.C4thOrder(1, 3, 2, 1) = tmp;
                objout.C4thOrder(1, 3, 1, 2) = tmp;
                objout.C4thOrder(3, 1, 2, 1) = tmp;
                objout.C4thOrder(3, 1, 1, 2) = tmp;
                objout.C4thOrder(2, 1, 1, 3) = tmp;
                objout.C4thOrder(2, 1, 3, 1) = tmp;
                objout.C4thOrder(1, 2, 1, 3) = tmp;
                objout.C4thOrder(1, 2, 3, 1) = tmp;
            end
        end
        % obj has C 4th, objout will get C voigt
        function objout = SetC4th2CVoigt(obj)
            objout = obj;
            objout.C(1, 1) = obj.C4thOrder(1, 1, 1, 1);
            objout.C(2, 2) = obj.C4thOrder(2, 2, 2, 2);
            objout.C(3, 3) = obj.C4thOrder(3, 3, 3, 3);
            objout.C(1, 2) = obj.C4thOrder(1, 1, 2, 2);
            objout.C(1, 3) = obj.C4thOrder(1, 1, 3, 3);
            objout.C(2, 3) = obj.C4thOrder(2, 2, 3, 3);
            objout.C(1, 4) = obj.C4thOrder(1, 1, 2, 3);
            objout.C(1, 5) = obj.C4thOrder(1, 1, 3, 1);
            objout.C(1, 6) = obj.C4thOrder(1, 1, 1, 2);
            objout.C(2, 4) = obj.C4thOrder(2, 2, 2, 3);
            objout.C(2, 5) = obj.C4thOrder(2, 2, 3, 1);
            objout.C(2, 6) = obj.C4thOrder(2, 2, 1, 2);
            objout.C(3, 4) = obj.C4thOrder(3, 3, 2, 3);
            objout.C(3, 5) = obj.C4thOrder(3, 3, 3, 1);
            objout.C(3, 6) = obj.C4thOrder(3, 3, 1, 2);
            objout.C(4, 4) = obj.C4thOrder(2, 3, 2, 3);
            objout.C(4, 5) = obj.C4thOrder(2, 3, 3, 1);
            objout.C(4, 6) = obj.C4thOrder(2, 3, 1, 2);
            objout.C(5, 5) = obj.C4thOrder(3, 1, 3, 1);
            objout.C(5, 6) = obj.C4thOrder(3, 1, 1, 2);
            objout.C(6, 6) = obj.C4thOrder(1, 2, 1, 2);
            for i = 1:6
                for j = i + 1:6
                    objout.C(j, i) = objout.C(i, j);
                end
            end
        end
        % obj is in original coordinate, objout is in rotated coordinate
        % r3DIn must have its Q set
        function [objout, r3Dout] = RotateDirectVoigt(obj, r3DIn)
            objout = obj;
            r3Dout = r3DIn;
            Qexp = [r3DIn.Q, r3DIn.Q; r3DIn.Q, r3DIn.Q];
            r3Dout.T = zeros(6, 6);
            for i = 1:3
                for j = 1:3
                    r3Dout.T(i, j) = Qexp(i, j) * Qexp(i, j);
                    r3Dout.T(i, j + 3) = 2 * Qexp(i, j + 1) * Qexp(i, j + 2);
                    r3Dout.T(i + 3, j) = Qexp(i + 1, j) * Qexp(i + 2, j);
                    r3Dout.T(i + 3, j + 3) = Qexp(i + 1, j + 1) * Qexp(i + 2, j + 2) + Qexp(i + 1, j + 2) * Qexp(i + 2, j + 1);
                end
            end
            objout.C = r3Dout.T * obj.C * transpose(r3Dout.T);
        end
        % this is just for debugging above and is not used in practice
        % Q in 3r3DIn is set
        function objout = RotateUsing4thOrderC(obj, r3DIn)
            objout = obj;
            obj = SetC(obj, obj.C, 1);
            objout.C4thOrder = zeros(3, 3, 3, 3);
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        for l = 1:3
                            for p = 1:3
                                for q = 1:3
                                    for r = 1:3
                                        for s = 1:3
                                            objout.C4thOrder(i, j, k, l) = ...
                                                objout.C4thOrder(i, j, k, l) + ...
                                                r3DIn.Q(i, p) * r3DIn.Q(j, q) * r3DIn.Q(k, r) * r3DIn.Q(l, s) * ...
                                                obj.C4thOrder(p, q, r, s);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
           	objout = SetC4th2CVoigt(objout);
        end
        function [objout, r3Dout] = RotateDirectVoigt_wAngles(obj, polarTheta, azimuthalPhi, alpha, isRad)
            r3Dout = Rotation3D;
            r3Dout = r3Dout.SetAngles(polarTheta, azimuthalPhi, alpha, isRad);
        	[objout, r3Dout] = RotateDirectVoigt(obj, r3Dout);
        end
        function [objout, r3Dout] = RotateUsing4thOrderC_wAngles(obj, polarTheta, azimuthalPhi, alpha, isRad)
            r3Dout = Rotation3D;
            r3Dout = r3Dout.SetAngles(polarTheta, azimuthalPhi, alpha, isRad);
        	objout = RotateUsing4thOrderC(obj, r3Dout);
        end
        function L2normFromVoigt = getL2NormFromVoigt(obj)
            nrm = obj.C(1, 1) * obj.C(1, 1) + obj.C(2, 2) * obj.C(2, 2) + obj.C(3, 3) * obj.C(3, 3) + ...
               2 * (obj.C(1, 2) * obj.C(1, 2) + obj.C(2, 3) * obj.C(2, 3) + obj.C(3, 1) * obj.C(3, 1)) + ...
               4 * (obj.C(1, 4) * obj.C(1, 4) + obj.C(1, 5) * obj.C(1, 5) + obj.C(1, 6) * obj.C(1, 6) + ...
               obj.C(2, 4) * obj.C(2, 4) + obj.C(2, 5) * obj.C(2, 5) + obj.C(2, 6) * obj.C(2, 6) + ...
               obj.C(3, 4) * obj.C(3, 4) + obj.C(3, 5) * obj.C(3, 5) + obj.C(3, 6) * obj.C(3, 6)) + ...
               4 * (obj.C(4, 4) * obj.C(4, 4) + obj.C(5, 5) * obj.C(5, 5) + obj.C(6, 6) * obj.C(6, 6)) + ...
               8 * (obj.C(4, 5) * obj.C(4, 5) + obj.C(4, 6) * obj.C(4, 6) + obj.C(5, 6) * obj.C(5, 6));
           L2normFromVoigt = sqrt(nrm);
        end
        function L2normFrom4th = getL2NormFrom4th(obj)
            if (length(obj.C4thOrder) == 0)
                obj = SetC(obj, obj.C, 1);
            end
            nrm = 0.0;
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        for l = 1:3
                            nrm = nrm + obj.C4thOrder(i, j, k, l) * obj.C4thOrder(i, j, k, l);
                        end
                    end
                end
            end
            L2normFrom4th = sqrt(nrm);
        end
        % computes transverse isotropic with rotation around 3rd axis
        % transverIsotropicNormalizedDiff: ||CTransvereIso - C||/||C||,
        % CTransverseIso = objout, C = obj
        % b_transverIsotropicNormalizedDiff: determines if this is computed 
        function [objout, outOfPlaneStiffness, outOfPlane2InplaneNormalStiffness, transverIsotropicNormalizedDiff, diffNorm, CNorm] = Compute_TransverseIsotropic_around3rdAxis(obj, b_transverIsotropicNormalizedDiff)
            transverIsotropicNormalizedDiff  = 0.0;
            diffNorm = 0;
            CNorm = 0;
            objout = obj;
            objout.C = zeros(6, 6);
            C11 = obj.C(1, 1);
            C22 = obj.C(2, 2);
            C12 = obj.C(1, 2);
            C66 = obj.C(6, 6);

            objout.C(3, 3) = obj.C(3, 3);
            outOfPlaneStiffness = objout.C(3, 3);
            objout.C(4, 4) = 0.5 * (obj.C(4, 4) + obj.C(5, 5));
            objout.C(5, 5) = objout.C(4, 4);

            objout.C(1, 3) = 0.5 * (obj.C(1, 3) + obj.C(2, 3));
            objout.C(2, 3) = objout.C(1, 3);
            objout.C(3, 1) = objout.C(1, 3);
            objout.C(3, 2) = objout.C(1, 3);

            Cp11 = 0.375 * (C11 + C22) + 0.25 * C12 + 0.5 * C66;
            Cp12 = 0.125 * (C11 + C22) + 0.75 * C12 - 0.5 * C66;
            Cp66 = 0.5 * (Cp11 - Cp12);
            objout.C(1, 1) = Cp11;
            objout.C(2, 2) = Cp11;
            objout.C(1, 2) = Cp12;
            objout.C(2, 1) = Cp12;
            objout.C(6, 6) = Cp66;
            outOfPlane2InplaneNormalStiffness = objout.C(3, 3) / objout.C(1, 1);
            if (b_transverIsotropicNormalizedDiff)
                c3DDiff = C3D;
                c3DDiff.C = objout.C - obj.C;
                diffNorm = c3DDiff.getL2NormFromVoigt();
                CNorm = getL2NormFromVoigt(obj);
                transverIsotropicNormalizedDiff = diffNorm / CNorm;
            end
        end
        function [rotatedC, rotatedCIso, outOfPlaneStiffness, outOfPlane2InplaneNormalStiffness, transverIsotropicNormalizedDiff, r3Dout, diffNorm, CNorm] = ComputeCRotated_CRotatedIso_etc(obj, polarTheta, azimuthalPhi, isRad, b_transverIsotropicNormalizedDiff)
            r3Dout = Rotation3D;
            alpha = 0;
            r3Dout = r3Dout.SetAngles(polarTheta, azimuthalPhi, alpha, isRad);
        	[rotatedC, r3Dout] = RotateDirectVoigt(obj, r3Dout);
            [rotatedCIso, outOfPlaneStiffness, outOfPlane2InplaneNormalStiffness, transverIsotropicNormalizedDiff, diffNorm, CNorm] = ...
                rotatedC.Compute_TransverseIsotropic_around3rdAxis(b_transverIsotropicNormalizedDiff);
        end
    end
end