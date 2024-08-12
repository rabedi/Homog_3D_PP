classdef OSU_3DData
    properties
        stiffnessNormalizer = 1;
        interfaceStrength = 1;
        % mohr-COulomb versus sqrt root model
        strengthModelIsMC = 1;
        % effective stress = tau + k sn
        MC_k = 0.3;
        square_beta = 1;
        strengthType = 5;
        % 1->matrix
        % traction of matrix is used for strength calculation
        % 2->fiber
        % traction of fiber is used for strength calculation
        % 3->aveTrac
        % average traction of matrix and fiber is used for strength calculation
        % 4->aveEffStress
        % effective stress is calculated separately from fiber and matrix. Finally, the average of the two are taken
        % 5->minEffStress
        % minimum of strengths calculated from fiber and matrix is taken 
        
        
        
        num_fibers2Read = -1; % means reading all
        num_polarThetaSegments = 30;
        num_azimuthalPhiSegments = 120;
        AnisoIndexEq_max2minMinus1 = 0;
        printScalarVals = 1;
        polarFromPole = 1;
        % 0->no,1->yes,2->onlyForFreshC,snComputation
        b_plot_print = 1;
        plot_pdfs = 1;
        root_path = '../3DOSU_data';
        root_folder = 'fiber_y';
        output_root_path = '../3DOSU_data_out';
        SVE_startNo = 0;
        SVE_endNo = 1000;

        %% computed
        num_polarThetaPoints;
        polarThetas;
        num_azimuthalPhiPoints;
        azimuthalPhis;
        
        completeRoot_input;
        completeRoot_output;
    end
    methods
        function objout = Read_config_OSU_3DData(obj, configName, root_folderIn, config_enNoIn)
            objout = obj;
            fid = fopen(configName, 'r');
            buf = fscanf(fid, '%s', 1);
            objout.stiffnessNormalizer = fscanf(fid, '%f', 1);
            buf = fscanf(fid, '%s', 1);
            objout.interfaceStrength = fscanf(fid, '%f', 1);
            buf = fscanf(fid, '%s', 1);
            objout.strengthModelIsMC = fscanf(fid, '%d', 1);
            buf = fscanf(fid, '%s', 1);
            objout.MC_k = fscanf(fid, '%f', 1);
            buf = fscanf(fid, '%s', 1);
            objout.square_beta = fscanf(fid, '%f', 1);
                
            buf = fscanf(fid, '%s', 1);
            objout.strengthType = fscanf(fid, '%d', 1);

            buf = fscanf(fid, '%s', 1);
            objout.num_polarThetaSegments = fscanf(fid, '%d', 1);

            objout.num_polarThetaPoints = objout.num_polarThetaSegments + 1;
            objout.polarThetas = zeros(objout.num_polarThetaPoints, 1);
            for i = 1:objout.num_polarThetaPoints
                objout.polarThetas(i) = (i - 1) * 0.5 * pi / objout.num_polarThetaSegments;
            end            
            
            buf = fscanf(fid, '%s', 1);
            objout.num_azimuthalPhiSegments = fscanf(fid, '%d', 1);

            objout.num_azimuthalPhiPoints = objout.num_azimuthalPhiSegments;
            objout.azimuthalPhis = zeros(objout.num_azimuthalPhiPoints, 1);
            for j = 1:objout.num_azimuthalPhiPoints
                objout.azimuthalPhis(j) = (j - 1) * 2.0 * pi / objout.num_azimuthalPhiSegments;
            end            
            
            buf = fscanf(fid, '%s', 1);
            objout.num_fibers2Read = fscanf(fid, '%d', 1);
            
            buf = fscanf(fid, '%s', 1);
            objout.AnisoIndexEq_max2minMinus1 = fscanf(fid, '%d', 1);
            buf = fscanf(fid, '%s', 1);
            objout.printScalarVals = fscanf(fid, '%d', 1);
            buf = fscanf(fid, '%s', 1);
            objout.polarFromPole = fscanf(fid, '%d', 1);
            buf = fscanf(fid, '%s', 1);
            objout.b_plot_print = fscanf(fid, '%d', 1);
            buf = fscanf(fid, '%s', 1);
            objout.plot_pdfs = fscanf(fid, '%d', 1);
            buf = fscanf(fid, '%s', 1);
            objout.root_path = fscanf(fid, '%s', 1);
            buf = fscanf(fid, '%s', 1);
            objout.root_folder	 = fscanf(fid, '%s', 1);
            if (strcmp(root_folderIn, 'none') == 0)
                objout.root_folder = root_folderIn;
            end
            buf = fscanf(fid, '%s', 1);
            objout.output_root_path = fscanf(fid, '%s', 1);
            buf = fscanf(fid, '%s', 1);
            objout.SVE_startNo = fscanf(fid, '%d', 1);
            buf = fscanf(fid, '%s', 1);
            objout.SVE_endNo = fscanf(fid, '%d', 1);
            if (config_enNoIn > 0)
                objout.SVE_endNo = config_enNoIn;
            end
            objout.completeRoot_input = [objout.root_path, '/', objout.root_folder, '_result/result/', objout.root_folder, '_'];
            objout.completeRoot_output = [objout.output_root_path, '/', objout.root_folder];
            [status,msg,msgID] = mkdir(objout.completeRoot_output);
            objout.completeRoot_output = [objout.completeRoot_output, '/', objout.root_folder, '_'];
            fclose(fid);
        end
        function validSVENos = ComputeStat_PlotSummarizedResults(obj)
            rootPlusStartFileName = '../3DOSU_data_out/fiber_3060/fiber_3060_';
            % plotting summary files, ...
            extensions = {'sum', 'bsum', 'csum', 'vcAllC6', 'vcAllC21'};
            printText = 1;
            pdf_theta_phiLookFromTop = obj.polarFromPole;            
            for i = 1:length(extensions)
                ComputeWrite_mean_sdiv_min_max_Values_AllFiles_sameExtension(obj.completeRoot_output, extensions{i}, printText, obj.plot_pdfs, pdf_theta_phiLookFromTop);
            end
            % compute overall stats from matrix data
            Compute_Statplot2D_thetaPolar_phiAzimuthalData_AllFiles(obj.completeRoot_output, ...
                obj.polarFromPole, obj.SVE_endNo, obj.SVE_startNo);%, validIndices)
        end        
        function validSVENos = ComputeSavePlot_AllSVE(obj)
            validSVENos = [];
            cntr = 0;
            for sveNo = obj.SVE_startNo:obj.SVE_endNo
            	[success_mode, C, sns] = obj.ComputeSavePlot_OneSVE(sveNo);
                if (success_mode ~= 0)
                    cntr = cntr + 1;
                    validSVENos(cntr) = sveNo;
                end
            end
        end        
        % success_mode: 0 input does not exist, 1 exists and processed, 2
        % existed and already processed results read
        function [success_mode, C, sns] = ComputeSavePlot_OneSVE(obj, sveNo)
            if(sveNo < 0)
                sveNo = obj.SVE_startNo;
            end
            [success_mode, C, sns] = obj.Compute_OneSVE(sveNo);
            if ((success_mode == 0) || ((success_mode == 2) && (obj.b_plot_print == 2)))
                return;
            end
            
            Csn = C3D_sns;
            AnisoIndexEq_max2minMinus1 = 0;
            %Csn = Csn.Compute_C3D_sns(CIn, sns, AnisoIndexEq_max2minMinus1);
            sn = num2str(sveNo);
            saveRootName = obj.completeRoot_output;
            plotAddedName = ['_', sn];
            b_plot = 1;
            b_print = 1;
            Csn = Csn.ComputePlotPrint_C3D_sns(C, sns, obj.AnisoIndexEq_max2minMinus1, saveRootName, plotAddedName, obj.printScalarVals, obj.polarFromPole, b_plot, b_print);
            fclose('all');
            close('all');
        end
        % success_mode: 0 input does not exist, 1 exists and processed, 2
        % existed and already processed results read
        function [success_mode, C, sns] = Compute_OneSVE(obj, sveNo)
            success_mode = 0;
            C = zeros(6, 6);
            sn = num2str(sveNo);
            for i = 1:6
                fileName = [obj.completeRoot_input, sn, '_', num2str(i), '.txt'];
                fids{i} = fopen(fileName, 'r');
                if (fids{i} < 0)
                    sns = zeros(obj.num_polarThetaPoints, obj.num_azimuthalPhiPoints);
                    return;
                else
                    fclose(fids{i});
                end
            end
            fileName_C = [obj.completeRoot_output, '', sn, '_C.ascii'];
            fileName_sn = [obj.completeRoot_output, '', sn, '_sn', '.ascii'];

            fidC = fopen(fileName_C, 'r');
            if (fidC > 0)
                C = load(fileName_C, '-ascii');
                fclose(fidC);
                sns = load(fileName_sn, '-ascii');
                success_mode = 2;
            else
                command = ['OSU_3DDataCpp.exe  -sven ', num2str(sveNo)];
                status = system(command);
                fidC = fopen(fileName_C, 'r');
                if (fidC < 0)
                    fprintf(1, 'cannot open file\t%s\n', fileName_C);
                end
                C = load(fileName_C, '-ascii');
                fclose(fidC);
                sns = load(fileName_sn, '-ascii');
                success_mode = 1;
            end
            C = C / obj.stiffnessNormalizer;
            sns = sns * obj.interfaceStrength;
        end
        function [voigtStresses4_polar_azimuthalAngles] = Generate_voigtStresses4_polar_azimuthalAngles(obj)
            fileName = ['_voigtStresses4_polar_azimuthalAngles', num2str(obj.num_polarThetaSegments), ...
                '_', num2str(obj.num_azimuthalPhiSegments), '.txt'];
            voigtStresses4_polar_azimuthalAngles = cell(obj.num_polarThetaPoints, obj.num_azimuthalPhiPoints);
            fidr = fopen(fileName, 'r');
            if (fidr > 0)
                for i = 1:obj.num_polarThetaPoints
                    for j = 1:obj.num_azimuthalPhiPoints
                        buf = fscanf(fidr, '%d', 1);
                        buf = fscanf(fidr, '%d', 1);
                        for k = 1:6
                            tmp(k) = fscanf(fidr,'%f', 1);
                        end
                        voigtStresses4_polar_azimuthalAngles{i,j} = tmp;
                    end
                end
                fclose(fidr);
            end
            fidw = fopen(fileName, 'w');
            for i = 1:obj.num_polarThetaPoints
                polarTheta = obj.polarThetas(i);
                for j = 1:obj.num_azimuthalPhiPoints
                    azimuthalPhi = obj.azimuthalPhis(j);
                    r3D = Rotation3D;
                    [stressVoigt, stress, r3D] = r3D.getVoigtStressCorersponding2Sigma_n_polar_azimuthal_dir(polarTheta, azimuthalPhi, 1);
                    voigtStresses4_polar_azimuthalAngles{i,j} = stressVoigt;
                    fprintf(fidw, '%d\t%d', i, j);
                    for k = 1:6
                        fprintf(fidw, '\t%f', stressVoigt(k));
                    end
                    fprintf(fidw, '\n');
                end
            end
            fclose(fidw);
        end
            
            
            
        % success_mode: 0 input does not exist, 1 exists and processed, 2
        % existed and already processed results read
        function [success_mode, C, sns] = Compute_OneSVE_obsoleteMatlabImplementationSuperSlow(obj, sveNo)
            success_mode = 0;
            C = zeros(6, 6);
            sn = num2str(sveNo);
            for i = 1:6
                fileName = [obj.completeRoot_input, sn, '_', num2str(i), '.txt'];
                fids{i} = fopen(fileName, 'r');
                if (fids{i} < 0)
                    sns = zeros(obj.num_polarThetaPoints, obj.num_azimuthalPhiPoints);
                    return;
                end
            end

            %   actual strength calculations
            % mat -> 1, fiber -> 2, ave -> 3
            sns_tractypes = cell(3, 1);
            for tracttype = 1:3
                sns_tractypes{tracttype} = inf * ones(obj.num_polarThetaPoints, obj.num_azimuthalPhiPoints);
            end
            
            fileName_C = [obj.completeRoot_output, '', sn, '_C.ascii'];
            for tracttype = 1:3
                fileName_sn{tracttype} = [obj.completeRoot_output, '', sn, '_sn', num2str(tracttype), '.ascii'];
            end
            if (obj.strengthType < 4)
                tracttypeMin = obj.strengthType;
                tracttypeMax = obj.strengthType;
            else
                tracttypeMin = 1;
                tracttypeMax = 2;
            end                
                
            
            fidC = fopen(fileName_C, 'r');
            if (fidC > 0)
                C = load(fileName_C, '-ascii');
                fclose(fidC);
                for tracttype = tracttypeMin:tracttypeMax
                    sns_tractypes{tracttype} = load(fileName_sn{tracttype}, '-ascii');
                end
                success_mode = 2;
            else
                success_mode = 1;
                % order of stress and strain is:
                % s11, s22, s33, s23, s31, s12
                % eps11, eps22, eps33, 2eps23, 2eps31, 2eps12
                % the way strains are written, strains of 6 load cases generate
                % strains = 0.001 * I_6x6
                % stresses are read and it's noted that 
                % Pengfei I check my code and find out that the stress component 5 and 6 should be swapped.
                aveStress = zeros(6, 6);
                for i = 1:6 
                    buf = fgetl(fids{i});
                    buf = fgetl(fids{i});
                    aveStress(1, i) = fscanf(fids{i}, '%f', 1); % sigma11
                    aveStress(2, i) = fscanf(fids{i}, '%f', 1); % sigma22
                    aveStress(3, i) = fscanf(fids{i}, '%f', 1); % sigma33
                    aveStress(6, i) = fscanf(fids{i}, '%f', 1); % tau12
                    aveStress(5, i) = fscanf(fids{i}, '%f', 1); % tau13 (incorrectly labeled as tau23)
                    aveStress(4, i) = fscanf(fids{i}, '%f', 1); % tau23 (incorrectly labeled as tau13)
                    buf = fscanf(fids{i}, '%s', 1);
                    buf = fgetl(fids{i});
                    buf = fgetl(fids{i});
                    buf = fgetl(fids{i});
                    buf = fgetl(fids{i});
%                    buf = fgetl(fids{i});
                end
                
                % computing strengths
                %   load factors...
                aveStressInv = inv(aveStress);

                [voigtStresses4_polar_azimuthalAngles] = obj.Generate_voigtStresses4_polar_azimuthalAngles();
                stressLoadFactors = cell(obj.num_polarThetaPoints, obj.num_azimuthalPhiPoints);
                for i = 1:obj.num_polarThetaPoints
    %                polarTheta = obj.polarThetas(i);
                    for j = 1:obj.num_azimuthalPhiPoints
    %                    azimuthalPhi = obj.azimuthalPhis(j);
                         voigtStresses4_polar_azimuthal = voigtStresses4_polar_azimuthalAngles{i, j};
                         stressLoadFactors{i, j} = aveStressInv * voigtStresses4_polar_azimuthal;
                    end
                end

                strengthModelIsMC = obj.strengthModelIsMC;
                MC_k = obj.MC_k;
                MC_kInv = 1.0 / MC_k;
                square_beta = obj.square_beta;

                if (obj.num_fibers2Read < 0)
                    for fi = 1:6
                        vlRaw{fi} = fscanf(fids{fi}, '%f', [9, inf])';
                        fclose(fids{fi});
                    end
                else
                    for fi = 1:6
                        vlRaw{fi} = fscanf(fids{fi}, '%f', [9, obj.num_fibers2Read])';
                        fclose(fids{fi});
                    end
                end
                
                [numFibers, b] = size(vlRaw{1});
                
%                while (true)
                 for fiberCntr = 1:numFibers
                    for fi = 1:6
                        vals = vlRaw{fi}(fiberCntr, :);
%                         buf = fgetl(fids{fi});
%                         vals = str2num(buf);
%                         if (length(vals) == 0)
%                             break;
%                         end
                         matrixTraction_x = vals(1);
                         matrixTraction_y = vals(2);
                         matrixTraction_z = vals(3);

                         fiberTraction_x = vals(4);
                         fiberTraction_y = vals(5);
                         fiberTraction_z = vals(6);

                         normal_x = vals(7);
                         normal_y = vals(8);
                         normal_z = vals(9);
                         
%                         [matrixTraction_x, sucRead] = fscanf(fids{fi}, '%f', 1);
%                         if (~sucRead)
%                             break;
%                         end
%                         matrixTraction_y = fscanf(fids{fi}, '%f', 1);
%                         matrixTraction_z = fscanf(fids{fi}, '%f', 1);
% 
%                         fiberTraction_x = fscanf(fids{fi}, '%f', 1);
%                         fiberTraction_y = fscanf(fids{fi}, '%f', 1);
%                         fiberTraction_z = fscanf(fids{fi}, '%f', 1);
% 
%                         normal_x = fscanf(fids{fi}, '%f', 1);
%                         normal_y = fscanf(fids{fi}, '%f', 1);
%                         normal_z = fscanf(fids{fi}, '%f', 1);

                        matrix_trac = [matrixTraction_x, matrixTraction_y, matrixTraction_z]';
                        fiber_trac = [fiberTraction_x, fiberTraction_y, fiberTraction_z]';
                        ave_trac = 0.5 * (matrix_trac + fiber_trac);

                        trac{fi}{1} = matrix_trac;
                        trac{fi}{2} = fiber_trac;
                        trac{fi}{3} = ave_trac;

                        normal = [normal_x, normal_y, normal_z]';
                    end

                    for i = 1:obj.num_polarThetaPoints
                        for j = 1:obj.num_azimuthalPhiPoints
                             stressLoadFactor = stressLoadFactors{i, j};
                             tracFinal = cell(3, 1);
                             for tracttype = tracttypeMin:tracttypeMax
                                 tracFinal{tracttype} = zeros(3, 1);
                             end
                             for fi = 1:6
                                 lc = stressLoadFactor(fi);

                                 for tracttype = tracttypeMin:tracttypeMax
                                     tmpVec = trac{fi}{tracttype};
                                     for k = 1:3
                                        tracFinal{tracttype}(k) = tracFinal{tracttype}(k) + lc * tmpVec(k);
                                     end
                                 end
                             end

                             % compute strength
                             for tracttype = tracttypeMin:tracttypeMax
                                 traction = tracFinal{tracttype};
                                 sn = traction' * normal;
                                 tau = traction - sn * normal;
                                 abs_tau = norm(tau);

                                 % fracture model
                                 if (strengthModelIsMC)
                                     effStress = sn + abs_tau * MC_kInv;
    %                                 effStress = MC_k * sn + abs_tau;
                                 else
                                     effStress = square_beta * abs_tau;
                                     effStress = effStress * effStress;
                                     if (sn > 0)
                                         effStress = effStress + sn * sn;
                                     end
                                     effStress = sqrt(effStress);
                                 end
                                 if (effStress > 0)
                                     strengthThisAngle = 1 / effStress;
                                     strengthStoredSoFar = sns_tractypes{tracttype}(i, j);
                                     if (strengthThisAngle < strengthStoredSoFar)
                                         sns_tractypes{tracttype}(i, j) = strengthThisAngle;
                                     end
                                 end
                             end
                        end
                    end
                    if (mod(fiberCntr, 100000) == 0)
                        fprintf(1, '%d\t', fiberCntr);
                    end
                end
            end
            if (obj.strengthType == 1) % matrix values
                sns = sns_tractypes{1};
            elseif (obj.strengthType == 2) % fiber values
                sns = sns_tractypes{2};
            elseif (obj.strengthType == 3) % average traction is used
                sns = sns_tractypes{3};
            elseif (obj.strengthType == 4) % average of strengths
                sns = 0.5 * (sns_tractypes{1} + sns_tractypes{2});
            elseif (obj.strengthType == 5) % min of strengths
                sns = min(sns_tractypes{1}, sns_tractypes{2});
            end
            
            % stiffness
            if (success_mode == 1)
                C = aveStress / 0.001;
                save(fileName_C, 'C', '-ascii');
            end
            C = C / obj.stiffnessNormalizer;
            if (success_mode == 1)
                for tracttype = tracttypeMin:tracttypeMax
                    tmpmat = sns_tractypes{tracttype};
                    save(fileName_sn{tracttype}, 'tmpmat', '-ascii');
                end
            end
            sns = sns * obj.interfaceStrength;
        end
        % end of compute SVE obsolete (super slow, ~ 5 hour for 30 x 120
        % mesh strength calculation)
    end
end

