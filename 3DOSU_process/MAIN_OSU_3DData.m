function validSVENos = MAIN_OSU_3DData(configNameIn, root_folderIn, config_enNoIn)
addpath('../3DCsn');
if nargin < 1
    configNameIn = 'config/config.txt';
end
if nargin < 2
    root_folderIn = 'fiber_y'; %'fiber_y_result';
    root_folderIn = 'none'; %'fiber_y_result';
end
if nargin < 3
    config_enNoIn = -1;
end
validSVENos = [];

cnf = OSU_3DData;
cnf = cnf.Read_config_OSU_3DData(configNameIn, root_folderIn, config_enNoIn); 
%sveNo = -1;
%tic
%[success_mode, C, sns] = cnf.ComputeSavePlot_OneSVE(sveNo);
%toc


tic
validSVENos = cnf.ComputeSavePlot_AllSVE();
toc

cnf.ComputeStat_PlotSummarizedResults();
