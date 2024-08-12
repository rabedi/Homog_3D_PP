
addpath('../3DCsn');
fileName = '../3DOSU_data_out/fiber_3060/fiber_3060_sn_normal_0';
extension = 'mt';
polarFromPole = 1;
[X, Y, val, successful] = plot2D_thetaPolar_phiAzimuthalData(fileName, extension, polarFromPole);


ExtractPlotSphericalAngleData()
a = 12;

% gridx1 = -0.25:.05:1.25;
% gridx2 = 0:.1:15;
% [x1,x2] = meshgrid(gridx1, gridx2);
% x1 = x1(:);
% x2 = x2(:);
% xi = [x1 x2];
% 
% rng('default')  % For reproducibility
% x = [0+.5*rand(20,1) 5+2.5*rand(20,1);
%             .75+.25*rand(10,1) 8.75+1.25*rand(10,1)];
%         
% figure
% %ksdensity(x,xi);        
% ksdensity(x);        


angleInDeg = 1;
polarSphericalAngles = extractSphericalAngleData();


gridx1 = 0:5:360;
gridx2 = 0:5:90;
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
densityVals = polarSphericalAngles{1};
densityVals(:, 1) = polarSphericalAngles{1}(:, 2);
densityVals(:, 2) = polarSphericalAngles{1}(:, 1);

ksdensity(densityVals, xi, 'PlotFcn', 'surf');
view(0,90)
colorbar;
fs = 18;
hx = xlabel('$$\phi$$ (Azimuthal)', 'interpreter', 'latex', 'FontSize', fs);
set(hx, 'VerticalAlignment','Top');
hy = ylabel('$$\theta$$ (Polar)', 'interpreter', 'latex', 'FontSize', fs);
set(hy, 'VerticalAlignment','Bottom');
xlim([0, 360]);
ylim([0, 90]);
%fn = [fileName, '.', 'png'];
%print('-dpng', fn);
%fn = [fileName, '.', 'fig'];
%save(fn);

[meanAngles, r, rBased_sdiv, marginal_pdf_polar, marginal_pdf_azimuthal] = SphericalAngleStatisticalAnalysis(polarSphericalAngles{1}, angleInDeg);
meanAngles, r, rBased_sdiv
figure(1);
plot(marginal_pdf_polar{1}, marginal_pdf_polar{2});
figure(2);
plot(marginal_pdf_azimuthal{1}, marginal_pdf_azimuthal{2});


a = 12;