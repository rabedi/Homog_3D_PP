function ExtractPlotSphericalAngleData(fileName, extension, fileName4StatWOExt, angleInDeg, startAngleCols, addedAngleNames, plotPDFs, pdf_theta_phiLookFromTop)
fs = 18;

if nargin < 1
%    fileName = '../3DOSU_data_out/fiber_z/fiber_z_sn_normal';
    fileName = '../3DOSU_data_out/fiber_3060/fiber_3060_sn_normal';
end
if nargin < 2
    extension = 'sum';
end

if nargin < 3
    fileName4StatWOExt = 'none';
end

if nargin < 4
    angleInDeg = 1;
end

if nargin < 5
    startAngleCols = [5, 12];
end

if nargin < 6
    addedAngleNames = { 'min', 'Max'};
end

if nargin < 7
    plotPDFs = 1;
end

if nargin < 8
    pdf_theta_phiLookFromTop = 1;
end

[a, fieldName] = fileparts(fileName);


polarSphericalAngles = extractSphericalAngleData(fileName, extension, startAngleCols, angleInDeg);
sz = length(polarSphericalAngles);
if (sz == 0)
    return;
end

for ds = 1:sz
    added_s = addedAngleNames{ds};
    fnfld = [fileName, '_', added_s];
    polarSphericalAngleSet = polarSphericalAngles{ds};
    [meanAngles, r, rBased_sdiv, marginal_pdf_polar, marginal_pdf_azimuthal] = SphericalAngleStatisticalAnalysis(polarSphericalAngleSet, angleInDeg, plotPDFs);

    
    if (strcmp(fileName4StatWOExt, 'none') == 1)
        fileNameStat = [fileName, '_circular_stat.stat'];
    else
        fileNameStat = [fileName4StatWOExt , '_', added_s, '_circular_stat.stat'];
    end
    fids = fopen(fileNameStat, 'a');
    fprintf(fids, '%s\t%s\t', fieldName, added_s);
    fprintf(fids, 'mean\tazimuthal\t%f\tpolar\t%f\t', meanAngles(2), meanAngles(1));
    fprintf(fids, 'r\t%f\trBased_sdiv \t%f\n', r, rBased_sdiv );
    fclose(fids);
    
    if (~plotPDFs)
        continue;
    end
    gridx1 = 0:5:360;
    gridx2 = 0:5:90;
    if (~angleInDeg)
        gridx1 = pi / 180.0 * gridx1;
        gridx2 = pi / 180.0 * gridx2;
    end
    alim = [gridx1(1), gridx1(length(gridx1))];
    plim = [gridx2(1), gridx2(length(gridx2))];
    [x1,x2] = meshgrid(gridx1, gridx2);
    x1 = x1(:);
    x2 = x2(:);
    xi = [x1 x2];
    densityVals = polarSphericalAngleSet;
    densityVals(:, 1) = polarSphericalAngleSet(:, 2);
    densityVals(:, 2) = polarSphericalAngleSet(:, 1);
    ksdensity(densityVals, xi, 'PlotFcn', 'surf');
    if (pdf_theta_phiLookFromTop)
        view(0,90);
    end
    colorbar;
    hx = xlabel('$$\phi$$ (Azimuthal)', 'interpreter', 'latex', 'FontSize', fs);
    set(hx, 'VerticalAlignment','Top');
    hy = ylabel('$$\theta$$ (Polar)', 'interpreter', 'latex', 'FontSize', fs);
    set(hy, 'VerticalAlignment','Bottom');
    xlim(alim);
    ylim(plim);
    fnbase = [fnfld, '_pdf2_azimuthal_polar']; 
    fn = [fnbase, '.', 'png'];
    print('-dpng', fn);
    fn = [fnbase, '.', 'fig'];
    save(fn);

    plot(marginal_pdf_azimuthal{1}, marginal_pdf_azimuthal{2});
    hx = xlabel('$$\phi$$ (Azimuthal)', 'interpreter', 'latex', 'FontSize', fs);
    set(hx, 'VerticalAlignment','Top');
    hy = ylabel('PDF', 'interpreter', 'latex', 'FontSize', fs);
    set(hy, 'VerticalAlignment','Bottom');
    xlim(alim);
    fnbase = [fnfld, '_pdf1_azimuthal']; 
    fn = [fnbase, '.', 'png'];
    print('-dpng', fn);
    fn = [fnbase, '.', 'fig'];
    save(fn);

    plot(marginal_pdf_polar{1}, marginal_pdf_polar{2});
    hx = xlabel('$$\theta$$ (Polar)', 'interpreter', 'latex', 'FontSize', fs);
    set(hx, 'VerticalAlignment','Top');
    hy = ylabel('PDF', 'interpreter', 'latex', 'FontSize', fs);
    set(hy, 'VerticalAlignment','Bottom');
    xlim(plim);
    fnbase = [fnfld, '_pdf2_polar']; 
    fn = [fnbase, '.', 'png'];
    print('-dpng', fn);
    fn = [fnbase, '.', 'fig'];
    save(fn);
end



