function [X, Y, val, successful] = plot2D_thetaPolar_phiAzimuthalData(fileName, extension, polarFromPole)
fn = [fileName, '.', extension];
fid = fopen(fn, 'r');
if (fid < 0)
    successful = 0;
    X = [];
    Y = [];
    val = [];
end
successful = 1;
scalarDat = load(fn, '-ascii'); 
fclose(fid);
[num_polarThetaPoints, num_azimuthalPhiPoints] = size(scalarDat);
num_polarThetaSegments = num_polarThetaPoints - 1;
num_azimuthalPhiSegments = num_azimuthalPhiPoints;
azimuthal = zeros(num_azimuthalPhiSegments, 1);
factor = 360 / num_azimuthalPhiPoints;
for i = 1:num_azimuthalPhiPoints + 1
    azimuthal(i) = (i - 1) * factor; 
end
polar = zeros(num_polarThetaPoints, 1);
factor = 90 / num_polarThetaSegments;
for i = 1:num_polarThetaPoints
    polar(i) = (i - 1) * factor; 
end
[X, Y] = meshgrid(azimuthal, polar);
for i = 1:num_azimuthalPhiPoints + 1
    ai = i;
    if (ai == num_azimuthalPhiPoints + 1)
        ai = 1;
    end
    for j = 1:num_polarThetaPoints
        pj = j;
        if (~polarFromPole)
            pj = num_polarThetaPoints + 1 - j;
        end
        val(j, i) = scalarDat(pj, ai);
    end
end
contourf(X,Y,val, 'EdgeColor', 'none');
colorbar;
fs = 18;
hx = xlabel('$$\phi$$ (Azimuthal)', 'interpreter', 'latex', 'FontSize', fs);
set(hx, 'VerticalAlignment','Top');
hy = ylabel('$$\theta$$ (Polar)', 'interpreter', 'latex', 'FontSize', fs);
set(hy, 'VerticalAlignment','Bottom');
fn = [fileName, '.', 'png'];
print('-dpng', fn);
fn = [fileName, '.', 'fig'];
save(fn);