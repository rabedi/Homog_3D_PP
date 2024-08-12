function polarSphericalAngles = extractSphericalAngleData(fileName, extension, startAngleCols, angleInDeg)
if nargin < 1
%    fileName = '../3DOSU_data_out/fiber_z/fiber_z_sn_normal';
    fileName = '../3DOSU_data_out/fiber_3060/fiber_3060_sn_normal';
end
if nargin < 2
    extension = 'sum';
end
if nargin < 3
    startAngleCols = [5, 12];
end

if nargin < 4
    angleInDeg = 1;
end
fn = [fileName, '.', extension];
dat = load(fn, '-ascii');
sz = length(startAngleCols);
polarSphericalAngles = cell(sz, 1);
for i = 1:sz
    st = startAngleCols(i);
    angles = dat(:, st:st + 1);
    polarSphericalAngles{i} = angles;
end
