%% Contour plot
clear;

%% Load Cell border data
load('cell001xy.mat')
im = squeeze(sum(BW,3));
im = im/max(im(:));

% load('cellxzyz.mat')
% im = squeeze(sum(BW,2));

%% Format Displacements

load('resultsFIDVC.mat')
timePt = 1;
slicePlane = 30;
sizeI = [512,512,64];
dm = 1;
n = cell(3,1);
[n{:}] = meshgrid(1:dm:sizeI(1),1:dm:sizeI(2),slicePlane);

for i = 1:3
    u{timePt}{i} = permute(u{timePt}{i},[2,1,3]);
    
    driftRegion = u{timePt}{i}(:);
    u{timePt}{i} = u{timePt}{i} - median(driftRegion(:));

    F = scatteredInterpolant(m{1}(:),m{2}(:),m{3}(:),u{timePt}{i}(:),'natural','none');
    U{i} = F(n{1},n{2},n{3});
end
U{4} = sqrt(U{1}.^2 + U{2}.^2 + U{3}.^2);
% umag_max = max(U{4}(:));


%% Format Tractions





%% Cones for plotting

% Random position for cone plot
bw = ones(size(BW,1),size(BW,2));
bw = logical(bw);
% Density of cone plot 
density = 6.5;
coneSize = 0.025;
density = 1/density^2;

% Density_mask;
cone_mask = rand(size(bw));
cone_mask(cone_mask>density) = 0;
cone_mask(cone_mask>0) = 1;

% Finding index of cone position;
idx = find(cone_mask);

[vert{2},vert{1},vert{3}] = ind2sub(size(cone_mask),idx);
% U{4} = sqrt(U{1}.^2 + U{2}.^2 + U{3}.^2);


% Calibration for x,y,z directions
v{1}{1} = 0.155*U{2}(idx);
v{1}{2} = 0.155*U{1}(idx);
v{1}{3} = 0.3*U{3}(idx);
mag = sqrt(v{1}{1}.^2 + v{1}{2}.^2 + v{1}{3}.^2);
umag_max = max(mag);
%% Colormap setting
figure;
tt = zeros(128,3);
ttt = [179,136,255];
for i = 1:128
    tt(i,:) = [0,i/128,0];
    tt(i,:) = ttt/255*i/128/3+ttt/255/3*2;
end
tt(1,:) = [0,0,0];
for i = 2:40
    tt(i,:) = tt(41,:)*i/41;
end

colormap([tt;jet(128)]);
IM = umag_max*double(im)-umag_max;

[~,hC]= contourf((IM),128);
set(hC,'LineStyle','none');
axis image
freezeColors
hold on

hc = coneplot(vert{1},vert{2},vert{3},v{1}{2}(),v{1}{1},v{1}{3},coneSize,'nointerp'); %Size 0.01 for zoom and 0.015 for full size
caxis([-umag_max,umag_max])
fvc = repmat(mag(:)',[42 1]);
set(hc, 'facecolor', 'flat', 'facevertexcdata', fvc(:));
hc.EdgeColor = 'none';
hc.AmbientStrength = 0.6;
hc.DiffuseStrength = 0.75;
hc.SpecularStrength = 0.4;
hc.SpecularExponent = 3;
%whole
%xlim([100,400])
%ylim([100,400])
%colorbar
axis off;
lighting phong;
export_fig('Soft_1.7_Disp.png',  '-png', '-opengl', '-r1200');



