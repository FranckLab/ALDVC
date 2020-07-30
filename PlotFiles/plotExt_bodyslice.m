%%%%%%%%%%%%% Extensions for body and slice plottings %%%%%%%%%%%%%%%%
% Extension of ALDVC code: to plot solved displacements and strains
%
%   INPUT: main_ALDVC.m saved matfile after Section 8.
%          Feel free to modify this code based on your purpose. 
%          I comment "#TOMODIFY" where you can easily modify the code.
%
%   OUTPUT: body plots and slice plots for solved displacements and strains
%
%          %%%%%%%%%%% PLOT DISPLACEMENTS %%%%%%%%%%%%%
%          (I-i)   Body view of solved displacements
%          (I-ii)  Slice view of solved displacements
%          (I-iii) Quiver plot of solved displacements
%          (I-iv)  Cone plot of solved displacements
%          (I-v)   Streamline plot of solvd displacements
%
%          %%%%%%%%%%% PLOT STRAINS %%%%%%%%%%%%%
%          (II-i)  Body view of solved strains
%          (II-ii) Slice view of solved strains
%          
%
% ****** ATTENTION ******
% The "x,y,z" or "1-,2-,3-" coordinates in this exchange file correspond to 
% the 1st, 2nd and 3rd indices of Matlab workspace variable. For example, 
% p_meas(:,1) and p_meas(:,2) are the x- & y-coordinates of scattered points. 
%
% This is a little different from some MATLAB image processing functions. 
% For example, if a 3D image has size M*N*L, in this code, we always have
% the image size_x=M, size_y=N, size_z=L. If you use some Matlab computer
% vision/image post-processing function, for example, 'imagesc3D', or
% 'imagesc', or 'imshow', or 'surf', it will reads size_x=N, size_y=M, size_z=L.
%
% Please pay attention to this.  
% 
%
% -------------------------------------------
% Author: Jin Yang, aldicdvc@gmail.com
% Date: 2020.07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Initialization
addpath('./func','./PlotFiles');
% If your current Folder is inside "./3D_ALDVC/PlotFiles", please use:
% addpath('./func','../PlotFiles');

% You can directly start these plottings from saved results matfile
% #TOMODIFY:
ImgSeqNum = 2; % Define frame sequence number, starting from 2, the first ref image always has "ImgSeqNum=1"


%%%%%%%%%%% PLOT DISPLACEMENTS %%%%%%%%%%%%%
%% (I-i) Body view of solved displacements
Plotdisp03(ResultDisp{ImgSeqNum-1}.U,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM,DVCpara.PlotComponentEachOrAll);
% #TOMODIFY:
% Change "DVCpara.PlotComponentEachOrAll" to "Individual" to plot each displacement component in a single figure 
Plotdisp03(ResultDisp{ImgSeqNum-1}.U,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM,'Individual');


%% (I-ii) Slice view of solved displacements
% To view slice of solved displacements
xyz0 = DVCmesh.xyz0; MNL = size(xyz0.x); % xyz0 is the xyz coordinates of nodal points
U = ResultDisp{ImgSeqNum-1}.U; % U is a long vector for measured U = [Ux_pt1, Uy_pt1, Uz_pt1,  Ux_pt2, Uy_pt2, Uz_pt2, ... ]';
Ux = reshape(U(1:3:end),MNL); Uy = reshape(U(2:3:end),MNL); Uz = reshape(U(3:3:end),MNL);
% Ux,Uy,Uz are x-, y-, and z- components of solved displacements, whose sizes are same with xyz0.x or xyz0.y or xyz0.z

figure, imagesc3(Ux); % The units in this plot is the neighboring subset distance "winstepsize" 
% Coordinates in this plot: going right is x+ direction, going downwards is y+ direction

% #TOMODIFY:
sliceNum = 6; % z-slice #

figure, surf(squeeze(xyz0.x(:,:,sliceNum)),squeeze(xyz0.y(:,:,sliceNum)),squeeze(Ux(:,:,sliceNum))); 
view([0,-90]); axis equal; axis tight; set(gca,'fontsize',18); cb = colorbar;
title(['Solved x-disp at z-slice #',num2str(sliceNum),'  (unit:vx)'],'fontweight','normal');


%% (I-iii) Quiver plot of solved displacements
xyz0 = DVCmesh.xyz0; MNL = size(xyz0.x); % xyz0 is the xyz coordinates of nodal points
uxGrid = reshape(U(1:3:end),MNL); uyGrid = reshape(U(2:3:end),MNL);  uzGrid = reshape(U(3:3:end),MNL);
figure, quiver3(xyz0.x, xyz0.y, xyz0.z, uxGrid, uyGrid, uzGrid);
xlabel('x axis'); ylabel('y axis'); zlabel('z axis'); set(gca,'fontsize',16); 


%% (I-iv) Cone plot of solved displacements
plotCone3(DVCmesh.coordinatesFEM, ResultDisp{ImgSeqNum-1}.U);


%% (I-v) Streamline plot of solvd displacements
xyz0 = DVCmesh.xyz0; MNL = size(xyz0.x); % xyz0 is the xyz coordinates of nodal points
uxGrid = reshape(U(1:3:end),MNL); uyGrid = reshape(U(2:3:end),MNL);  uzGrid = reshape(U(3:3:end),MNL);

% Define streamline starting points
[yGridSlice,xGridSlice,zGridSlice] = meshgrid( xyz0.y(1):30:xyz0.y(end), xyz0.x(1):30:xyz0.x(end), xyz0.z(1):10:xyz0.z(end));

% Plot streamline
plotStreamline3( xyz0.y, xyz0.x, xyz0.z, uxGrid, uyGrid, uzGrid, yGridSlice,xGridSlice,zGridSlice);
 


%%%%%%%%%%% PLOT STRAINS %%%%%%%%%%%%%
%% (ii-i) Body view of solved strains
xyz0 = DVCmesh.xyz0; % xyz0 is the xyz coordinates of nodal points
F = ResultDefGrad{ImgSeqNum-1}.F; % F is a long vector for solved F = (grad_X x - I) =  ...
% [F11_pt1, F21_pt1, F31_pt1, F12_pt1, F22_pt1, F32_pt1, F13_pt1, F23_pt1, F33_pt1, ...
% [F11_pt2, F21_pt2, F31_pt2, F12_pt2, F22_pt2, F32_pt2, F13_pt2, F23_pt2, F33_pt2, ... ]';
Plotstrain03(full(F),xyz0.x, xyz0.y, xyz0.z,0,DVCpara.PlotComponentEachOrAll);
% #TOMODIFY:
% Change "DVCpara.PlotComponentEachOrAll" to "Individual" to plot each strain component in a single figure 
% Plotstrain03(full(F),xyz0.x, xyz0.y, xyz0.z,0,'Individual');

% ===================================
% If you want to plot strains after executing main_ALDVC.m Section 8 (where additional modifications can be done)
xyz0 = DVCmesh.xyz0; MNL = size(xyz0.x);
FStraintemp = ResultStrain{ImgSeqNum-1}.Strain; coordinatesFEMStrain = ResultStrain{ImgSeqNum-1}.coordinatesFEMStrain;
tempRad=0; 
while abs(prod(MNL-tempRad)-size(coordinatesFEMStrain,1))>0 && tempRad<size(coordinatesFEMStrain,1)
    tempRad = tempRad+2;
end
Rad=repmat(tempRad/2,1,3); % Or manually assign the value of "Rad" used in previous Section 8
 
% each displacement component in a single figure
Plotstrain03(full(FStraintemp),xyz0.x(1+Rad(1):MNL(1)-Rad(1),1+Rad(2):MNL(2)-Rad(2),1+Rad(3):MNL(3)-Rad(3)), ...
    xyz0.y(1+Rad(1):MNL(1)-Rad(1),1+Rad(2):MNL(2)-Rad(2),1+Rad(3):MNL(3)-Rad(3)), ...
    xyz0.z(1+Rad(1):MNL(1)-Rad(1),1+Rad(2):MNL(2)-Rad(2),1+Rad(3):MNL(3)-Rad(3)),0,DVCpara.PlotComponentEachOrAll);
% #TOMODIFY:
% Change "DVCpara.PlotComponentEachOrAll" to "Individual" to plot each strain component in a single figure
% Plotstrain03(full(FStraintemp),xyz0.x(1+Rad(1):MNL(1)-Rad(1),1+Rad(2):MNL(2)-Rad(2),1+Rad(3):MNL(3)-Rad(3)), ...
%     xyz0.y(1+Rad(1):MNL(1)-Rad(1),1+Rad(2):MNL(2)-Rad(2),1+Rad(3):MNL(3)-Rad(3)), ...
%     xyz0.z(1+Rad(1):MNL(1)-Rad(1),1+Rad(2):MNL(2)-Rad(2),1+Rad(3):MNL(3)-Rad(3)),0,'Individual');
  

%% (II-ii) Slice view of solved strains
xyz0 = DVCmesh.xyz0; MNL = size(xyz0.x); % xyz0 is the xyz coordinates of nodal points
F = ResultDefGrad{ImgSeqNum-1}.F; % F is a long vector for solved F = (grad_X x - I) =  ...
% [F11_pt1, F21_pt1, F31_pt1, F12_pt1, F22_pt1, F32_pt1, F13_pt1, F23_pt1, F33_pt1, ...
% [F11_pt2, F21_pt2, F31_pt2, F12_pt2, F22_pt2, F32_pt2, F13_pt2, F23_pt2, F33_pt2, ... ]';
 
F11 = reshape(F(1:9:end),MNL); F21 = reshape(F(2:9:end),MNL); F31 = reshape(F(3:9:end),MNL);
F12 = reshape(F(4:9:end),MNL); F22 = reshape(F(5:9:end),MNL); F32 = reshape(F(6:9:end),MNL);
F13 = reshape(F(7:9:end),MNL); F23 = reshape(F(8:9:end),MNL); F33 = reshape(F(9:9:end),MNL);
% Individual components of solved strains, whose sizes are same with xyz0.x or xyz0.y or xyz0.z

figure, imagesc3(F11); 

% #TOMODIFY:
sliceNum = 6; % z-slice #

figure, surf(squeeze(xyz0.x(:,:,sliceNum)),squeeze(xyz0.y(:,:,sliceNum)),squeeze(F11(:,:,sliceNum))); 
view([0,-90]); axis equal; axis tight; set(gca,'fontsize',18); cb = colorbar;
title(['Strain exx at z-slice #',num2str(sliceNum)],'fontweight','normal');


% ===================================
% If you want to plot strains after executing main_ALDVC.m Section 8 (where additional modifications can be done)
FStraintemp = ResultStrain{ImgSeqNum-1}.Strain; coordinatesFEMStrain = ResultStrain{ImgSeqNum-1}.coordinatesFEMStrain;
xList=unique(coordinatesFEMStrain(:,1)); yList=unique(coordinatesFEMStrain(:,2)); zList=unique(coordinatesFEMStrain(:,3));
MNLtemp=[length(xList),length(yList),length(zList)]; [yGrid,xGrid]=meshgrid(yList,xList);
 
F11 = reshape(FStraintemp(1:9:end),MNLtemp); F21 = reshape(FStraintemp(2:9:end),MNLtemp); F31 = reshape(FStraintemp(3:9:end),MNLtemp);
F12 = reshape(FStraintemp(4:9:end),MNLtemp); F22 = reshape(FStraintemp(5:9:end),MNLtemp); F32 = reshape(FStraintemp(6:9:end),MNLtemp);
F13 = reshape(FStraintemp(7:9:end),MNLtemp); F23 = reshape(FStraintemp(8:9:end),MNLtemp); F33 = reshape(FStraintemp(9:9:end),MNLtemp);

figure, imagesc3(F11); 

% #TOMODIFY:
sliceNum = 6; % z-slice #
figure, surf(xGrid,yGrid,squeeze(F11(:,:,sliceNum))); 
view([0,-90]); axis equal; axis tight; set(gca,'fontsize',18); cb = colorbar;
title(['Strain exx at z-slice #',num2str(sliceNum)],'fontweight','normal');





