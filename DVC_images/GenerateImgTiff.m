% Generate image tiffs from matfile
%
% === ATTENTION ===
% Please pay attention to the discrepancy between image intrinsic coordinates 
% and the Matlab variable cooridinates we are using in ALDVC code.
%
% For example, each image stack has the following xy coordinates
%   -------- x plus direction ------->
%   |
%   |
%   |y
%   |plus
%   |direction
%   |
%   |
%   V
% So if you use imagesc3D (similar to imagesc) function to directly plot 
% these volumetric images, The first index of the image corresponds to the 
% y axis, and the second index corresponds to the x axis.
%
% However, in all the variables used in ALDVC code, the 1st, 2nd, 3rd indices
% are corresponding to the x, y, and z coordinates, respectively. 
% For example, Img{1} in the ALDVC code is a 3-dimensional data variable, the
% 1st, 2nd, 3rd indicies are corresponding to x, y, and z coordinates;
% Another example, in ALDVC, coordinatesFEM is a matrix with size [MNL,3],
% coordinatesFEM(:,1), coordinatesFEM(:,2), coordinatesFEM(:,3) are x, y,
% and z coordinates of FE-mesh nodes.
%
% So in all the ALDVC code, we use imagesc3 function to flip the x and y axes. 
% If you want to use imagesc3D, please use permute(xxx, [2,1,3]) to switch 
% the x and y axes. If you just follow the code, please ignore these comments.
% 
%
% For more information please contact: Jin Yang, Email: aldicdvc@gmail.com 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Load matlab vol matfiles
vol1 = load('vol_stretch_1001.mat'); vol1=permute(vol1.vol{1},[2,1,3]); 
vol2 = load('vol_stretch_1002.mat'); vol2=permute(vol2.vol{1},[2,1,3]);

%% Generate ref tiffs
cd('./vol_stretch_1001_tiff'); 
fileName = 'vol_stretch_1001';

for tempi=1:size(vol1,3)
    
    f=(squeeze(vol1(:,:,tempi)));  
    if tempi<10
        frameName = [fileName,'_000',num2str(tempi),'.tif'];
    elseif tempi<100
        frameName = [fileName,'_00',num2str(tempi),'.tif'];
    elseif tempi<1000
        frameName = [fileName,'_0',num2str(tempi),'.tif'];
    elseif tempi<10000
        frameName = [fileName,'_',num2str(tempi),'.tif'];
    else
        disp('Please modify GenerateImgTiff.m'); pause;
    end
    imwrite(f,frameName);
    
end
cd('../');

%% Generate def tiffs
cd('./vol_stretch_1002_tiff'); 
fileName = 'vol_stretch_1002';

for tempi=1:size(vol1,3)
    
    f=(squeeze(vol2(:,:,tempi)));  
    if tempi<10
        frameName = [fileName,'_000',num2str(tempi),'.tif'];
    elseif tempi<100
        frameName = [fileName,'_00',num2str(tempi),'.tif'];
    elseif tempi<1000
        frameName = [fileName,'_0',num2str(tempi),'.tif'];
    elseif tempi<10000
        frameName = [fileName,'_',num2str(tempi),'.tif'];
    else
        disp('Please modify GenerateImgTiff.m'); pause;
    end
    imwrite(f,frameName);
    
end
cd('../');