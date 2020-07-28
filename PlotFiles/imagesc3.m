function imagesc3(vol)
% IMAGESC3 function is the plot function for ALDVC code
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

imagesc3D(permute(vol,[2,1,3]));
