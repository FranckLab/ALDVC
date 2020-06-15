%% ***************************
% FEM Mesh Set Up for DVC
% This code is to set up mesh for both Global FEM DIC and Local Subset DIC
% Author: Jin Yang
% Date: Oct 2017
% ****************************
function [DVCmesh] = MeshSetUp3(xyz,DVCpara)

winstepsize = DVCpara.winstepsize;
x = xyz.x; y = xyz.y; z = xyz.z;

%% ========== mesh for global method ==========
M = size(x,1);  N = size(x,2); L = size(x,3); % N is vertically in image; M is horizontally in image; L is in z direction
coordinatesFEM = zeros(M*N*L,3);
 
% I have transpose x and y because Matlab is read matrix in column direction
% x = premute(x,[2,1,3]); y = premute(y,[2,1,3]); z = premute(z,[2,1,3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xtemp = zeros(M,N,L); ytemp = zeros(M,N,L);
% for tempi = 1:size(x,3), xtemp(:,:,tempi) = x(:,:,tempi)'; ytemp(:,:,tempi) = y(:,:,tempi)'; end
% x = xtemp; y = ytemp;
 
for i = 1:size(coordinatesFEM,1)
    coordinatesFEM(i,:)  = [x(i),y(i),z(i)];
    % x is horizontal position in the image
    % y is vertical position in the image
    % z is z-direction position in the image
end

elementsFEM = zeros((M-1)*(N-1)*(L-1),8);
for k = 1:L-1
    for j = 1:N-1
        for i = 1:M-1
            elementsFEM( (k-1)*(N-1)*(M-1)+(j-1)*(M-1)+i ,:) = ...
                [ (k-1)*N*M+(j-1)*(M)+i (k-1)*N*M+(j-1)*(M)+i+1 (k-1)*N*M+j*(M)+i+1 (k-1)*N*M+j*(M)+i ...
                  (k)*N*M+(j-1)*(M)+i (k)*N*M+(j-1)*(M)+i+1 (k)*N*M+j*(M)+i+1 (k)*N*M+j*(M)+i];
        end
    end
end


%% mesh for local method
% N is vertically in image; M is horizontally in image;

% coordinates = zeros((M+1)*(N+1)*(L+1),3);
% gridxtemp = [x(1,1,1)-0.5*winstepsize; reshape(x(:,1,1),M,1)+0.5*winstepsize];
% gridytemp = [y(1,1,1)-0.5*winstepsize reshape(y(1,:,1),1,N)+0.5*winstepsize]';
% gridztemp = [z(1,1,1)-0.5*winstepsize; reshape(z(1,1,:),L,1)];
% 
% clear gridx gridy gridz
% [gridx,gridy,gridz]=meshgrid(gridxtemp,gridytemp,gridztemp);
% 
% gridxtemp = zeros(size(gridx,2),size(gridx,1),L); 
% gridytemp = zeros(size(gridx,2),size(gridx,1),L);
% for tempi = 1:size(gridx,3)
%     gridxtemp(:,:,tempi) = gridx(:,:,tempi)'; 
%     gridytemp(:,:,tempi) = gridy(:,:,tempi)'; 
% end
% gridx = gridxtemp; gridy = gridytemp;
% 
% for i = 1:size(coordinates,1)
%     coordinates(i,:)  = [gridx(i),gridy(i),gridz(i)];
%     % x is horizontal position in the image
%     % y is vertical position in the image
%     % z is the z direction position in the image
% end
% 
% elements = zeros( M*N*L ,8);
% for k = 1:L
%     for j = 1:N
%         for i = 1:M
%             elements( (k-1)*(N)*(M)+(j-1)*(M)+i ,:) = ...
%             [ (k-1)*(N+1)*(M+1)+(j-1)*(M+1)+i  (k-1)*(N+1)*(M+1)+(j-1)*(M+1)+i+1 ...
%               (k-1)*(N+1)*(M+1)+j*(M+1)+i+1   (k-1)*(N+1)*(M+1)+j*(M+1)+i ...
%               (k)*(N+1)*(M+1)+(j-1)*(M+1)+i  (k)*(N+1)*(M+1)+(j-1)*(M+1)+i+1 ...
%               (k)*(N+1)*(M+1)+j*(M+1)+i+1   (k)*(N+1)*(M+1)+j*(M+1)+i ];
%         end
%     end
% end

% ======== Assign BC values ==========
% % -------- dirichlet BC ----- Need to change---
% % dirichlet = [linspace(1,(M-1),(M-1))', linspace(2,M,(M-1))' ;
% %             linspace(1,1+M*(N-2),(N-1))', linspace((1+M),(1+M*(N-1)),(N-1))' ;
% %             linspace(M,M*(N-1),(N-1))', linspace(2*M,M*N,(N-1))' ;
% %             linspace((M*(N-1)+1), (M*N-1), (M-1))', linspace((M*(N-1)+2), (M*N), (M-1))' ];
% %dirichlet = [];
% % Find all the points with boundary points
% [row1,~] = find(coordinatesFEM(:,1)==min(x(:))); [row2,~] = find(coordinatesFEM(:,1)==max(x(:)));
% [row3,~] = find(coordinatesFEM(:,2)==min(y(:))); [row4,~] = find(coordinatesFEM(:,2)==max(y(:)));
% [row5,~] = find(coordinatesFEM(:,3)==min(z(:))); [row6,~] = find(coordinatesFEM(:,3)==max(z(:)));
% 
% dirichlet = unique([row1;row2;row3;row4;row5;row6]);
dirichlet = [];

% -------- neumann BC ----Need to change ----       
neumann = [linspace(1,(M-1),(M-1))', linspace(2,M,(M-1))', zeros(M-1,1), -ones(M-1,1) ;
             linspace(1,1+M*(N-2),(N-1))', linspace((1+M),(1+M*(N-1)),(N-1))', -ones(N-1,1), zeros(N-1,1) ;
             linspace(M,M*(N-1),(N-1))', linspace(2*M,M*N,(N-1))', ones(N-1,1), zeros(N-1,1) ;
             linspace((M*(N-1)+1), (M*N-1), (M-1))', linspace((M*(N-1)+2), (M*N), (M-1))', zeros(M-1,1), ones(M-1,1) ];
% neumann = [];
 
xyz.x = x; xyz.y = y; xyz.z = z;

%% Assign variables
DVCmesh.coordinatesFEM = coordinatesFEM;
DVCmesh.elementsFEM = elementsFEM;
DVCmesh.dirichlet = dirichlet;
DVCmesh.neumann = neumann;
DVCmesh.xyz0 = xyz;


