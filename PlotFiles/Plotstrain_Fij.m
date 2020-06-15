%close all;

clear u_ dm m2vx

u_ = cell(1);
temp1 = reshape(ULocal(1:3:end),M,N,L); temp10 = zeros(N,M,L);
for tempkkk = 1:L
    temp10(:,:,tempkkk) = temp1(:,:,tempkkk)';
end
u_{1}{1} = temp10;
temp1 = reshape(ULocal(2:3:end),M,N,L); temp10 = zeros(N,M,L);
for tempkkk = 1:L
    temp10(:,:,tempkkk) = temp1(:,:,tempkkk)';
end
u_{2}{1} = temp10;
 
temp1 = reshape(ULocal(3:3:end),M,N,L);  temp10 = zeros(N,M,L);
for tempkkk = 1:L
    temp10(:,:,tempkkk) = temp1(:,:,tempkkk)';
end
u_{3}{1} = temp10;
temp1 = reshape(sqrt((ULocal(1:3:end).^2 + ULocal(2:3:end).^2 + ULocal(3:3:end).^2)),M,N,L);  temp10 = zeros(N,M,L);
for tempkkk = 1:L
    temp10(:,:,tempkkk) = temp1(:,:,tempkkk)';
end
u_{4}{1} = temp10;
 
dm  = [winstepsize];

m2vx  = [1,1,1];

[Fij, J] = calculateFij(u_,dm,m2vx,'optimal9'); 

%% Calculate Lagrangian strain (Eij), eulerian strain (eij), Strain energy density (U), and stress (Sij)
[Eij, eij] = calculateEij(Fij);              
%[Sij, U] = calculateSij(Fij, J, eij{1,1}, Eij{1,1}, [E v], mechModel); 
%save('mechanicsvariables_dyedcell_xy5.mat', 'Fij','J','Eij','eij','U','Sij');

%%
FSubpb3 = zeros(M*N*L*9,1);
FSubpb3(1:9:end) = reshape(Fij{1}{5},M*N*L,1)-1;
FSubpb3(5:9:end) = reshape(Fij{1}{1},M*N*L,1)-1;
FSubpb3(9:9:end) = reshape(Fij{1}{9},M*N*L,1)-1;
FSubpb3(2:9:end) = reshape(Fij{1}{4},M*N*L,1);
FSubpb3(4:9:end) = reshape(Fij{1}{2},M*N*L,1);
FSubpb3(3:9:end) = reshape(Fij{1}{6},M*N*L,1);
FSubpb3(7:9:end) = reshape(Fij{1}{8},M*N*L,1);
FSubpb3(6:9:end) = reshape(Fij{1}{3},M*N*L,1);
FSubpb3(8:9:end) = reshape(Fij{1}{7},M*N*L,1);
%for index=[1:9]
%FSubpb3(index:9:end) = reshape(Fij{1}{index},M*N*L,1);
%end 
%Plotstrain_show3(FSubpb3,ResultcoordinatesFEM{1}.coordinatesFEM,ResultelementsFEM{1}.elementsFEM);


%%
Rad = 2; % FilterSizeInput-2;
% Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
temp = 1:1:size(coordinatesFEM,1); temp = temp';
temp = reshape(temp,M,N,L); temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad, Rad+1:L-Rad);
temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad)*(L-2*Rad),1);

temp3 = zeros(4*(M-2*Rad)*(N-2*Rad)*(L-2*Rad),1);
for i = 1:(M-2*Rad)*(N-2*Rad)*(L-2*Rad)
    temp3(9*i-8:9*i) = 9*temp2(i)*ones(9,1)+[-8:1:0]';
end
FStraintemp = FSubpb3(temp3);
  
%[x,y,z] = ndgrid(1:52,1:52,1:19);
%Plotstrain03(FStraintemp,xyz0.x(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad),xyz0.y(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad), ...
%    xyz0.z(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad),size(Img{1}),'Individual');





