function Vq = fungInterpolation_g3(Xq,Yq,Zq,I)

x0 = floor(Xq) - 1;
y0 = floor(Yq) - 1;
z0 = floor(Zq) - 1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I = f(x0:(x0+3), y0:(y0+3), z0:(z0+3));
Dxtemp = [1 1 1 1 1; 
      2 2 2 2 2;
      Xq-x0+1 Xq-x0+1 Xq-x0+1 Xq-x0+1 Xq-x0+1; 
      3 3 3 3 3; 
      4 4 4 4 4];
Dx = zeros(5,5,5);
for tempi = 1:5
    Dx(:,:,tempi) = Dxtemp;
end

Dytemp = [1 2 Yq-y0+1 3 4;
      1 2 Yq-y0+1 3 4;
      1 2 Yq-y0+1 3 4;
      1 2 Yq-y0+1 3 4;
      1 2 Yq-y0+1 3 4];
Dy = zeros(5,5,5);
for tempj = 1:5
    Dy(:,:,tempj) = Dytemp;
end

Dz = zeros(5,5,5);
Dz(:,:,1) = ones(5,5);
Dz(:,:,2) = 2*ones(5,5);
Dz(:,:,3) = (Zq-z0+1)*ones(5,5);
Dz(:,:,4) = 3*ones(5,5);
Dz(:,:,5) = 4*ones(5,5);

I2 = ba_interp3(I,Dy,Dx,Dz,'cubic');
Vq = I2(3,3,3);



%%%%%%%%%% figure;
%%%%%%%%% vi(I2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
