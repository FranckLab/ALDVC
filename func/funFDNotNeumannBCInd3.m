% Find indices of finite difference scheme Neumann boundary points 
function [notNeumannBCInd_F, notNeumannBCInd_u]  = funFDNotNeumannBCInd3(size1coordinatesFEM,MNL,Rad)

if length(Rad)==1, Rad=Rad*ones(1,3); end

temp = 1:1:size1coordinatesFEM; % temp = temp'; 
temp = reshape(temp,MNL); 
temp2 = temp(Rad(1)+1:MNL(1)-Rad(1), Rad(2)+1:MNL(2)-Rad(2), Rad(3)+1:MNL(3)-Rad(3)); 
temp2 = temp2(:); % reshape(temp2, (MNL(1)-2*Rad)*(MNL(2)-2*Rad)*(MNL(3)-2*Rad),1);
temp3 = zeros(9*(MNL(1)-2*Rad(1))*(MNL(2)-2*Rad(2))*(MNL(3)-2*Rad(3)),1);

for i = 1:(MNL(1)-2*Rad(1))*(MNL(2)-2*Rad(2))*(MNL(3)-2*Rad(3)), 
    temp3(9*i-8:9*i) = 9*temp2(i)*ones(1,9)+[-8:1:0]; 
end
% atemp = a(temp3);

temp4 = zeros(3*(MNL(1)-2*Rad(1))*(MNL(2)-2*Rad(2))*(MNL(3)-2*Rad(3)),1);
for i = 1:(MNL(1)-2*Rad(1))*(MNL(2)-2*Rad(2))*(MNL(3)-2*Rad(3)), 
    temp4(3*i-2:3*i) = 3*temp2(i)*ones(1,3)+[-2:1:0]; 
end
% btemp = b(temp4);

notNeumannBCInd_F = temp3(:);
notNeumannBCInd_u = temp4(:);