%% Correct and Smooth initial displacements
function [U0] =Init3(uvw,xyz)

% in inpaint_nans part, I use Spring Model for u and v initial guesses
uInit = uvw.u; vInit = uvw.v; wInit = uvw.w;
x = xyz.x; y = xyz.y; z = xyz.z;

%% Initial Value
U0 = zeros(3*size(x,1)*size(y,2)*size(z,3),1); 

% Here transpose bcz following give valus in column order
% M = size(x,1); N = size(y,2); L = size(z,3);
% uInittemp = zeros(M,N,L); vInittemp = zeros(M,N,L); wInittemp = zeros(M,N,L); % Phitemp = zeros(M,N,L);
% for tempi = 1:size(x,3)
%     uInittemp(:,:,tempi) = uInit(:,:,tempi)'; 
%     vInittemp(:,:,tempi) = vInit(:,:,tempi)'; 
%     wInittemp(:,:,tempi) = wInit(:,:,tempi)'; 
% end
 

for i = 1:(size(x,1)*size(y,2)*size(z,3))
    U0(3*i-2) = uInit(i); %u(i);
    U0(3*i-1) = vInit(i); %v(i);
    U0(3*i)   = wInit(i); %w(i);
end
disp('--- Finish setting up mesh and assigning initial value! ---')
 

end