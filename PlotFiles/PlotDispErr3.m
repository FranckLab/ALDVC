function [xAxis,avgDispErr,avgDispxx,stdDispxx,ExactDispxx] = PlotDispErr3( ...
                                coordinatesFEM,USubpb2,winstepsize,ImgSeqNum)

M = 1 + (coordinatesFEM(end,1)-coordinatesFEM(1,1))/winstepsize(1);
N = 1 + (coordinatesFEM(end,2)-coordinatesFEM(1,2))/winstepsize(2);
L = 1 + (coordinatesFEM(end,3)-coordinatesFEM(1,3))/winstepsize(3);
    
ErrDisp = zeros(size(coordinatesFEM,1),1); avgDispErr = 0; % Initialize

load('Sample14Exact.mat');

for i = 1:size(coordinatesFEM,1)
%     ErrDisp(i) = sqrt((-USubpb2(2*i-1) - (1e-4*coordinatesFEM(i,1)^2 + 2e-4*coordinatesFEM(i,1)+2))^2  + ...
%     (USubpb2(2*i) - (2e-5*coordinatesFEM(i,2)^2 + 2e-4*coordinatesFEM(i,2)+4))^2 );

%     ErrDisp(i) = sqrt((USubpb2(2*i-1) )^2  + ...
%     (USubpb2(2*i) - 8*heaviside(coordinatesFEM(i,2)-200) )^2 );

% Sample-1
%     ErrDisp(i) = sqrt((USubpb2(2*i-1) - 1 )^2  + ...
%     ( USubpb2(2*i) - 1 )^2 );

% Sample-14
    ErrDisp(i) = sqrt( (-USubpb2(3*i-2) - Sample14Exact(coordinatesFEM(i,1),2*ImgSeqNum-2))^2 + ...
                            (USubpb2(3*i-1) - 0)^2 + (USubpb2(3*i-0) - 0)^2 );
    
    avgDispErr = avgDispErr + ErrDisp(i)^2;
    
end

avgDispErr = sqrt(avgDispErr/size(coordinatesFEM,1));

% figure; surf(reshape(ErrDisp,M,N)'); axis equal; colorbar; view(2); title('||u-u_0||_{L_2}^{1/2}')
% figure; mesh(x,y,reshape(ErrDisp,M,N)); axis tight; set(gca,'fontSize',18); view(-20, 50); title('Displacement Absolute Error')
 
Dispxx = zeros(size(coordinatesFEM,1),1); avgDispxx = 0; stdDispxx = 0; xAxis = 0; ExactDispxx = 0;
for i = 1:size(coordinatesFEM,1)
    Dispxx(i) =  USubpb2(3*i-2)  ;
end
Dispxx = reshape(Dispxx,M,N,L);
for i = 1:M
    avgDispxx(i) = sum(Dispxx(i,:))/(N*L);
    stdDispxx(i) = std(Dispxx(i,:));
    % ExactDispxx(i) = 0;
    ExactDispxx(i) = Sample14Exact(coordinatesFEM(i,1),2*ImgSeqNum-2);
end
avgDispxx=avgDispxx';
stdDispxx=stdDispxx';
ExactDispxx=ExactDispxx';
xAxis = [coordinatesFEM(1,1):winstepsize(1):coordinatesFEM(end,1)];%x(:,1);

 