function [xAxis,avgDispxx,stdDispxx,ExactDispxx, yAxis,avgDispyy,stdDispyy,ExactDispyy]=PlotDispErr2(coordinatesFEM,USubpb2,x,y,M,N )
  
ErrDisp = zeros(size(coordinatesFEM,1),1); avgDispErr = 0;

load('Sample14Exact.mat');

for i = 1:size(coordinatesFEM,1)
%     ErrDisp(i) = sqrt((-USubpb2(2*i-1) - (1e-4*coordinatesFEM(i,1)^2 + 2e-4*coordinatesFEM(i,1)+2))^2  + ...
%     (USubpb2(2*i) - (2e-5*coordinatesFEM(i,2)^2 + 2e-4*coordinatesFEM(i,2)+4))^2 );

%     ErrDisp(i) = sqrt( (-USubpb2(2*i-1) - sin(0.05)*(coordinatesFEM(i,2)-211) - (cos(0.05)-1)*(coordinatesFEM(i,1)-211))^2  + ...
%     (USubpb2(2*i) - (cos(0.05)-1)*(coordinatesFEM(i,2)-211) + sin(0.05)*(coordinatesFEM(i,1)-211))^2 );

%     ErrDisp(i) = sqrt((-USubpb2(2*i-1) - (2*sin(1/200*coordinatesFEM(i,1))) )^2  + ...
%     (USubpb2(2*i) - (1.5*sin(1/100*coordinatesFEM(i,2))) )^2 );

%     ErrDisp(i) = sqrt((-USubpb2(2*i-1) - (2*sin(1/100*coordinatesFEM(i,1))) )^2  + ...
%         (USubpb2(2*i) - (1.5*sin(1/50*coordinatesFEM(i,2))) )^2 );

%     ErrDisp(i) = sqrt((USubpb2(2*i-1) - (4*sin(1/200*coordinatesFEM(i,1))) )^2  + ...
%     (USubpb2(2*i) - (3*sin(1/100*coordinatesFEM(i,2))) )^2 );

%     ErrDisp(i) = sqrt((USubpb2(2*i-1) )^2  + ...
%     (USubpb2(2*i) - 8*heaviside(coordinatesFEM(i,2)-200) )^2 );

% Sample-1
%     ErrDisp(i) = sqrt((USubpb2(2*i-1) - 1 )^2  + ...
%     ( USubpb2(2*i) - 1 )^2 );

% Sample-14
    ErrDisp(i) = sqrt( (USubpb2(2*i-1) - Sample14Exact(coordinatesFEM(i,1) ,2))^2 + (USubpb2(2*i) - 0)^2 );
    
    avgDispErr = avgDispErr + ErrDisp(i)^2;
end

avgDispErr = sqrt(avgDispErr/size(coordinatesFEM,1))

% figure; surf(reshape(ErrDisp,M,N)'); axis equal; colorbar; view(2); title('||u-u_0||_{L_2}^{1/2}')
% figure; mesh(x,y,reshape(ErrDisp,M,N)); axis tight; set(gca,'fontSize',18); view(-20, 50); title('Displacement Absolute Error')
 
Dispxx = zeros(size(coordinatesFEM,1),1); avgDispxx = 0; stdDispxx = 0; xAxis = 0; ExactDispxx = 0;
for i = 1:size(coordinatesFEM,1)
    Dispxx(i) =  USubpb2(2*i-1)  ;
end
Dispxx = reshape(Dispxx,M,N);
for i = 1:M
    avgDispxx(i) = sum(Dispxx(i,:))/N;
    stdDispxx(i) = std(Dispxx(i,:));
%     ExactDispxx(i) = -(1e-4*x(i,1)^2 + 2e-4*x(i,1)+2) ;
%     ExactDispxx(i) =  -(2*sin(1/200*x(i,1)));
%     ExactDispxx(i) = -(2*sin(1/100*x(i,1)));
%     ExactDispxx(i) = -(4*sin(1/200*x(i,1)));
%     ExactDispxx(i) = 0;
ExactDispxx(i) = Sample14Exact(coordinatesFEM(i,1),2);
 
end
avgDispxx=avgDispxx';
stdDispxx=stdDispxx';
ExactDispxx=ExactDispxx';
xAxis = x(:,1);


Dispyy = zeros(size(coordinatesFEM,1),1); avgDispyy = 0; stdDispyy = 0; yAxis = 0; ExactDispyy = 0;
for i = 1:size(coordinatesFEM,1)
    Dispyy(i) =  USubpb2(2*i)  ;
end
Dispyy = reshape(Dispyy,M,N);
for i = 1:N
    avgDispyy(i) = sum(Dispyy(:,i))/M;
    stdDispyy(i) = std(Dispyy(:,i));
%     ExactDispyy(i) = (2e-5*y(1,i)^2 + 2e-4*y(1,i)+4) ;
%     ExactDispyy(i) = (1.5*sin(1/100*y(1,i)));
%     ExactDispyy(i) = (1.5*sin(1/50*y(1,i)));
%     ExactDispyy(i) = (3*sin(1/100*y(1,i)));
%     ExactDispyy(i) = 8*heaviside(y(1,i)-200);
%     ExactDispyy(i) = 0;
ExactDispyy(i) = 0;

end
avgDispyy=avgDispyy';
stdDispyy=stdDispyy';
ExactDispyy=ExactDispyy';
yAxis = y(1,:)';