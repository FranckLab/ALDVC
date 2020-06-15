function [xAxis,avgStrainxx,stdStrainxx,ExactStrainxx, yAxis,avgStrainyy,stdStrainyy,ExactStrainyy]= PlotStrainErr2(coordinatesFEM,FSubpb2,x,y,M,N, Exact)
 
load('Sample14Exact.mat');

ErrStrain = zeros(size(coordinatesFEM,1),1); avgErrStrain = 0;
for i = 1:size(coordinatesFEM,1)
%     ErrStrain(i) = sqrt( (-FSubpb2(4*i-3) - (2e-4*coordinatesFEM(i,1) + 2e-4 ) )^2  + ...
%         (FSubpb2(4*i) - (4e-5*coordinatesFEM(i,2) + 2e-4))^2 + (FSubpb2(4*i-2) )^2 + (FSubpb2(4*i-1) )^2  );
% 
%     ErrStrain(i) = sqrt( (-FSubpb2(4*i-3) - (cos(0.05)-1))^2  + ...
%         (FSubpb2(4*i) - (cos(0.05)-1))^2 + (FSubpb2(4*i-2)+sin(0.05))^2 + (-FSubpb2(4*i-1)-sin(0.05) )^2  );
% 
%     ErrStrain(i) = sqrt((-FSubpb2(4*i-3) - (2/200*cos(1/200*coordinatesFEM(i,1))) )^2  + ...
%         (FSubpb2(4*i) - (1.5/100*cos(1/100*coordinatesFEM(i,2))) )^2 + (FSubpb2(4*i-2) )^2 + (FSubpb2(4*i-1) )^2 );
%
%     ErrStrain(i) = sqrt((-FSubpb2(4*i-3) - (2/100*cos(1/100*coordinatesFEM(i,1))) )^2  + ...
%        (FSubpb2(4*i) - (1.5/50*cos(1/50*coordinatesFEM(i,2))) )^2 + (FSubpb2(4*i-2) )^2 + (FSubpb2(4*i-1) )^2 );
%
%     ErrStrain(i) = sqrt((FSubpb2(4*i-3) - (4/200*cos(1/200*coordinatesFEM(i,1))) )^2  + ...
%        (FSubpb2(4*i) - (3/100*cos(1/100*coordinatesFEM(i,2))) )^2 + (FSubpb2(4*i-2) )^2 + (FSubpb2(4*i-1) )^2 );
% 
%     ErrStrain(i) = sqrt((FSubpb2(4*i-3) )^2  + ...
%       (FSubpb2(4*i) )^2 + (FSubpb2(4*i-2) )^2 + (FSubpb2(4*i-1) )^2 );

% Sample-14
	ErrStrain(i) = sqrt( (FSubpb2(4*i-3)-Sample14Exact(coordinatesFEM(i,1),3)*1e-6)^2 + (FSubpb2(4*i))^2 + (FSubpb2(4*i-2))^2 + (FSubpb2(4*i-1))^2 );
   
	avgErrStrain = avgErrStrain + ErrStrain(i)^2; 
end
avgErrStrain = sqrt(avgErrStrain/size(coordinatesFEM,1))
% figure; surf(reshape(ErrDisp,M,N)'); axis equal; colorbar; view(2); title('||u-u_0||_{L_2}^{1/2}')
% figure; mesh(x,y,reshape(ErrStrain,M,N)); axis tight; set(gca,'fontSize',20); view(-20, 50); title('Strain Absolute Error')
  
  Strainxx = zeros(size(coordinatesFEM,1),1); avgStrainxx = 0; stdStrainxx = 0; xAxis = 0; ExactStrainxx = 0;
  for i = 1:size(coordinatesFEM,1)
  Strainxx(i) =  FSubpb2(4*i-3)  ;
  end
  Strainxx = reshape(Strainxx,M,N); 
  for i = 1:M
      avgStrainxx(i) = sum(Strainxx(i,:))/N;
      stdStrainxx(i) = std(Strainxx(i,:));
%     ExactStrainxx(i) = -(2e-4*x(i,1)  + 2e-4 ) ;
%     ExactStrainxx(i) = -(2/200*cos(1/200*x(i,1)));
%     ExactStrainxx(i) = -(2/100*cos(1/100*x(i,1)));
%     ExactStrainxx(i) = (4/200*cos(1/200*x(i,1)));
%     ExactStrainxx(i) = 0;
	  ExactStrainxx(i) = Sample14Exact(coordinatesFEM(i,1),3)*1e-6;
  end
  avgStrainxx=avgStrainxx';
  stdStrainxx=stdStrainxx';
    ExactStrainxx=ExactStrainxx';
    xAxis = x(:,1);  
  %figure; errorbar(xAxis, avgStrainxx  , stdStrainxx );
 
Strainyy = zeros(size(coordinatesFEM,1),1); avgStrainyy = 0; stdStrainyy = 0; yAxis = 0; ExactStrainyy = 0;
for i = 1:size(coordinatesFEM,1)
    Strainyy(i) =  FSubpb2(4*i)  ;
end
Strainyy = reshape(Strainyy,M,N);
for i = 1:N
    avgStrainyy(i) = sum(Strainyy(:,i))/M;
    stdStrainyy(i) = std(Strainyy(:,i));
%     ExactStrainyy(i) = (4e-5*y(1,i) + 2e-4 ) ;
%     ExactStrainyy(i) =  (1.5/100*cos(1/100*y(1,i)));
%     ExactStrainyy(i) = (1.5/50*cos(1/50*y(1,i)));
%     ExactStrainyy(i) = (3/100*cos(1/100*y(1,i)));
	ExactStrainyy(i) =0;
end
avgStrainyy=avgStrainyy';
stdStrainyy=stdStrainyy';
ExactStrainyy=ExactStrainyy';
yAxis = y(1,:)';