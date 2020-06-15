function [xAxis,avgStrainErr,avgStrainxx,stdStrainxx,ExactStrainxx] = PlotStrainErr3( ...
                                coordinatesFEM,FSubpb2,winstepsize,ImgSeqNum)
 
M = 1 + (coordinatesFEM(end,1)-coordinatesFEM(1,1))/winstepsize(1);
N = 1 + (coordinatesFEM(end,2)-coordinatesFEM(1,2))/winstepsize(2);
L = 1 + (coordinatesFEM(end,3)-coordinatesFEM(1,3))/winstepsize(3);

ErrStrain = zeros(size(coordinatesFEM,1),1); avgErrStrain = 0;

load('Sample14Exact.mat');

for i = 1:size(coordinatesFEM,1)
%     ErrStrain(i) = sqrt( (-FSubpb2(4*i-3) - (2e-4*coordinatesFEM(i,1) + 2e-4 ) )^2  + ...
%         (FSubpb2(4*i) - (4e-5*coordinatesFEM(i,2) + 2e-4))^2 + (FSubpb2(4*i-2) )^2 + (FSubpb2(4*i-1) )^2  );

% Sample-14
	ErrStrain(i) = sqrt( (-FSubpb2(9*i-8) - Sample14Exact(coordinatesFEM(i,1),2*ImgSeqNum-1)*1e-6)^2 + ...
        (FSubpb2(9*i-7))^2 + (FSubpb2(9*i-6))^2 + (FSubpb2(9*i-5))^2 + (FSubpb2(9*i-4))^2 + ...
        (FSubpb2(9*i-3))^2 + (FSubpb2(9*i-2))^2 + (FSubpb2(9*i-1))^2 + (FSubpb2(9*i-0))^2 );
   
	avgErrStrain = avgErrStrain + ErrStrain(i)^2; 
end
avgStrainErr = sqrt(avgErrStrain/size(coordinatesFEM,1));
% figure; surf(reshape(ErrDisp,M,N)'); axis equal; colorbar; view(2); title('||u-u_0||_{L_2}^{1/2}')
% figure; mesh(x,y,reshape(ErrStrain,M,N)); axis tight; set(gca,'fontSize',20); view(-20, 50); title('Strain Absolute Error')
  
  Strainxx = zeros(size(coordinatesFEM,1),1); avgStrainxx = 0; stdStrainxx = 0; xAxis = 0; ExactStrainxx = 0;
  for i = 1:size(coordinatesFEM,1)
  Strainxx(i) =  FSubpb2(9*i-8)  ;
  end
  Strainxx = reshape(Strainxx,M,N,L); 
  for i = 1:M
      avgStrainxx(i) = sum(Strainxx(i,:))/(N*L);
      stdStrainxx(i) = std(Strainxx(i,:));
%     ExactStrainxx(i) = -(2e-4*x(i,1)  + 2e-4 ) ;
%     ExactStrainxx(i) = -(2/200*cos(1/200*x(i,1)));
%     ExactStrainxx(i) = -(2/100*cos(1/100*x(i,1)));
%     ExactStrainxx(i) = (4/200*cos(1/200*x(i,1)));
%     ExactStrainxx(i) = 0;
	  ExactStrainxx(i) = Sample14Exact(coordinatesFEM(i,1),2*ImgSeqNum-1)*1e-6;
  end
  avgStrainxx=avgStrainxx';
  stdStrainxx=stdStrainxx';
    ExactStrainxx=ExactStrainxx';
    xAxis = [coordinatesFEM(1,1):winstepsize(1):coordinatesFEM(end,1)];%x(:,1);
  %figure; errorbar(xAxis, avgStrainxx  , stdStrainxx );
  