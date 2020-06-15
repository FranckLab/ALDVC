function [U, F, stepwithinwhile, HGlobal] = funICGN3(U0,Coords,Df,ImgEle,Img2,winsize,tol,method,interpmethod,MaxIterNum)

warning('off','all');

DIM = 3; x0 = Coords(1); y0 = Coords(2); z0 = Coords(3);
% imgfNormalizedbc = Img{1}; imggNormalizedbc = Img{2};
imgfNormalizedbc = ImgEle.Imgf(4:end-3, 4:end-3, 4:end-3);
imggNormalizedbc = Img2; %imggNormalizedbc = ImgEle.Imgg;

% DfDxStartx = Df.DfAxis(1); DfDxStarty = Df.DfAxis(3); DfDxStartz = Df.DfAxis(5);
% imgSize = Df.imgSize;
%ImggStartx = ImgEle.ImggAxis(1); ImggStarty = ImgEle.ImggAxis(3); ImggStartz = ImgEle.ImggAxis(5);
ImggStartx = 0; ImggStarty = 0; ImggStartz = 0;
imgSize = Df.imgSize;

% ---------------------------
% Find local subset region 
x = [x0-floor(winsize(1)/2); x0+floor(winsize(1)/2); x0+floor(winsize(1)/2); x0-floor(winsize(1)/2); 
    x0-floor(winsize(1)/2); x0+floor(winsize(1)/2); x0+floor(winsize(1)/2); x0-floor(winsize(1)/2)]; % [coordinates(elements(j,:),1)];
y = [y0-floor(winsize(2)/2); y0-floor(winsize(2)/2); y0+floor(winsize(2)/2); y0+floor(winsize(2)/2); 
    y0-floor(winsize(2)/2); y0-floor(winsize(2)/2); y0+floor(winsize(2)/2); y0+floor(winsize(2)/2)];  % [coordinates(elements(j,:),2)];
z = [z0-floor(winsize(3)/2); z0-floor(winsize(3)/2); z0-floor(winsize(3)/2); z0-floor(winsize(3)/2); 
    z0+floor(winsize(3)/2); z0+floor(winsize(3)/2); z0+floor(winsize(3)/2); z0+floor(winsize(3)/2)];    % [coordinates(elements(j,:),3)];

% ---------------------------
% Initialization: Get P0
P0 = [zeros(1,9), U0(1) U0(2) U0(3)]'; P = P0;  

% ---------------------------
% Find region for f
[XX,YY,ZZ] = ndgrid([x(1):1:x(7)],[y(1):1:y(7)],[z(1):1:z(7)]);

%tempf = imgfNormalizedbc.eval(XX,YY);
%DfDx = imgfNormalizedbc.eval_Dx(XX,YY);
%DfDy = imgfNormalizedbc.eval_Dy(XX,YY);

% tempf = imgfNormalizedbc([x(1):1:x(7)],[y(1):1:y(7)],[z(1):1:z(7)]);
% DfDx = Df.DfDx((x(1)-DfDxStartx):1:(x(7)-DfDxStartx), (y(1)-DfDxStarty):1:(y(7)-DfDxStarty), (z(1)-DfDxStartz):1:(z(7)-DfDxStartz)); 
% DfDy = Df.DfDy((x(1)-DfDxStartx):1:(x(7)-DfDxStartx), (y(1)-DfDxStarty):1:(y(7)-DfDxStarty), (z(1)-DfDxStartz):1:(z(7)-DfDxStartz)); 
% DfDz = Df.DfDz((x(1)-DfDxStartx):1:(x(7)-DfDxStartx), (y(1)-DfDxStarty):1:(y(7)-DfDxStarty), (z(1)-DfDxStartz):1:(z(7)-DfDxStartz)); 
tempf = imgfNormalizedbc;
% DfDx = DfEle.DfDx; DfDy = DfEle.DfDy; DfDz = DfEle.DfDz;
tempDf = funImgGradient3(ImgEle.Imgf,'stencil7'); 
DfDx = tempDf.DfDx; DfDy = tempDf.DfDy; DfDz = tempDf.DfDz;

meanf = mean(tempf(:));
bottomf = sqrt((length(tempf(:))-1)*var(tempf(:)));

% DfDxStartx = x(1)-1; DfDxStarty = y(1)-1;
% tempf = zeros(winsize+1,winsize+1);
% tempf(1:(winsize+1),1:(winsize+1)) = f(x(1):x(3), y(1):y(3));
%[tempCoordy, tempCoordx] = meshgrid(y(1):y(3),x(1):x(3));
%tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
 
 
H2 = zeros(12,12); DfDxSq = (DfDx.^2); DfDySq = (DfDy.^2); DfDzSq = (DfDz.^2); 
DfDxDfDy = DfDx.*DfDy; DfDxDfDz = DfDx.*DfDz;  DfDyDfDz = DfDy.*DfDz;
XXSq = (XX-x0).^2; YYSq = (YY-y0).^2; ZZSq = (ZZ-z0).^2; 
XXYY = (XX-x0).*(YY-y0); XXZZ = (XX-x0).*(ZZ-z0); YYZZ = (YY-y0).*(ZZ-z0);
H2(1,1) = sum(sum(sum(XXSq.*DfDxSq)));       H2(1,2) = sum(sum(sum(XXSq.*DfDxDfDy)));         H2(1,3) = sum(sum(sum(XXSq.*DfDxDfDz)));
H2(1,4) = sum(sum(sum(DfDxSq.*XXYY)));       H2(1,5) = sum(sum(sum(DfDxDfDy.*XXYY)));         H2(1,6) = sum(sum(sum(DfDxDfDz.*XXYY)));
H2(1,7) = sum(sum(sum(DfDxSq.*XXZZ)));       H2(1,8) = sum(sum(sum(DfDxDfDy.*XXZZ)));         H2(1,9) = sum(sum(sum(DfDxDfDz.*XXZZ))); 
H2(1,10) = sum(sum(sum(DfDxSq.*(XX-x0))));   H2(1,11) = sum(sum(sum(DfDxDfDy.*(XX-x0))));     H2(1,12) = sum(sum(sum(DfDxDfDz.*(XX-x0))));
H2(2,2) = sum(sum(sum(DfDySq.*XXSq)));       H2(2,3) = sum(sum(sum(DfDyDfDz.*XXSq)));      	H2(2,4) = sum(sum(sum(DfDxDfDy.*XXYY)));
H2(2,5) = sum(sum(sum(DfDySq.*XXYY)));       H2(2,6) = sum(sum(sum(DfDyDfDz.*XXYY)));         H2(2,7) = sum(sum(sum(DfDxDfDy.*XXZZ)));
H2(2,8) = sum(sum(sum(DfDySq.*XXZZ)));       H2(2,9) = sum(sum(sum(DfDyDfDz.*XXZZ)));         H2(2,10) = sum(sum(sum(DfDxDfDy.*(XX-x0))));
H2(2,11) = sum(sum(sum(DfDySq.*(XX-x0))));    H2(2,12) = sum(sum(sum(DfDyDfDz.*(XX-x0))));
H2(3,3) = sum(sum(sum(DfDzSq.*XXSq)));  H2(3,4) = sum(sum(sum(DfDxDfDz.*XXYY)));  H2(3,5) = sum(sum(sum(DfDyDfDz.*XXYY)));
H2(3,6) = sum(sum(sum(DfDzSq.*XXYY)));  H2(3,7) = sum(sum(sum(DfDxDfDz.*XXZZ)));  H2(3,8) = sum(sum(sum(DfDyDfDz.*XXZZ)));
H2(3,9) = sum(sum(sum(DfDzSq.*XXZZ)));  H2(3,10) = sum(sum(sum(DfDxDfDz.*(XX-x0))));  H2(3,11) = sum(sum(sum(DfDyDfDz.*(XX-x0))));
H2(3,12) = sum(sum(sum(DfDzSq.*(XX-x0))));
H2(4,4) = sum(sum(sum(DfDxSq.*YYSq)));  H2(4,5) = sum(sum(sum(DfDxDfDy.*YYSq)));  H2(4,12) = sum(sum(sum(DfDxDfDz.*(YY-y0))));
H2(4,6) = sum(sum(sum(DfDxDfDz.*YYSq)));  H2(4,7) = sum(sum(sum(DfDxSq.*YYZZ)));  H2(4,8) = sum(sum(sum(DfDxDfDy.*YYZZ)));
H2(4,9) = sum(sum(sum(DfDxDfDz.*YYZZ)));  H2(4,10) = sum(sum(sum(DfDxSq.*(YY-y0))));  H2(4,11) = sum(sum(sum(DfDxDfDy.*(YY-y0))));
H2(5,5) = sum(sum(sum(DfDySq.*YYSq)));   H2(5,12) = sum(sum(sum(DfDyDfDz.*(YY-y0))));
H2(5,6) = sum(sum(sum(DfDyDfDz.*YYSq)));  H2(5,7) = sum(sum(sum(DfDxDfDy.*YYZZ)));  H2(5,8) = sum(sum(sum(DfDySq.*YYZZ)));
H2(5,9) = sum(sum(sum(DfDyDfDz.*YYZZ)));  H2(5,10) = sum(sum(sum(DfDxDfDy.*(YY-y0))));  H2(5,11) = sum(sum(sum(DfDySq.*(YY-y0))));
H2(6,6) = sum(sum(sum(DfDzSq.*YYSq)));
H2(6,7) = sum(sum(sum(DfDxDfDz.*YYZZ)));  H2(6,8) = sum(sum(sum(DfDyDfDz.*YYZZ)));  H2(6,9) = sum(sum(sum(DfDzSq.*YYZZ)));
H2(6,10) = sum(sum(sum(DfDxDfDz.*(YY-y0))));  H2(6,11) = sum(sum(sum(DfDyDfDz.*(YY-y0))));  H2(6,12) = sum(sum(sum(DfDzSq.*(YY-y0))));
H2(7,7) = sum(sum(sum(DfDxSq.*ZZSq)));  H2(7,8) = sum(sum(sum(DfDxDfDy.*ZZSq))); H2(7,12) = sum(sum(sum(DfDxDfDz.*(ZZ-z0))));
H2(7,9) = sum(sum(sum(DfDxDfDz.*ZZSq)));  H2(7,10) = sum(sum(sum(DfDxSq.*(ZZ-z0))));  H2(7,11) = sum(sum(sum(DfDxDfDy.*(ZZ-z0))));
H2(8,8) = sum(sum(sum(DfDySq.*ZZSq))); H2(8,12) = sum(sum(sum(DfDyDfDz.*(ZZ-z0))));
H2(8,9) = sum(sum(sum(DfDyDfDz.*ZZSq)));  H2(8,10) = sum(sum(sum(DfDxDfDy.*(ZZ-z0))));  H2(8,11) = sum(sum(sum(DfDySq.*(ZZ-z0))));
H2(9,9) = sum(sum(sum(DfDzSq.*ZZSq)));  H2(9,10) = sum(sum(sum(DfDxDfDz.*(ZZ-z0))));  H2(9,11) = sum(sum(sum(DfDyDfDz.*(ZZ-z0))));
H2(9,12) = sum(sum(sum(DfDzSq.*(ZZ-z0))));
H2(10,10) = sum(sum(sum(DfDxSq)));  H2(10,11) = sum(sum(sum(DfDxDfDy))); H2(10,12) = sum(sum(sum(DfDxDfDz)));
H2(11,11) = sum(sum(sum(DfDySq))); H2(11,12) = sum(sum(sum(DfDyDfDz)));
H2(12,12) = sum(sum(sum(DfDzSq)));

H = H2 + H2' - diag(diag(H2));

% H = zeros(12,12);
% tempCoordx = XX(:); tempCoordy = YY(:); tempCoordz = ZZ(:);
% for tempij = 1:size(tempCoordx,1)
%  
%         H = H + ([Df.DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty,tempCoordz(tempij)-DfDxStartz), ...
%                   Df.DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty,tempCoordz(tempij)-DfDxStartz), ...
%                   Df.DfDz(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty,tempCoordz(tempij)-DfDxStartz)]*...
%             [tempCoordx(tempij)-x0 0 0 tempCoordy(tempij)-y0 0 0 tempCoordz(tempij)-z0 0 0 1 0 0;  
%             0 tempCoordx(tempij)-x0 0 0 tempCoordy(tempij)-y0 0 0 tempCoordz(tempij)-z0 0  0 1 0;
%             0 0 tempCoordx(tempij)-x0 0 0 tempCoordy(tempij)-y0 0 0 tempCoordz(tempij)-z0 0  0 1])'* ...
%             ([Df.DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty,tempCoordz(tempij)-DfDxStartz), ...
%               Df.DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty,tempCoordz(tempij)-DfDxStartz), ...
%               Df.DfDz(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty,tempCoordz(tempij)-DfDxStartz)]*...
%             [tempCoordx(tempij)-x0 0 0 tempCoordy(tempij)-y0 0 0 tempCoordz(tempij)-z0 0 0 1 0 0;  
%             0 tempCoordx(tempij)-x0 0 0 tempCoordy(tempij)-y0 0 0 tempCoordz(tempij)-z0 0  0 1 0;
%             0 0 tempCoordx(tempij)-x0 0 0 tempCoordy(tempij)-y0 0 0 tempCoordz(tempij)-z0 0  0 1]);
% 
% end



% --------------------------
% Initialize while loop
normOfWOld = 2; normOfWNew = 1; stepwithinwhile=0; normOfWNewInit = 1e20;
switch method   % For Gauss-Newton method
    case 'GaussNewton'
        delta = 0;
    case 'LevenbergMarquardt'
        delta = 0.001; % For Levenberg-Marquardt method 
        KappaOld = 1e10; KappaNew = 1e10; KappaStore = zeros(10,1); PStore = zeros(10,12);
    otherwise
        delta = 0;
end
 
try 
    
    while(  stepwithinwhile <= MaxIterNum && normOfWNew > tol  )
           
    stepwithinwhile = stepwithinwhile+1;

    % Find region for g
      
    %[tempCoordy, tempCoordx] = meshgrid(y(1):y(3),x(1):x(3));
    tempCoordxMat = XX - x0*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
    tempCoordyMat = YY - y0*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
    tempCoordzMat = ZZ - z0*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
    
    u22 = (1+P(1))*tempCoordxMat + P(4)*tempCoordyMat + P(7)*tempCoordzMat + (x0+P(10))*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
    v22 = P(2)*tempCoordxMat + (1+P(5))*tempCoordyMat + P(8)*tempCoordzMat + (y0+P(11))*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
    w22 = P(3)*tempCoordxMat + P(6)*tempCoordyMat + (1+P(9))*tempCoordzMat + (z0+P(12))*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
    
    row1 = find(u22<1); row2 = find(u22>imgSize(1)); 
    row3 = find(v22<1); row4 = find(v22>imgSize(2));
    row5 = find(w22<1); row6 = find(w22>imgSize(3));
    
    if ~isempty([row1; row2; row3; row4; row5; row6])
        normOfWNew = nan;
         % warning('Out of image boundary!')
        break;
    else
        
        switch interpmethod
            case 'linear'
                %tempg = mirt3D_mexinterp(imggNormalizedbc, v22, u22, w22);
                tempg = mirt3D_mexinterp(imggNormalizedbc, v22-ImggStarty, u22-ImggStartx, w22-ImggStartz);
            case 'cubicspline'
                % tempg = interp3(imggNormalizedbc, v22-ImggStarty, u22-ImggStartx, w22-ImggStartz, 'spline'); % 400 times slower than tempg
                tempg = interp3(imggNormalizedbc, v22-ImggStarty, u22-ImggStartx, w22-ImggStartz, 'spline');
            case 'cubic'
                % tempg = ba_interp3(imggNormalizedbc, v22, u22, w22, 'cubic');
                tempg = ba_interp3(imggNormalizedbc, v22-ImggStarty, u22-ImggStartx, w22-ImggStartz, 'cubic');
            otherwise
                % tempg = interp3(imggNormalizedbc, v22, u22, w22, 'linear');
                tempg = interp3(imggNormalizedbc, v22-ImggStarty, u22-ImggStartx, w22-ImggStartz, 'linear');
        end
          
        % ====== Old version codes ======
        % tempg = zeros(size(tempf,1)*size(tempf,2),1);
        % [tempCoordy, tempCoordx] = meshgrid(1:winsize+1,1:winsize+1);
        % tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
        % 
        % for tempij = 1:size(tempCoordx,1)
        %     tempg(tempij)= ...
        %         fungInterpolation_g(u22(tempCoordx(tempij),tempCoordy(tempij)), v22(tempCoordx(tempij),tempCoordy(tempij)), ...
        %         g(floor(u22(tempCoordx(tempij),tempCoordy(tempij)))-1:floor(u22(tempCoordx(tempij),tempCoordy(tempij)))+2, ...
        %         floor(v22(tempCoordx(tempij),tempCoordy(tempij)))-1:floor(v22(tempCoordx(tempij),tempCoordy(tempij)))+2));
        % end
        % 
        % tempg = reshape(tempg, winsize+1, winsize+1);
        % ===============================

        % A = [1+P(1) P(2) 0; P(3) 1+P(4) 0; P(5) P(6) 1];
        % tform = affine2d((A));
        % 
        % tempg2 = g((x(1)-winsize/2):(x(3)+winsize/2), (y(1)-winsize/2):(y(3)+winsize/2));
        % tempg3 = imwarp(tempg2,tform,'cubic');
        % 
        % figure; imshow(tempf,[]);
        % figure; imshow(tempg2,[]);
        % figure; imshow(tempg3,[]);
        % 
        % [M,N] = size(tempg3) 
        % tempg = tempg3(ceil((M+1)/2)-winsize/2:ceil((M+1)/2)+winsize/2, ceil((N+1)/2)-winsize/2:ceil((N+1)/2)+winsize/2);
        % figure; imshow(tempg,[]);
        
        meang = mean(tempg(:));
        bottomg = sqrt((length(tempg(:))-1)*var(tempg(:)));

        % ============ For Levenberg-Marquardt method ============    
        switch method
            case 'LevenbergMarquardt'
            % Compute functinoal error
            KappaOld = KappaNew;
            Kappatemp = (tempf-meanf)/bottomf - (tempg-meang)/bottomg;
            Kappatemp = Kappatemp.*Kappatemp;
            KappaNew = sum(Kappatemp(:));  

               if KappaNew < 1.02*KappaOld
                   delta = delta/10; 
               else
                   delta = delta*10;
                   % Perform P inverse 
                   DP = -DP;
                   
                   tempP1 = (DP(3)*(DP(7) + DP(5)*DP(7) - DP(4)*DP(8)) + DP(2)*(DP(4) - DP(6)*DP(7) + DP(4)*DP(9)) - ...
                       DP(1)*(1 + DP(5) - DP(6)*DP(8) + DP(9) + DP(5)*DP(9)))/detDP;
                   tempP4 = (DP(6)*DP(7) - DP(4)*(1 + DP(9)))/detDP;
                   tempP7 = (-(1 + DP(5))*DP(7) + DP(4)*DP(8))/detDP;
                   tempP10 = (DP(12)*(DP(7) + DP(5)*DP(7) - DP(4)*DP(8)) + DP(11)*(DP(4) - DP(6)*DP(7) + DP(4)*DP(9)) - ...
                       DP(10)*(1 + DP(5) - DP(6)*DP(8) + DP(9) + DP(5)*DP(9)))/detDP;
                   tempP2 = (DP(3)*DP(8) - DP(2)*(1 + DP(9)))/detDP;
                   tempP5 = ((-DP(3)*DP(4) + DP(6) + DP(1)*DP(6))*DP(8) - DP(5)*(1 + DP(1) - DP(3)*DP(7) + DP(9) + DP(1)*DP(9)) + ...
                       DP(2)*(DP(4) - DP(6)*DP(7) + DP(4)*DP(9)))/detDP;
                   tempP8 = (DP(2)*DP(7) - (1 + DP(1))*DP(8))/detDP;
                   tempP11 = (DP(12)*(-DP(2)*DP(7) + DP(8) + DP(1)*DP(8)) - DP(11)*(1 + DP(1) - DP(3)*DP(7) + DP(9) + DP(1)*DP(9)) + ...
                       DP(10)*(DP(2) - DP(3)*DP(8) + DP(2)*DP(9)))/detDP;
                   tempP3 = (-DP(3)*(1 + DP(5)) + DP(2)*DP(6))/detDP;
                   tempP6 = (DP(3)*DP(4) - (1 + DP(1))*DP(6))/detDP;
                   tempP9 = (DP(3)*(DP(7) + DP(5)*DP(7) - DP(4)*DP(8)) + ...
                       DP(2)*(-DP(6)*DP(7) + DP(4)*DP(9)) + (1 + DP(1))*(DP(6)*DP(8) - (1 + DP(5))*DP(9)))/detDP;
                   tempP12 = (-DP(12)*(1 + DP(1) - ...
                       DP(2)*DP(4) + DP(5) + DP(1)*DP(5)) + DP(11)*(-DP(3)*DP(4) + DP(6) + DP(1)*DP(6)) + ...
                       DP(10)*(DP(3) + DP(3)*DP(5) - DP(2)*DP(6)))/detDP;
                   
                   tempMatrix = [1+P(1) P(4) P(7) P(10); P(2) 1+P(5) P(8) P(11); P(3) P(6) 1+P(9) P(12); 0 0 0 1] * ...
                       [1+tempP1 tempP4 tempP7 tempP10; tempP2 1+tempP5 tempP8 tempP11; tempP3 tempP6 1+tempP9 tempP12; 0 0 0 1];
                   
                   P1 = tempMatrix(1,1)-1; P4 = tempMatrix(1,2);   P7 = tempMatrix(1,3);   P10 = tempMatrix(1,4);
                   P2 = tempMatrix(2,1);   P5 = tempMatrix(2,2)-1; P8 = tempMatrix(2,3);   P11 = tempMatrix(2,4);
                   P3 = tempMatrix(3,1);   P6 = tempMatrix(3,2);   P9 = tempMatrix(3,3)-1; P12 = tempMatrix(3,4);
                   P = [P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12]';
                 
               end

               % Find region for g
               u22 = (1+P(1))*tempCoordxMat + P(4)*tempCoordyMat + P(7)*tempCoordzMat + (x0+P(10))*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
               v22 = P(2)*tempCoordxMat + (1+P(5))*tempCoordyMat + P(8)*tempCoordzMat + (y0+P(11))*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
               w22 = P(3)*tempCoordxMat + P(6)*tempCoordyMat + (1+P(9))*tempCoordzMat + (z0+P(12))*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
               
               switch interpmethod
                   case 'linear'
                       %tempg = mirt3D_mexinterp(imggNormalizedbc, v22, u22, w22);
                       tempg = mirt3D_mexinterp(imggNormalizedbc, v22-ImggStarty, u22-ImggStartx, w22-ImggStartz);
                   case 'cubicspline'
                       % tempg = interp3(imggNormalizedbc, v22-ImggStarty, u22-ImggStartx, w22-ImggStartz, 'spline'); % 400 times slower than tempg
                       tempg = interp3(imggNormalizedbc, v22-ImggStarty, u22-ImggStartx, w22-ImggStartz, 'spline');
                   case 'cubic'
                       % tempg = ba_interp3(imggNormalizedbc, v22, u22, w22, 'cubic');
                       tempg = ba_interp3(imggNormalizedbc, v22-ImggStarty, u22-ImggStartx, w22-ImggStartz, 'cubic');
                   otherwise
                       % tempg = interp3(imggNormalizedbc, v22, u22, w22, 'linear');
                       tempg = interp3(imggNormalizedbc, v22-ImggStarty, u22-ImggStartx, w22-ImggStartz, 'linear');
               end
        
               % ====== Old version codes ======
               % tempg = zeros(size(tempf,1)*size(tempf,2),1);
               %
               % [tempCoordy, tempCoordx] = meshgrid(1:winsize+1,1:winsize+1);
               % tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
               %
               % parfor tempij = 1:size(tempCoordx,1)
               %     tempg(tempij)= ...
               %         fungInterpolation_g(u22(tempCoordx(tempij),tempCoordy(tempij)), v22(tempCoordx(tempij),tempCoordy(tempij)), ...
               %         g(floor(u22(tempCoordx(tempij),tempCoordy(tempij)))-1:floor(u22(tempCoordx(tempij),tempCoordy(tempij)))+2, ...
               %         floor(v22(tempCoordx(tempij),tempCoordy(tempij)))-1:floor(v22(tempCoordx(tempij),tempCoordy(tempij)))+2));
               % end
               %
               % tempg = reshape(tempg, winsize+1, winsize+1);
               % ==================================
               meang = mean(tempg(:));
               bottomg = sqrt((length(tempg(:))-1)*var(tempg(:)));
               
            otherwise
         end
        
        % % ============ End of Levenberg-Marquardt method ============ 

        % Assemble b vector
%         b = zeros(6,1);
%         
%         %[tempCoordy, tempCoordx] = meshgrid(y(1):y(3),x(1):x(3));
%         %tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
%          
%         for tempij = 1:size(tempCoordx,1)
%             b = b + bottomf*([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
%                     [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
%                     ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
%                     (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
%         end    
         
         
        b2 = zeros(12,1); 
        % tempfMinustempgOrig = tempf-tempg;
        tempfMinustempg = (tempf-meanf*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1))/bottomf - ...
                        (tempg-meang*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1))/bottomg;
        b2(1) = sum(sum(sum( (XX-x0).*DfDx.*tempfMinustempg )));
        b2(2) = sum(sum(sum( (XX-x0).*DfDy.*tempfMinustempg )));
        b2(3) = sum(sum(sum( (XX-x0).*DfDz.*tempfMinustempg )));
        b2(4) = sum(sum(sum( (YY-y0).*DfDx.*tempfMinustempg )));
        b2(5) = sum(sum(sum( (YY-y0).*DfDy.*tempfMinustempg )));
        b2(6) = sum(sum(sum( (YY-y0).*DfDz.*tempfMinustempg )));
        b2(7) = sum(sum(sum( (ZZ-z0).*DfDx.*tempfMinustempg )));
        b2(8) = sum(sum(sum( (ZZ-z0).*DfDy.*tempfMinustempg )));
        b2(9) = sum(sum(sum( (ZZ-z0).*DfDz.*tempfMinustempg )));
        b2(10) = sum(sum(sum( DfDx.*tempfMinustempg )));
        b2(11) = sum(sum(sum( DfDy.*tempfMinustempg )));
        b2(12) = sum(sum(sum( DfDz.*tempfMinustempg )));
         
        b = bottomf * b2;
        
        normOfWOld = normOfWNew;
        normOfWNew = norm(b(:));
        
        if stepwithinwhile == 1
            normOfWNewInit = normOfWNew;
        end
        
        normOfWNew = normOfWNew/normOfWNewInit;
          
        if (normOfWNew < tol) || (normOfWNew*normOfWNewInit < 1e-5)
            break
        else
             
            % DeltaP = [0 0 0 0 0 0];
            % tempH = (H + delta*diag(diag(H)));
            % tempb = b;
            % DeltaP(5:6) = -tempH(5:6,5:6)\tempb(5:6);
            DP = -(H + delta*diag(diag(H))) \ b;
            % temp = ((1+DeltaP(1))*(1+DeltaP(4)) - DeltaP(2)*DeltaP(3));
            detDP = 1+DP(5)+(-1).*DP(3).*DP(7)+(-1).*DP(3).*DP(5).*DP(7)+DP(3).*DP(4).*DP(8)+...
                (-1).*DP(6).*DP(8)+DP(9)+DP(5).*DP(9)+(-1).*DP(2).*(DP(4)+(-1).*DP(6).*DP(7)+DP(4).* ...
                DP(9))+DP(1).*(1+DP(5)+(-1).*DP(6).*DP(8)+DP(9)+DP(5).*DP(9));

            if (detDP ~= 0)
                % tempP1 =  (-DeltaP(1)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/temp;
                % tempP2 =  -DeltaP(2)/temp;
                % tempP3 =  -DeltaP(3)/temp;
                % tempP4 =  (-DeltaP(4)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/temp;
                % tempP5 =  (-DeltaP(5)-DeltaP(4)*DeltaP(5)+DeltaP(3)*DeltaP(6))/temp;
                % tempP6 =  (-DeltaP(6)-DeltaP(1)*DeltaP(6)+DeltaP(2)*DeltaP(5))/temp;
                tempP1 = (DP(3)*(DP(7) + DP(5)*DP(7) - DP(4)*DP(8)) + DP(2)*(DP(4) - DP(6)*DP(7) + DP(4)*DP(9)) - ...
                    DP(1)*(1 + DP(5) - DP(6)*DP(8) + DP(9) + DP(5)*DP(9)))/detDP;
                tempP4 = (DP(6)*DP(7) - DP(4)*(1 + DP(9)))/detDP;
                tempP7 = (-(1 + DP(5))*DP(7) + DP(4)*DP(8))/detDP;
                tempP10 = (DP(12)*(DP(7) + DP(5)*DP(7) - DP(4)*DP(8)) + DP(11)*(DP(4) - DP(6)*DP(7) + DP(4)*DP(9)) - ...
                    DP(10)*(1 + DP(5) - DP(6)*DP(8) + DP(9) + DP(5)*DP(9)))/detDP;
                tempP2 = (DP(3)*DP(8) - DP(2)*(1 + DP(9)))/detDP;
                tempP5 = ((-DP(3)*DP(4) + DP(6) + DP(1)*DP(6))*DP(8) - DP(5)*(1 + DP(1) - DP(3)*DP(7) + DP(9) + DP(1)*DP(9)) + ...
                    DP(2)*(DP(4) - DP(6)*DP(7) + DP(4)*DP(9)))/detDP;
                tempP8 = (DP(2)*DP(7) - (1 + DP(1))*DP(8))/detDP;
                tempP11 = (DP(12)*(-DP(2)*DP(7) + DP(8) + DP(1)*DP(8)) - DP(11)*(1 + DP(1) - DP(3)*DP(7) + DP(9) + DP(1)*DP(9)) + ...
                    DP(10)*(DP(2) - DP(3)*DP(8) + DP(2)*DP(9)))/detDP;
                tempP3 = (-DP(3)*(1 + DP(5)) + DP(2)*DP(6))/detDP;
                tempP6 = (DP(3)*DP(4) - (1 + DP(1))*DP(6))/detDP;
                tempP9 = (DP(3)*(DP(7) + DP(5)*DP(7) - DP(4)*DP(8)) + ...
                    DP(2)*(-DP(6)*DP(7) + DP(4)*DP(9)) + (1 + DP(1))*(DP(6)*DP(8) - (1 + DP(5))*DP(9)))/detDP;
                tempP12 = (-DP(12)*(1 + DP(1) - ...
                    DP(2)*DP(4) + DP(5) + DP(1)*DP(5)) + DP(11)*(-DP(3)*DP(4) + DP(6) + DP(1)*DP(6)) + ...
                    DP(10)*(DP(3) + DP(3)*DP(5) - DP(2)*DP(6)))/detDP;
   
                % tempMatrix = [1+P(1) P(3) P(5); P(2) 1+P(4) P(6); 0 0 1]*...
                %    [1+tempP1 tempP3 tempP5; tempP2 1+tempP4 tempP6; 0 0 1];
                tempMatrix = [1+P(1) P(4) P(7) P(10); P(2) 1+P(5) P(8) P(11); P(3) P(6) 1+P(9) P(12); 0 0 0 1] * ...
                   [1+tempP1 tempP4 tempP7 tempP10; tempP2 1+tempP5 tempP8 tempP11; tempP3 tempP6 1+tempP9 tempP12; 0 0 0 1];
                
%                 tempMatrix = [1+P(1) P(4) P(7) P(10); P(2) 1+P(5) P(8) P(11); P(3) P(6) 1+P(9) P(12); 0 0 0 1]  ...
%                     /([1+DP(1), DP(4), DP(7) DP(10); DP(2), 1+DP(5), DP(8) DP(11); DP(3), DP(6), 1+DP(9) DP(12); 0,0,0,1]);
            
%             [1+P(1) P(4) P(7) P(10); P(2) 1+P(5) P(8) P(11); P(3) P(6) 1+P(9) P(12); 0 0 0 1]  ...
%                     *inv([1+DP(1), DP(4), DP(7) DP(10); DP(2), 1+DP(5), DP(8) DP(11); DP(3), DP(6), 1+DP(9) DP(12); 0,0,0,1]);

                % P1 = tempMatrix(1,1)-1;
                % P2 = tempMatrix(2,1);
                % P3 = tempMatrix(1,2);
                % P4 = tempMatrix(2,2)-1;
                % P5 = tempMatrix(1,3);
                % P6 = tempMatrix(2,3);
                % P = [P1 P2 P3 P4 P5 P6]';
                P1 = tempMatrix(1,1)-1; P4 = tempMatrix(1,2);   P7 = tempMatrix(1,3);   P10 = tempMatrix(1,4);
                P2 = tempMatrix(2,1);   P5 = tempMatrix(2,2)-1; P8 = tempMatrix(2,3);   P11 = tempMatrix(2,4);
                P3 = tempMatrix(3,1);   P6 = tempMatrix(3,2);   P9 = tempMatrix(3,3)-1; P12 = tempMatrix(3,4);
                P = [P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12]';
                 
            else
                disp(['Det(DeltaP)==0!'])
                break
            end

        end
    end  
end % end of while


catch
    
    normOfWNew = nan; P = zeros(12,1);
    
    
end

U = [0,0,0]'; F = [0,0,0,0,0,0,0,0,0]';
U(1) = P(10); U(2) = P(11); U(3) = P(12);
F(1) = P(1); F(2) = P(2); F(3) = P(3); F(4) = P(4); F(5) = P(5); F(6) = P(6); F(7) = P(7); F(8) = P(8); F(9) = P(9); 


HGlobal = [H(1:12) H(12*1+2:12*2) H(12*2+3:12*3) H(12*3+4:12*4) H(12*4+5:12*5) H(12*5+6:12*6) H(12*6+7:12*7) ...
    H(12*7+8:12*8) H(12*8+9:12*9) H(12*9+10:12*10) H(12*10+11:12*11) H(12*12)];

 
if (normOfWNew < tol) % || (normOfWNew*normOfWNewInit < 1e-5)
    % elementsLocalMethodConvergeOrNot = 1;
else
    stepwithinwhile = MaxIterNum+1 ;
end
 
if (isnan(normOfWNew)==1)
    stepwithinwhile = MaxIterNum+2 ;
end

 
end