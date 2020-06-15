%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function AL-DIC Subproblem 1 3D DVC
% Object: to find deformation field using local methods
% Author: Jin Yang
% Last date modified: 2018.03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U, stepwithinwhile, HGlobal] = funICGN_Subpb13(Coords,DfEle,ImgEle,Img2,winsize,H,beta,mu,udual,vdual,UOld,FOld,tol,method,interpmethod,MaxIterNum)

warning('off','all');

DIM = 3; x0 = Coords(1); y0 = Coords(2); z0 = Coords(3);
imgfNormalizedbc = ImgEle.Imgf(4:end-3, 4:end-3, 4:end-3); 
imggNormalizedbc = Img2;

% DfDxStartx = Df.DfAxis(1); DfDxStarty = Df.DfAxis(3); DfDxStartz = Df.DfAxis(5);
% imgSize = Df.imgSize;
ImggStartx = 0; ImggStarty = 0; ImggStartz = 0; 
imgSize = DfEle.imgSize;

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
P0 = full([reshape(FOld,1,DIM^2), reshape(UOld,1,DIM)]'); P = P0;  
 
% ---------------------------
% Find region for f
[XX,YY,ZZ] = ndgrid([x(1):1:x(7)],[y(1):1:y(7)],[z(1):1:z(7)]);

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

if norm(H,2)<abs(eps)
    
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
    
end

HGlobal = [H(1:12) H(12*1+2:12*2) H(12*2+3:12*3) H(12*3+4:12*4) H(12*4+5:12*5) H(12*5+6:12*6) H(12*6+7:12*7) ...
    H(12*7+8:12*8) H(12*8+9:12*9) H(12*9+10:12*10) H(12*10+11:12*11) H(12*12)];


H2 = H(10:12,10:12)*2/(bottomf^2) + [mu 0 0; 0 mu 0; 0 0 mu];

%%%%%%%%%%%%%%%% Tried below, not succeed %%%%%%%%%%%%%%%%%%%
% if CrackOrNot == 0
%     H2 = H(5:6,5:6)*2/(bottomf^2) + [mu 0; 0 mu];
% else
%     H2 = H*2/(bottomf^2) + [beta 0 0 0 0 0; 0 beta 0 0 0 0; 0 0 beta 0 0 0; 0 0 0 beta 0 0 ; 0 0 0 0 mu 0; 0 0 0 0 0 mu];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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
 
while(  stepwithinwhile <= MaxIterNum && normOfWNew > tol  )
           
    stepwithinwhile = stepwithinwhile+1;

    % Find region for g
    tempCoordxMat = XX - x0*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
    tempCoordyMat = YY - y0*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
    tempCoordzMat = ZZ - z0*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
    
    
    u22 = (1+P(1))*tempCoordxMat + P(4)*tempCoordyMat + P(7)*tempCoordzMat + (x0+P(10))*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
    v22 = P(2)*tempCoordxMat + (1+P(5))*tempCoordyMat + P(8)*tempCoordzMat + (y0+P(11))*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
    w22 = P(3)*tempCoordxMat + P(6)*tempCoordyMat + (1+P(9))*tempCoordzMat + (z0+P(12))*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1);
   
    row1 = find(u22<1); row2 = find(u22>imgSize(1) ); 
    row3 = find(v22<1); row4 = find(v22>imgSize(2) );
    row5 = find(w22<1); row6 = find(w22>imgSize(3) );
   
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
        %     tempg(tempij,1)= ...
        %         fungInterpolation_g(  u22(tempCoordx(tempij),tempCoordy(tempij)) , v22(tempCoordx(tempij),tempCoordy(tempij)), ...
        %         g((floor(u22(tempCoordx(tempij),tempCoordy(tempij)))-1):(floor(u22(tempCoordx(tempij),tempCoordy(tempij)))+2), ...
        %         (floor(v22(tempCoordx(tempij),tempCoordy(tempij)))-1):(floor(v22(tempCoordx(tempij),tempCoordy(tempij)))+2))  );
        % end
        % 
        % tempg = reshape(tempg, winsize+1, winsize+1);
        
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
               u22 = (1+P(1))*tempCoordxMat + P(4)*tempCoordyMat + P(7)*tempCoordzMat + (x0+P(10))*ones(winsize+1,winsize+1,winsize+1);
               v22 = P(2)*tempCoordxMat + (1+P(5))*tempCoordyMat + P(8)*tempCoordzMat + (y0+P(11))*ones(winsize+1,winsize+1,winsize+1);
               w22 = P(3)*tempCoordxMat + P(6)*tempCoordyMat + (1+P(9))*tempCoordzMat + (z0+P(12))*ones(winsize+1,winsize+1,winsize+1);
               
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
        
        % ============ End of Levenberg-Marquardt method ============ 

        % Assemble b vector
%         b = zeros(6,1);
%         [tempCoordy, tempCoordx] = meshgrid(y(1):y(3),x(1):x(3));
%         tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
%  
%         % Compute for all the pixels in the master part
%         for tempij = 1:size(tempCoordx,1)
%             if PassCrackOrNot == 1 % Crack pass through this local subset
%                 % if x0 < CrackTip(1) % Left part is master, and right part will be discarded
%                 if CrackPathCen(1)*x0 + CrackPathCen(2)*y0 + 1 > 0  % Bottom part is master, and top part will be discarded
%                     % if tempCoordx(tempij) > CrackTip(1) && tempCoordy(tempij) < CrackTip(2)
%                     %     tempCoordx(tempij) = tempCoordx(tempij)-winsize;
%                     % end
%                     % if tempCoordx(tempij)+1-x(1)>1 && tempCoordx(tempij)>x0-0.8*winsize
%                     %     b = b+ bottomf*([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
%                     %         [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
%                     %         ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
%                     %         (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
%                     % end
%                     if tempCoordx(tempij)*CrackPathCen(1) + tempCoordy(tempij)*CrackPathCen(2) + 1 < 0
%                         tempCoordy(tempij) = tempCoordy(tempij)-winsize;
%                     end
%                     if (tempCoordy(tempij)+1-y(1) > 1) && (tempCoordy(tempij)+1-y(1) < size(tempf,2)) ...
%                         && (tempCoordx(tempij)*CrackPath2(1) + tempCoordy(tempij)*CrackPath2(2) + 1 > 0)
%                            b = b+ bottomf*([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
%                             [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
%                             ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
%                             (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
%                     end    
%                 % elseif (x0 >= CrackTip(1)) % Right part is master, and left part will be discarded
%                 elseif CrackPathCen(1)*x0 + CrackPathCen(2)*y0 + 1 <= 0  % Top part is master, and bottom part will be discarded
%                     % if tempCoordx(tempij) < CrackTip(1) && tempCoordy(tempij) < CrackTip(2)
%                     %     tempCoordx(tempij) = tempCoordx(tempij)+winsize;
%                     % end
%                     % if tempCoordx(tempij)<size(DfDx,1) && tempCoordx(tempij)+1-x(1)<size(tempf,1) && tempCoordx(tempij)<x0+0.5*winsize
%                     %     b = b+ bottomf*([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
%                     %         [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
%                     %         ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
%                     %         (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
%                     % end
%                     if tempCoordx(tempij)*CrackPathCen(1) + tempCoordy(tempij)*CrackPathCen(2) + 1 > 0
%                         tempCoordy(tempij) = tempCoordy(tempij)+winsize;
%                     end
%                     if (tempCoordy(tempij)-DfDxStarty < size(DfDy,2)) && (tempCoordy(tempij)+1-y(1) < size(tempf,2)) ...
%                         && (tempCoordx(tempij)*CrackPath1(1) + tempCoordy(tempij)*CrackPath1(2) + 1 < 0)
%                             b = b + bottomf*([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
%                             [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
%                             ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
%                             (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
%                     end
%                 end
%             else
%                 b = b + bottomf*([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
%                     [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
%                     ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
%                     (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
%             end
%         end
        
        b2 = zeros(12,1); 
        tempfMinustempg = (tempf-meanf*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1))/bottomf - ...
                        (tempg-meang*ones(winsize(1)+1,winsize(2)+1,winsize(3)+1))/bottomg;
        % b2(1) = sum(sum(sum( (XX-x0).*DfDx.*tempfMinustempg )));
        % b2(2) = sum(sum(sum( (XX-x0).*DfDy.*tempfMinustempg )));
        % b2(3) = sum(sum(sum( (XX-x0).*DfDz.*tempfMinustempg )));
        % b2(4) = sum(sum(sum( (YY-y0).*DfDx.*tempfMinustempg )));
        % b2(5) = sum(sum(sum( (YY-y0).*DfDy.*tempfMinustempg )));
        % b2(6) = sum(sum(sum( (YY-y0).*DfDz.*tempfMinustempg )));
        % b2(7) = sum(sum(sum( (ZZ-z0).*DfDx.*tempfMinustempg )));
        % b2(8) = sum(sum(sum( (ZZ-z0).*DfDy.*tempfMinustempg )));
        % b2(9) = sum(sum(sum( (ZZ-z0).*DfDz.*tempfMinustempg )));
        b2(10) = sum(sum(sum( DfDx.*tempfMinustempg )));
        b2(11) = sum(sum(sum( DfDy.*tempfMinustempg )));
        b2(12) = sum(sum(sum( DfDz.*tempfMinustempg )));
        
        b = bottomf * b2;
        
        tempb = b(10:12)*2/(bottomf^2) + mu*[(P(10)-UOld(1)-vdual(1)); (P(11)-UOld(2)-vdual(2)); (P(12)-UOld(3)-vdual(3))];
        %%%%%%%%%%%%%%%% Tried below, not succeed %%%%%%%%%%%%%%%%%%%
        % if CrackOrNot == 0
        %     tempb = b(5:6)*2/(bottomf^2)  + [mu*(P(5)-UOld(1)-vdual(1)); mu*(P(6)-UOld(2)-vdual(2))];
        % else
        %     tempb = b*2/(bottomf^2) +[beta*(P(1)-FOld(1)-udual(1)); beta*(P(2)-FOld(2)-udual(2));
        %                               beta*(P(3)-FOld(3)-udual(3)); beta*(P(4)-FOld(4)-udual(4));
        %                               mu*(P(5)-UOld(1)-vdual(1)); mu *(P(6)-UOld(2)-vdual(2))];
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
         
        %%%%%%%%%%%%%%%% Tried below, not succeed %%%%%%%%%%%%%%%%%%%
        % if CrackOrNot == 0
        %     DeltaP(5:6) = -tempH\tempb;
        % else
        %     DeltaP = -(H2 + delta*diag(diag(H2))) \ tempb;
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        normOfWOld = normOfWNew;
        normOfWNew = norm(tempb(:));
         
        if stepwithinwhile == 1
            normOfWNewInit = normOfWNew;
        end
        
        normOfWNew = normOfWNew/normOfWNewInit;
        
        if (normOfWNew < tol) || (normOfWNew*normOfWNewInit < mu*1e-4)
            break
        else
             
            DP = zeros(1,12);
            H2 = H(10:12,10:12)*2/(bottomf^2) + mu*eye(3);
            tempH = (H2 + delta*diag(diag(H2))*eye(3));

            DP(10:12) = -tempH\tempb;
         
            detDP = 1+DP(5)+(-1).*DP(3).*DP(7)+(-1).*DP(3).*DP(5).*DP(7)+DP(3).*DP(4).*DP(8)+ ...
                (-1).*DP(6).*DP(8)+DP(9)+DP(5).*DP(9)+(-1).*DP(2).*(DP(4)+(-1).*DP(6).*DP(7)+DP(4).* ...
                DP(9))+DP(1).*(1+DP(5)+(-1).*DP(6).*DP(8)+DP(9)+DP(5).*DP(9));
            
            if (detDP ~= 0)
                
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
                
            else
                disp( 'Det(DeltaP)==0!' );
                break
            end

        end
    end  
end % end of while
 
U(1) = P(10); U(2) = P(11); U(3) = P(12);
% F(1) = P(1); F(2) = P(2); F(3) = P(3); F(4) = P(4); F(5) = P(5); F(6) = P(6); F(7) = P(7); F(8) = P(8); F(9) = P(9); 

if (normOfWNew < tol) || (normOfWNew*normOfWNewInit < mu*1e-4)
else
    stepwithinwhile = MaxIterNum+1 ;
end
 
if (isnan(normOfWNew)==1)
    stepwithinwhile = MaxIterNum+2 ;
end

% if (elementsLocalMethodConvergeOrNot == 0)
% 	[row] = find(KappaStore == min(KappaStore(1:min(10,stepwithinwhile-1))));
% 	U(1) = PStore(row,5); U(2) = PStore(row,6);
% 	F(1) = PStore(row,1); F(2) = PStore(row,2); F(3) = PStore(row,3); F(4) = PStore(row,4);
% end

end


 