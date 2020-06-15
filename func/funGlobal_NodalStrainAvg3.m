%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function AL-DIC Subproblem 2  3d case                    %
% Object: compute strain field from Subpb2 disp results    %
% Author: Jin Yang                                         %
% Last date modified: 2019.03                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,StrainGaussPt,CoordsGaussPt] = funGlobal_NodalStrainAvg3(coordinatesFEM,elementsFEM,U,GaussPtOrder)

DIM = 3; NodesPerEle = 8;

% ------- Info of CrackTip path-line --------
% k_tip = -0.5*(CrackPath1(1)/CrackPath1(2)+CrackPath2(1)/CrackPath2(2));
% CrackPathCen = [k_tip/(CrackTip(2)-k_tip*CrackTip(1)), -1/(CrackTip(2)-k_tip*CrackTip(1))];
% CrackPathCen = [-1/CrackTip(1),0];

% CrackTipOrNot = 1-CrackTipOrNot; % Keep consistent with old version codes.
% if CrackOrNot == 0
% FEMSize = size(coordinatesFEM,1);
% elseif CrackOrNot>0 && CrackTipOrNot == 1
%     FEMSize = size(coordinates,1) + size(EnrHAndTipEleIndex,1) + 4*size(EnrTipEleIndex,1);
% else
%     FEMSize = size(coordinates,1) + size(EnrHAndTipEleIndex,1);
% end

% -------------------------------------------
U = [U;zeros(NodesPerEle*DIM,1)];
FStrainAvgTimes = zeros(9*size(coordinatesFEM,1),1); FStrain = zeros(9*size(coordinatesFEM,1),1);
% No of Gauss Pt: 2^DIM
StrainGaussPt = zeros(2^DIM*size(elementsFEM,1),9); CoordsGaussPt = zeros(2^DIM*size(elementsFEM,1),3);

% ====== Gaussian quadrature parameter ======
switch GaussPtOrder
    case 0 % ------ All pixels ------
    case 2 % ------ 2 Gauss points ------
        gqpt1 = -1/sqrt(3); gqpt2 = 1/sqrt(3); gqpt = [gqpt1,gqpt2];
        gqwt1 = 1; gqwt2 = 1; gqwt = [gqwt1,gqwt2];
    case 3 % ------ 3 Gauss points ------
        gqpt1 = 0; gqpt2 = 0.57735; gqpt3 = -0.57735; gqpt = [gqpt1,gqpt2,gqpt3];
        gqwt1 = 1; gqwt2 = 1; gqwt3 = 1; gqwt = [gqwt1,gqwt2,gqwt3];
    case 5 % ------ 5 Gauss points ------
        gqpt1 = 0; gqpt2 = 0.538469; gqpt3 = -0.538469; gqpt4 = 0.90618; gqpt5 = -0.90618;
        gqwt1 = 0.568889; gqwt2 = 0.478629; gqwt3 = 0.478629; gqwt4 = 0.236927; gqwt5 = 0.236927;
        gqpt = [gqpt1,gqpt2,gqpt3,gqpt4,gqpt5]; gqwt = [gqwt1,gqwt2,gqwt3,gqwt4,gqwt5];
end
% ===========================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:size(elementsFEM,1) % j is the element index
    
    % EleCrackTipOrNot = 0; EleCrackHOrNot = 0;
    
    % ------ Find corner pts ------
    pt1x = coordinatesFEM(elementsFEM(j,1),1); pt1y = coordinatesFEM(elementsFEM(j,1),2); pt1z = coordinatesFEM(elementsFEM(j,1),3);
    pt2x = coordinatesFEM(elementsFEM(j,2),1); pt2y = coordinatesFEM(elementsFEM(j,2),2); pt2z = coordinatesFEM(elementsFEM(j,2),3);
    pt3x = coordinatesFEM(elementsFEM(j,3),1); pt3y = coordinatesFEM(elementsFEM(j,3),2); pt3z = coordinatesFEM(elementsFEM(j,3),3);
    pt4x = coordinatesFEM(elementsFEM(j,4),1); pt4y = coordinatesFEM(elementsFEM(j,4),2); pt4z = coordinatesFEM(elementsFEM(j,4),3);
    pt5x = coordinatesFEM(elementsFEM(j,5),1); pt5y = coordinatesFEM(elementsFEM(j,5),2); pt5z = coordinatesFEM(elementsFEM(j,5),3);
    pt6x = coordinatesFEM(elementsFEM(j,6),1); pt6y = coordinatesFEM(elementsFEM(j,6),2); pt6z = coordinatesFEM(elementsFEM(j,6),3);
    pt7x = coordinatesFEM(elementsFEM(j,7),1); pt7y = coordinatesFEM(elementsFEM(j,7),2); pt7z = coordinatesFEM(elementsFEM(j,7),3);
    pt8x = coordinatesFEM(elementsFEM(j,8),1); pt8y = coordinatesFEM(elementsFEM(j,8),2); pt8z = coordinatesFEM(elementsFEM(j,8),3);
    
    % winstepsize = abs(pt1x-pt7x);
    
    % ------ Find mid pts 5/6/7/8 -------
    % DELETED;
    
    % ------ Calculate coefficients of ksi, eta and zeta ------
%     lMatrix = [ pt1x*pt1y*pt1z, pt1x*pt1y, pt1y*pt1z, pt1z*pt1x, pt1x, pt1y, pt1z, 1;
%         pt2x*pt2y*pt2z, pt2x*pt2y, pt2y*pt2z, pt2z*pt2x, pt2x, pt2y, pt2z, 1;
%         pt3x*pt3y*pt3z, pt3x*pt3y, pt3y*pt3z, pt3z*pt3x, pt3x, pt3y, pt3z, 1;
%         pt4x*pt4y*pt4z, pt4x*pt4y, pt4y*pt4z, pt4z*pt4x, pt4x, pt4y, pt4z, 1;
%         pt5x*pt5y*pt5z, pt5x*pt5y, pt5y*pt5z, pt5z*pt5x, pt5x, pt5y, pt5z, 1;
%         pt6x*pt6y*pt6z, pt6x*pt6y, pt6y*pt6z, pt6z*pt6x, pt6x, pt6y, pt6z, 1;
%         pt7x*pt7y*pt7z, pt7x*pt7y, pt7y*pt7z, pt7z*pt7x, pt7x, pt7y, pt7z, 1;
%         pt8x*pt8y*pt8z, pt8x*pt8y, pt8y*pt8z, pt8z*pt8x, pt8x, pt8y, pt8z, 1];
%     
%     lb = [-1;1;1;-1;-1;1;1;-1]; lCoeff = linsolve(lMatrix,lb);
%     mb = [-1;-1;1;1;-1;-1;1;1]; mCoeff = linsolve(lMatrix,mb);
%     nb = [-1;-1;-1;-1;1;1;1;1]; nCoeff = linsolve(lMatrix,nb);
    
    % ------ Input info for CrackPath -----
    % DELETED
    
    % ------ Find the element nodal indices ------
    tp = ones(1,DIM);
    tempIndexU = 3*elementsFEM(j,[tp,2*tp,3*tp,4*tp,5*tp,6*tp,7*tp,8*tp]);
    tempIndexU(1:3:end) = tempIndexU(1:3:end)-2;
    tempIndexU(2:3:end) = tempIndexU(2:3:end)-1;
    
    
    if GaussPtOrder > 0
        % ------ Set Gauss points ------
        pt1 = elementsFEM(j,1); pt2 = elementsFEM(j,7);
        ptOfx = zeros(length(gqwt),1); ptOfy = ptOfx; ptOfz = ptOfx;
        for tempi = 1:length(gqwt)
            ptOfx(tempi) = gqpt(tempi)*0.5*(coordinatesFEM(pt2,1)-coordinatesFEM(pt1,1))+0.5*(coordinatesFEM(pt2,1)+coordinatesFEM(pt1,1));
            ptOfy(tempi) = gqpt(tempi)*0.5*(coordinatesFEM(pt2,2)-coordinatesFEM(pt1,2))+0.5*(coordinatesFEM(pt2,2)+coordinatesFEM(pt1,2));
            ptOfz(tempi) = gqpt(tempi)*0.5*(coordinatesFEM(pt2,3)-coordinatesFEM(pt1,3))+0.5*(coordinatesFEM(pt2,3)+coordinatesFEM(pt1,3));
        end
        StrainWithinEachElementGausspoint = zeros(length(gqwt)^3,9); % each GP there are 9 F terms.
        
        for tempk = 1:size(ptOfz,1)
            for tempi = 1:length(ptOfx)
                for tempj = 1:length(ptOfy)
                    
                    % ------ Calculate ksi, eta and zeta ------
%                     ksi = lCoeff(1)*ptOfx(tempi)*ptOfy(tempj)*ptOfz(tempk) + lCoeff(2)*ptOfx(tempi)*ptOfy(tempj) + lCoeff(3)*ptOfy(tempj)*ptOfz(tempk) + ...
%                         lCoeff(4)*ptOfz(tempk)*ptOfx(tempi) + lCoeff(5)*ptOfx(tempi) + lCoeff(6)*ptOfy(tempj) + lCoeff(7)*ptOfz(tempk) + lCoeff(8);
%                     eta = mCoeff(1)*ptOfx(tempi)*ptOfy(tempj)*ptOfz(tempk) + mCoeff(2)*ptOfx(tempi)*ptOfy(tempj) + mCoeff(3)*ptOfy(tempj)*ptOfz(tempk) + ...
%                         mCoeff(4)*ptOfz(tempk)*ptOfx(tempi) + mCoeff(5)*ptOfx(tempi) + mCoeff(6)*ptOfy(tempj) + mCoeff(7)*ptOfz(tempk) + mCoeff(8);
%                     zeta = nCoeff(1)*ptOfx(tempi)*ptOfy(tempj)*ptOfz(tempk) + nCoeff(2)*ptOfx(tempi)*ptOfy(tempj) + nCoeff(3)*ptOfy(tempj)*ptOfz(tempk) + ...
%                         nCoeff(4)*ptOfz(tempk)*ptOfx(tempi) + nCoeff(5)*ptOfx(tempi) + nCoeff(6)*ptOfy(tempj) + nCoeff(7)*ptOfz(tempk) + nCoeff(8);
                    ksi = gqpt(tempi);
                    eta = gqpt(tempj);
                    zeta = gqpt(tempk);
                    
                    % ------ Calculate N ------
                    % N1 = 1/8*(1-ksi)*(1-eta)*(1-zeta); N2 = 1/8*(1+ksi)*(1-eta)*(1-zeta);
                    % N3 = 1/8*(1+ksi)*(1+eta)*(1-zeta); N4 = 1/8*(1-ksi)*(1+eta)*(1-zeta);
                    % N5 = 1/8*(1-ksi)*(1-eta)*(1+zeta); N6 = 1/8*(1+ksi)*(1-eta)*(1+zeta);
                    % N7 = 1/8*(1+ksi)*(1+eta)*(1+zeta); N8 = 1/8*(1-ksi)*(1+eta)*(1+zeta);
                    
                    % ------ Generate [N] shape function matrix ------
                    % tpN1 = diag([N1,N1,N1]); tpN2 = diag([N2,N2,N2]); tpN3 = diag([N3,N3,N3]); tpN4 = diag([N4,N4,N4]);
                    % tpN5 = diag([N5,N5,N5]); tpN6 = diag([N6,N6,N6]); tpN7 = diag([N7,N7,N7]); tpN8 = diag([N8,N8,N8]);
                    % NOrig = [tpN1,tpN2,tpN3,tpN4,tpN5,tpN6,tpN7,tpN8]; % N = NOrig;
                    
                    % NDiag = diag([N1*ones(1,DIM),N2*ones(1,DIM),N3*ones(1,DIM),N4*ones(1,DIM),N5*ones(1,DIM),N6*ones(1,DIM),N7*ones(1,DIM),N8*ones(1,DIM)]);
                    
                    
                    % ------ Build J matrix ------
                    % Comment: I didn't change Jacobian matrix J when enriched
                    % functions are added.
                    J = [funDN1Dksi(ksi,eta,zeta),funDN2Dksi(ksi,eta,zeta),funDN3Dksi(ksi,eta,zeta),funDN4Dksi(ksi,eta,zeta), ...
                        funDN5Dksi(ksi,eta,zeta),funDN6Dksi(ksi,eta,zeta),funDN7Dksi(ksi,eta,zeta),funDN8Dksi(ksi,eta,zeta);
                        funDN1Deta(ksi,eta,zeta),funDN2Deta(ksi,eta,zeta),funDN3Deta(ksi,eta,zeta),funDN4Deta(ksi,eta,zeta), ...
                        funDN5Deta(ksi,eta,zeta),funDN6Deta(ksi,eta,zeta),funDN7Deta(ksi,eta,zeta),funDN8Deta(ksi,eta,zeta);
                        funDN1Dzeta(ksi,eta,zeta),funDN2Dzeta(ksi,eta,zeta),funDN3Dzeta(ksi,eta,zeta),funDN4Dzeta(ksi,eta,zeta), ...
                        funDN5Dzeta(ksi,eta,zeta),funDN6Dzeta(ksi,eta,zeta),funDN7Dzeta(ksi,eta,zeta),funDN8Dzeta(ksi,eta,zeta)] * ...
                        [pt1x,pt1y,pt1z; pt2x,pt2y,pt2z; pt3x,pt3y,pt3z; pt4x,pt4y,pt4z; ...
                         pt5x,pt5y,pt5z; pt6x,pt6y,pt6z; pt7x,pt7y,pt7z; pt8x,pt8y,pt8z];
                    
                    % Jacobian = det(J);
                    InvJ = inv(J);
                    
                    % ------ Compute DN matrix ------
                    DNOrig = [InvJ zeros(3,3) zeros(3,3); zeros(3,3) InvJ zeros(3,3); zeros(3,3) zeros(3,3) InvJ] * ...
                        [funDN1Dksi(ksi,eta,zeta) 0 0 funDN2Dksi(ksi,eta,zeta) 0 0 funDN3Dksi(ksi,eta,zeta) 0 0 funDN4Dksi(ksi,eta,zeta) 0 0 ...
                        funDN5Dksi(ksi,eta,zeta) 0 0 funDN6Dksi(ksi,eta,zeta) 0 0 funDN7Dksi(ksi,eta,zeta) 0 0 funDN8Dksi(ksi,eta,zeta) 0 0;
                        funDN1Deta(ksi,eta,zeta) 0 0 funDN2Deta(ksi,eta,zeta) 0 0 funDN3Deta(ksi,eta,zeta) 0 0 funDN4Deta(ksi,eta,zeta) 0 0 ...
                        funDN5Deta(ksi,eta,zeta) 0 0 funDN6Deta(ksi,eta,zeta) 0 0 funDN7Deta(ksi,eta,zeta) 0 0 funDN8Deta(ksi,eta,zeta) 0 0;
                        funDN1Dzeta(ksi,eta,zeta) 0 0 funDN2Dzeta(ksi,eta,zeta) 0 0 funDN3Dzeta(ksi,eta,zeta) 0 0 funDN4Dzeta(ksi,eta,zeta) 0 0 ...
                        funDN5Dzeta(ksi,eta,zeta) 0 0 funDN6Dzeta(ksi,eta,zeta) 0 0 funDN7Dzeta(ksi,eta,zeta) 0 0 funDN8Dzeta(ksi,eta,zeta) 0 0;
                        0 funDN1Dksi(ksi,eta,zeta) 0 0 funDN2Dksi(ksi,eta,zeta) 0 0 funDN3Dksi(ksi,eta,zeta) 0 0 funDN4Dksi(ksi,eta,zeta) 0 ...
                        0 funDN5Dksi(ksi,eta,zeta) 0 0 funDN6Dksi(ksi,eta,zeta) 0 0 funDN7Dksi(ksi,eta,zeta) 0 0 funDN8Dksi(ksi,eta,zeta) 0 ;
                        0 funDN1Deta(ksi,eta,zeta) 0 0 funDN2Deta(ksi,eta,zeta) 0 0 funDN3Deta(ksi,eta,zeta) 0 0 funDN4Deta(ksi,eta,zeta) 0 ...
                        0 funDN5Deta(ksi,eta,zeta) 0 0 funDN6Deta(ksi,eta,zeta) 0 0 funDN7Deta(ksi,eta,zeta) 0 0 funDN8Deta(ksi,eta,zeta) 0 ;
                        0 funDN1Dzeta(ksi,eta,zeta) 0 0 funDN2Dzeta(ksi,eta,zeta) 0 0 funDN3Dzeta(ksi,eta,zeta) 0 0 funDN4Dzeta(ksi,eta,zeta) 0 ...
                        0 funDN5Dzeta(ksi,eta,zeta) 0 0 funDN6Dzeta(ksi,eta,zeta) 0 0 funDN7Dzeta(ksi,eta,zeta) 0 0 funDN8Dzeta(ksi,eta,zeta) 0 ;
                        0 0 funDN1Dksi(ksi,eta,zeta) 0 0 funDN2Dksi(ksi,eta,zeta) 0 0 funDN3Dksi(ksi,eta,zeta) 0 0 funDN4Dksi(ksi,eta,zeta) ...
                        0 0 funDN5Dksi(ksi,eta,zeta) 0 0 funDN6Dksi(ksi,eta,zeta) 0 0 funDN7Dksi(ksi,eta,zeta) 0 0 funDN8Dksi(ksi,eta,zeta) ;
                        0 0 funDN1Deta(ksi,eta,zeta) 0 0 funDN2Deta(ksi,eta,zeta) 0 0 funDN3Deta(ksi,eta,zeta) 0 0 funDN4Deta(ksi,eta,zeta) ...
                        0 0 funDN5Deta(ksi,eta,zeta) 0 0 funDN6Deta(ksi,eta,zeta) 0 0 funDN7Deta(ksi,eta,zeta) 0 0 funDN8Deta(ksi,eta,zeta) ;
                        0 0 funDN1Dzeta(ksi,eta,zeta) 0 0 funDN2Dzeta(ksi,eta,zeta) 0 0 funDN3Dzeta(ksi,eta,zeta) 0 0 funDN4Dzeta(ksi,eta,zeta) ...
                        0 0 funDN5Dzeta(ksi,eta,zeta) 0 0 funDN6Dzeta(ksi,eta,zeta) 0 0 funDN7Dzeta(ksi,eta,zeta) 0 0 funDN8Dzeta(ksi,eta,zeta)];
                    DN = DNOrig;
                    
                    StrainWithinEachElementGausspoint(length(gqpt)^2*(tempk-1)+length(gqpt)*(tempi-1)+tempj,1:9) = DN*U(tempIndexU);
                    
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%% Replace with 9-Gauss point %%%%%%%%%%%%%%%%%%%%%%
        switch GaussPtOrder
            case 5
                StrainGaussPt((length(gqpt)^3*(j-1)+1):length(gqpt)^3*j,1:9) = StrainWithinEachElementGausspoint;
                CoordsGaussPt((length(gqpt)^3*(j-1)+1):length(gqpt)^3*j,1:3) = ...
                    [ptOfx(1),ptOfy(1),ptOfz(1);ptOfx(1),ptOfy(2),ptOfz(1);ptOfx(1),ptOfy(3),ptOfz(1);ptOfx(1),ptOfy(4),ptOfz(1);ptOfx(1),ptOfy(5),ptOfz(1);
                    ptOfx(2),ptOfy(1),ptOfz(1);ptOfx(2),ptOfy(2),ptOfz(1);ptOfx(2),ptOfy(3),ptOfz(1);ptOfx(2),ptOfy(4),ptOfz(1);ptOfx(2),ptOfy(5),ptOfz(1);
                    ptOfx(3),ptOfy(1),ptOfz(1);ptOfx(3),ptOfy(2),ptOfz(1);ptOfx(3),ptOfy(3),ptOfz(1);ptOfx(3),ptOfy(4),ptOfz(1);ptOfx(3),ptOfy(5),ptOfz(1);
                    ptOfx(4),ptOfy(1),ptOfz(1);ptOfx(4),ptOfy(2),ptOfz(1);ptOfx(4),ptOfy(3),ptOfz(1);ptOfx(4),ptOfy(4),ptOfz(1);ptOfx(4),ptOfy(5),ptOfz(1);
                    ptOfx(5),ptOfy(1),ptOfz(1);ptOfx(5),ptOfy(2),ptOfz(1);ptOfx(5),ptOfy(3),ptOfz(1);ptOfx(5),ptOfy(4),ptOfz(1);ptOfx(5),ptOfy(5),ptOfz(1);
                    ptOfx(1),ptOfy(1),ptOfz(2);ptOfx(1),ptOfy(2),ptOfz(2);ptOfx(1),ptOfy(3),ptOfz(2);ptOfx(1),ptOfy(4),ptOfz(2);ptOfx(1),ptOfy(5),ptOfz(2);
                    ptOfx(2),ptOfy(1),ptOfz(2);ptOfx(2),ptOfy(2),ptOfz(2);ptOfx(2),ptOfy(3),ptOfz(2);ptOfx(2),ptOfy(4),ptOfz(2);ptOfx(2),ptOfy(5),ptOfz(2);
                    ptOfx(3),ptOfy(1),ptOfz(2);ptOfx(3),ptOfy(2),ptOfz(2);ptOfx(3),ptOfy(3),ptOfz(2);ptOfx(3),ptOfy(4),ptOfz(2);ptOfx(3),ptOfy(5),ptOfz(2);
                    ptOfx(4),ptOfy(1),ptOfz(2);ptOfx(4),ptOfy(2),ptOfz(2);ptOfx(4),ptOfy(3),ptOfz(2);ptOfx(4),ptOfy(4),ptOfz(2);ptOfx(4),ptOfy(5),ptOfz(2);
                    ptOfx(5),ptOfy(1),ptOfz(2);ptOfx(5),ptOfy(2),ptOfz(2);ptOfx(5),ptOfy(3),ptOfz(2);ptOfx(5),ptOfy(4),ptOfz(2);ptOfx(5),ptOfy(5),ptOfz(2);
                    ptOfx(1),ptOfy(1),ptOfz(3);ptOfx(1),ptOfy(2),ptOfz(3);ptOfx(1),ptOfy(3),ptOfz(3);ptOfx(1),ptOfy(4),ptOfz(3);ptOfx(1),ptOfy(5),ptOfz(3);
                    ptOfx(2),ptOfy(1),ptOfz(3);ptOfx(2),ptOfy(2),ptOfz(3);ptOfx(2),ptOfy(3),ptOfz(3);ptOfx(2),ptOfy(4),ptOfz(3);ptOfx(2),ptOfy(5),ptOfz(3);
                    ptOfx(3),ptOfy(1),ptOfz(3);ptOfx(3),ptOfy(2),ptOfz(3);ptOfx(3),ptOfy(3),ptOfz(3);ptOfx(3),ptOfy(4),ptOfz(3);ptOfx(3),ptOfy(5),ptOfz(3);
                    ptOfx(4),ptOfy(1),ptOfz(3);ptOfx(4),ptOfy(2),ptOfz(3);ptOfx(4),ptOfy(3),ptOfz(3);ptOfx(4),ptOfy(4),ptOfz(3);ptOfx(4),ptOfy(5),ptOfz(3);
                    ptOfx(5),ptOfy(1),ptOfz(3);ptOfx(5),ptOfy(2),ptOfz(3);ptOfx(5),ptOfy(3),ptOfz(3);ptOfx(5),ptOfy(4),ptOfz(3);ptOfx(5),ptOfy(5),ptOfz(3);
                    ptOfx(1),ptOfy(1),ptOfz(4);ptOfx(1),ptOfy(2),ptOfz(4);ptOfx(1),ptOfy(3),ptOfz(4);ptOfx(1),ptOfy(4),ptOfz(4);ptOfx(1),ptOfy(5),ptOfz(4);
                    ptOfx(2),ptOfy(1),ptOfz(4);ptOfx(2),ptOfy(2),ptOfz(4);ptOfx(2),ptOfy(3),ptOfz(4);ptOfx(2),ptOfy(4),ptOfz(4);ptOfx(2),ptOfy(5),ptOfz(4);
                    ptOfx(3),ptOfy(1),ptOfz(4);ptOfx(3),ptOfy(2),ptOfz(4);ptOfx(3),ptOfy(3),ptOfz(4);ptOfx(3),ptOfy(4),ptOfz(4);ptOfx(3),ptOfy(5),ptOfz(4);
                    ptOfx(4),ptOfy(1),ptOfz(4);ptOfx(4),ptOfy(2),ptOfz(4);ptOfx(4),ptOfy(3),ptOfz(4);ptOfx(4),ptOfy(4),ptOfz(4);ptOfx(4),ptOfy(5),ptOfz(4);
                    ptOfx(5),ptOfy(1),ptOfz(4);ptOfx(5),ptOfy(2),ptOfz(4);ptOfx(5),ptOfy(3),ptOfz(4);ptOfx(5),ptOfy(4),ptOfz(4);ptOfx(5),ptOfy(5),ptOfz(4);
                    ptOfx(1),ptOfy(1),ptOfz(5);ptOfx(1),ptOfy(2),ptOfz(5);ptOfx(1),ptOfy(3),ptOfz(5);ptOfx(1),ptOfy(4),ptOfz(5);ptOfx(1),ptOfy(5),ptOfz(5);
                    ptOfx(2),ptOfy(1),ptOfz(5);ptOfx(2),ptOfy(2),ptOfz(5);ptOfx(2),ptOfy(3),ptOfz(5);ptOfx(2),ptOfy(4),ptOfz(5);ptOfx(2),ptOfy(5),ptOfz(5);
                    ptOfx(3),ptOfy(1),ptOfz(5);ptOfx(3),ptOfy(2),ptOfz(5);ptOfx(3),ptOfy(3),ptOfz(5);ptOfx(3),ptOfy(4),ptOfz(5);ptOfx(3),ptOfy(5),ptOfz(5);
                    ptOfx(4),ptOfy(1),ptOfz(5);ptOfx(4),ptOfy(2),ptOfz(5);ptOfx(4),ptOfy(3),ptOfz(5);ptOfx(4),ptOfy(4),ptOfz(5);ptOfx(4),ptOfy(5),ptOfz(5);
                    ptOfx(5),ptOfy(1),ptOfz(5);ptOfx(5),ptOfy(2),ptOfz(5);ptOfx(5),ptOfy(3),ptOfz(5);ptOfx(5),ptOfy(4),ptOfz(5);ptOfx(5),ptOfy(5),ptOfz(5)];
                %%%%%%%%%%%%%%%%%%% Comment following 4-Gauss point %%%%%%%%%%%%%%%%%%%%%%
            case 3
                StrainGaussPt((length(gqpt)^3*(j-1)+1):length(gqpt)^3*j,1:9) = StrainWithinEachElementGausspoint;
                CoordsGaussPt((length(gqpt)^3*(j-1)+1):length(gqpt)^3*j,1:3) = ...
                    [ptOfx(1),ptOfy(1),ptOfz(1);ptOfx(1),ptOfy(2),ptOfz(1);ptOfx(1),ptOfy(3),ptOfz(1); 
                    ptOfx(2),ptOfy(1),ptOfz(1);ptOfx(2),ptOfy(2),ptOfz(1);ptOfx(2),ptOfy(3),ptOfz(1);
                    ptOfx(3),ptOfy(1),ptOfz(1);ptOfx(3),ptOfy(2),ptOfz(1);ptOfx(3),ptOfy(3),ptOfz(1);
                    ptOfx(1),ptOfy(1),ptOfz(2);ptOfx(1),ptOfy(2),ptOfz(2);ptOfx(1),ptOfy(3),ptOfz(2); 
                    ptOfx(2),ptOfy(1),ptOfz(2);ptOfx(2),ptOfy(2),ptOfz(2);ptOfx(2),ptOfy(3),ptOfz(2);
                    ptOfx(3),ptOfy(1),ptOfz(2);ptOfx(3),ptOfy(2),ptOfz(2);ptOfx(3),ptOfy(3),ptOfz(2);
                    ptOfx(1),ptOfy(1),ptOfz(3);ptOfx(1),ptOfy(2),ptOfz(3);ptOfx(1),ptOfy(3),ptOfz(3); 
                    ptOfx(2),ptOfy(1),ptOfz(3);ptOfx(2),ptOfy(2),ptOfz(3);ptOfx(2),ptOfy(3),ptOfz(3);
                    ptOfx(3),ptOfy(1),ptOfz(3);ptOfx(3),ptOfy(2),ptOfz(3);ptOfx(3),ptOfy(3),ptOfz(3)];
                %%%%%%%%%%%%%%%%%%% Comment following 4-Gauss point %%%%%%%%%%%%%%%%%%%%%%
            case 2
                StrainGaussPt((length(gqpt)^3*(j-1)+1):length(gqpt)^3*j,1:9) = StrainWithinEachElementGausspoint;
                CoordsGaussPt((length(gqpt)^3*(j-1)+1):length(gqpt)^3*j,1:3) = ...
                    [ptOfx(1),ptOfy(1),ptOfz(1);ptOfx(1),ptOfy(2),ptOfz(1);
                    ptOfx(2),ptOfy(1),ptOfz(1);ptOfx(2),ptOfy(2),ptOfz(1);
                    ptOfx(1),ptOfy(1),ptOfz(2);ptOfx(1),ptOfy(2),ptOfz(2);
                    ptOfx(2),ptOfy(1),ptOfz(2);ptOfx(2),ptOfy(2),ptOfz(2)];
                
                MatrixExtrapolation = [2.54904, -0.683013, 0.183013, -0.683013, -0.683013, 0.183013, -0.0490381, 0.183013;
                                        -0.683013, 2.54904, -0.683013, 0.183013, 0.183013, -0.683013, 0.183013, -0.0490381;
                                        0.183013, -0.683013, 2.54904, -0.683013, -0.0490381, 0.183013, -0.683013, 0.183013;
                                        -0.683013, 0.183013, -0.683013, 2.54904, 0.183013, -0.0490381, 0.183013, -0.683013;
                                        -0.683013, 0.183013, -0.0490381, 0.183013, 2.54904, -0.683013, 0.183013, -0.683013;
                                        0.183013, -0.683013, 0.183013, -0.0490381, -0.683013, 2.54904, -0.683013, 0.183013;
                                        -0.0490381, 0.183013, -0.683013, 0.183013, 0.183013, -0.683013, 2.54904, -0.683013;
                                        0.183013, -0.0490381, 0.183013, -0.683013, -0.683013, 0.183013, -0.683013, 2.54904];
               % ------ Nodal points strain extrapolation using Gauss points -----
               StrainWithinEachElementNodalpoint =  (MatrixExtrapolation * StrainWithinEachElementGausspoint)';
               
               % ------ Find the element nodal indices for strain ------
               tempStrainIndex = [9*elementsFEM(j,1)*ones(1,9)-[8:-1:0], 9*elementsFEM(j,2)*ones(1,9)-[8:-1:0], ...
                   9*elementsFEM(j,3)*ones(1,9)-[8:-1:0], 9*elementsFEM(j,4)*ones(1,9)-[8:-1:0], ...
                   9*elementsFEM(j,5)*ones(1,9)-[8:-1:0], 9*elementsFEM(j,6)*ones(1,9)-[8:-1:0], ...
                   9*elementsFEM(j,7)*ones(1,9)-[8:-1:0], 9*elementsFEM(j,8)*ones(1,9)-[8:-1:0]];
                 
               FStrain(tempStrainIndex) = FStrain(tempStrainIndex) + StrainWithinEachElementNodalpoint(:);
               FStrainAvgTimes(tempStrainIndex) = FStrainAvgTimes(tempStrainIndex) + ones(9*8,1);

        end
        %%%%%%%%%%%%%%%%%%%%%%%% End of Comment %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
end

switch GaussPtOrder
    case 2
        F = FStrain./FStrainAvgTimes; F=F(:);
    otherwise
        F = zeros(size(coordinatesFEM,1)*9,1); 
end

% % ------ get nodal strains ------
% Coordxnodes = unique(coordinatesFEM(:,1)); Coordynodes = unique(coordinatesFEM(:,2)); Coordznodes = unique(coordinatesFEM(:,3));
% [Xq,Yq,Zq] = ndgrid(Coordxnodes,Coordynodes,Coordznodes);
% M = length(Coordxnodes);  N = length(Coordynodes);  L = length(Coordznodes);
%  
% %Expensive, but for arbitrary mesh.
% F11temp = scatteredInterpolant(CoordsGaussPt(:,1),CoordsGaussPt(:,2),CoordsGaussPt(:,3),StrainGaussPt(:,1),'linear','linear');
% F21temp = scatteredInterpolant(CoordsGaussPt(:,1),CoordsGaussPt(:,2),CoordsGaussPt(:,3),StrainGaussPt(:,2),'linear','linear');
% F31temp = scatteredInterpolant(CoordsGaussPt(:,1),CoordsGaussPt(:,2),CoordsGaussPt(:,3),StrainGaussPt(:,3),'linear','linear');
% F12temp = scatteredInterpolant(CoordsGaussPt(:,1),CoordsGaussPt(:,2),CoordsGaussPt(:,3),StrainGaussPt(:,4),'linear','linear');
% F22temp = scatteredInterpolant(CoordsGaussPt(:,1),CoordsGaussPt(:,2),CoordsGaussPt(:,3),StrainGaussPt(:,5),'linear','linear');
% F32temp = scatteredInterpolant(CoordsGaussPt(:,1),CoordsGaussPt(:,2),CoordsGaussPt(:,3),StrainGaussPt(:,6),'linear','linear');
% F13temp = scatteredInterpolant(CoordsGaussPt(:,1),CoordsGaussPt(:,2),CoordsGaussPt(:,3),StrainGaussPt(:,7),'linear','linear');
% F23temp = scatteredInterpolant(CoordsGaussPt(:,1),CoordsGaussPt(:,2),CoordsGaussPt(:,3),StrainGaussPt(:,8),'linear','linear');
% F33temp = scatteredInterpolant(CoordsGaussPt(:,1),CoordsGaussPt(:,2),CoordsGaussPt(:,3),StrainGaussPt(:,9),'linear','linear');
% 
% F11 = F11temp(Xq,Yq,Zq); F21 = F21temp(Xq,Yq,Zq); F31 = F31temp(Xq,Yq,Zq);
% F12 = F12temp(Xq,Yq,Zq); F22 = F22temp(Xq,Yq,Zq); F32 = F32temp(Xq,Yq,Zq);
% F13 = F13temp(Xq,Yq,Zq); F23 = F23temp(Xq,Yq,Zq); F33 = F33temp(Xq,Yq,Zq);
%  
% 
% F = zeros(9*size(coordinatesFEM,1),1);
% for tempi = 1:size(coordinatesFEM,1)
%     
%     [row1,col1] = find(Coordxnodes==coordinatesFEM(tempi,1));
%     [row2,col2] = find(Coordynodes==coordinatesFEM(tempi,2));
%     [row3,col3] = find(Coordznodes==coordinatesFEM(tempi,3));
%      
%     F(9*tempi-8) = F11(row1,row2,row3); F(9*tempi-7) = F21(row1,row2,row3); F(9*tempi-6) = F31(row1,row2,row3);
%     F(9*tempi-5) = F12(row1,row2,row3); F(9*tempi-4) = F22(row1,row2,row3); F(9*tempi-3) = F32(row1,row2,row3);
%     F(9*tempi-2) = F13(row1,row2,row3); F(9*tempi-1) = F23(row1,row2,row3); F(9*tempi-0) = F33(row1,row2,row3);  
% end



%% plot Gaussian points strains

% dvcZOI = createpde(1);
% % Apply mesh
% DT = delaunayTriangulation(CoordsGaussPt(:,1), CoordsGaussPt(:,2), CoordsGaussPt(:,3)); 
% geometryFromMesh(dvcZOI,DT.Points',DT.ConnectivityList');
% figure, pdeplot3D(dvcZOI,'ColorMapData',StrainGaussPt(:,1),'FaceAlpha',0.5);
% figure, pdegplot(dvcZOI,'FaceLabels','on','FaceAlpha',0.5); set(gca,'fontsize',18);title('ZOI body')


end

%% ========= subroutines for  FEM shape function derivatives ========
function a = funDN1Dksi(ksi,eta,zeta)
a = 1/8*(-1)*(1-eta)*(1-zeta);
end
function a = funDN2Dksi(ksi,eta,zeta)
a = 1/8*( 1)*(1-eta)*(1-zeta);
end
function a = funDN3Dksi(ksi,eta,zeta)
a = 1/8*( 1)*(1+eta)*(1-zeta);
end
function a = funDN4Dksi(ksi,eta,zeta)
a = 1/8*(-1)*(1+eta)*(1-zeta);
end
function a = funDN5Dksi(ksi,eta,zeta)
a = 1/8*(-1)*(1-eta)*(1+zeta);
end
function a = funDN6Dksi(ksi,eta,zeta)
a = 1/8*( 1)*(1-eta)*(1+zeta);
end
function a = funDN7Dksi(ksi,eta,zeta)
a = 1/8*( 1)*(1+eta)*(1+zeta);
end
function a = funDN8Dksi(ksi,eta,zeta)
a = 1/8*(-1)*(1+eta)*(1+zeta);
end

% ----------------------------------------------------
function a = funDN1Deta(ksi,eta,zeta)
a = 1/8*(1-ksi)*(-1)*(1-zeta);
end
function a = funDN2Deta(ksi,eta,zeta)
a = 1/8*(1+ksi)*(-1)*(1-zeta);
end
function a = funDN3Deta(ksi,eta,zeta)
a = 1/8*(1+ksi)*( 1)*(1-zeta);
end
function a = funDN4Deta(ksi,eta,zeta)
a = 1/8*(1-ksi)*( 1)*(1-zeta);
end
function a = funDN5Deta(ksi,eta,zeta)
a = 1/8*(1-ksi)*(-1)*(1+zeta);
end
function a = funDN6Deta(ksi,eta,zeta)
a = 1/8*(1+ksi)*(-1)*(1+zeta);
end
function a = funDN7Deta(ksi,eta,zeta)
a = 1/8*(1+ksi)*( 1)*(1+zeta);
end
function a = funDN8Deta(ksi,eta,zeta)
a = 1/8*(1-ksi)*( 1)*(1+zeta);
end

% ----------------------------------------------------
function a = funDN1Dzeta(ksi,eta,zeta)
a = 1/8*(1-ksi)*(1-eta)*(-1);
end
function a = funDN2Dzeta(ksi,eta,zeta)
a = 1/8*(1+ksi)*(1-eta)*(-1);
end
function a = funDN3Dzeta(ksi,eta,zeta)
a = 1/8*(1+ksi)*(1+eta)*(-1);
end
function a = funDN4Dzeta(ksi,eta,zeta)
a = 1/8*(1-ksi)*(1+eta)*(-1);
end
function a = funDN5Dzeta(ksi,eta,zeta)
a = 1/8*(1-ksi)*(1-eta)*( 1);
end
function a = funDN6Dzeta(ksi,eta,zeta)
a = 1/8*(1+ksi)*(1-eta)*( 1);
end
function a = funDN7Dzeta(ksi,eta,zeta)
a = 1/8*(1+ksi)*(1+eta)*( 1);
end
function a = funDN8Dzeta(ksi,eta,zeta)
a = 1/8*(1-ksi)*(1+eta)*( 1);
end




