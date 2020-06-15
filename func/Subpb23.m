%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function AL-DIC Subproblem 2  3d case                    %
% Object: to find deformation field using global methods   %
% Author: Jin Yang                                         %
% Last date modified: 2019.03                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Uhat] = Subpb23(DVCmesh, beta, mu, U, F, W, v, alpha, GaussPtOrder)

coordinatesFEM = DVCmesh.coordinatesFEM;
elementsFEM = DVCmesh.elementsFEM;
dirichlet = DVCmesh.dirichlet;
neumann = DVCmesh.neumann;

DIM = 3; NodesPerEle = 8; alpha2=0; % Ignore alpha2 first for this code.
FEMSize = size(coordinatesFEM,1);
  
% ====== Info of CrackTip path-line ======
% k_tip = -0.5*(CrackPath1(1)/CrackPath1(2)+CrackPath2(1)/CrackPath2(2));
% CrackPathCen = [k_tip/(CrackTip(2)-k_tip*CrackTip(1)), -1/(CrackTip(2)-k_tip*CrackTip(1))];
% CrackPathCen = [0,0];

% ------ Pre-Subpb2 deal with cracks ------
% EnrHAndTipEleIndex is the collection of index for Heaviside function basis;
% EnrTipEleIndex is the collection of index for Crack tip singular basis;
% if CrackOrNot == 0
% EnrHAndTipEleIndex = []; EnrTipEleIndex = [];

% else
%     [EnrHAndTipEleIndex,EnrTipEleIndex] = EnrEleWithFrac(coordinatesFEM,elementsFEM,CrackPath1,CrackPath2,CrackTip);
%     if CrackTipOrNot == 1
%         FEMSize = size(coordinatesFEM,1) + size(EnrHAndTipEleIndex,1) + 4*size(EnrTipEleIndex,1);
%         udual = 0*FSubpb1; vdual = 0*USubpb1;
%         F = [FSubpb1;zeros(4*size(EnrHAndTipEleIndex,1),1);zeros(4*4*size(EnrTipEleIndex,1),1)];
%         U = [USubpb1;zeros(2*size(EnrHAndTipEleIndex,1),1);zeros(2*4*size(EnrTipEleIndex,1),1)];
%         W = [udual;zeros(4*size(EnrHAndTipEleIndex,1),1);zeros(4*4*size(EnrTipEleIndex,1),1)];
%         v = [vdual;zeros(2*size(EnrHAndTipEleIndex,1),1);zeros(2*4*size(EnrTipEleIndex,1),1)];
%     else
%         FEMSize = size(coordinatesFEM,1) + size(EnrHAndTipEleIndex,1);
%         udual = 0*FSubpb1; vdual = 0*USubpb1;
%         F = [FSubpb1;zeros(4*size(EnrHAndTipEleIndex,1),1) ];
%         U = [USubpb1;zeros(2*size(EnrHAndTipEleIndex,1),1) ];
%         W = [udual;zeros(4*size(EnrHAndTipEleIndex,1),1) ];
%         v = [vdual;zeros(2*size(EnrHAndTipEleIndex,1),1) ];
%     end
% end

% ====== Initialize variables ======
Uhat = U; U = [U;zeros(NodesPerEle*DIM,1)]; v = [v;zeros(NodesPerEle*DIM,1)];
F = [F;zeros(NodesPerEle*DIM^2,1)]; W = [W;zeros(NodesPerEle*DIM^2,1)];
% FMinusW1 = F(1:2:end)-W(1:2:end); FMinusW2 = F(2:2:end)-W(2:2:end); % Comment old codes
UMinusv = U-v; FMinusW = F-W;
% CrackTipOrNot = 1-CrackTipOrNot; % Keep consistent with old version codes

% ====== Initialize A matrix and b vector ======
% A = sparse(DIM*FEMSize+NodesPerEle*DIM, DIM*FEMSize+NodesPerEle*DIM);
% NodesPerEle*DIM is because I put all the zeros in elementsFEM to the end

% ====== Gaussian quadrature parameter ======
switch GaussPtOrder
    case 2 % ------ 2 Gauss points ------
        gqpt1 = -1/sqrt(3); gqpt2 = 1/sqrt(3); gqpt = [gqpt1,gqpt2];
        gqwt1 = 1; gqwt2 = 1; gqwt = [gqwt1,gqwt2];
    case 3 % ------ 3 Gauss points ------
        gqpt1 = 0; gqpt2 = 0.57735; gqpt3 = -0.57735; gqpt = [gqpt1,gqpt2,gqpt3];
        gqwt1 = 1; gqwt2 = 1; gqwt3 = 1; gqwt = [gqwt1,gqwt2,gqwt3];
    case 4 % ------ 4 Gaussian pts ------
        gqpt1 = 0.339981; gqpt2 = -0.339981; gqpt3 = 0.861136; gqpt4 = -0.861136;
        gqwt1 = 0.652145; gqwt2 = 0.652145; gqwt3 = 0.347855; gqwt4 = 0.347855;
        gqpt = [gqpt1,gqpt2,gqpt3,gqpt4]; gqwt = [gqwt1,gqwt2,gqwt3,gqwt4];
    case 5 % ------ 5 Gauss points ------
        gqpt1 = 0; gqpt2 = 0.538469; gqpt3 = -0.538469; gqpt4 = 0.90618; gqpt5 = -0.90618;
        gqwt1 = 0.568889; gqwt2 = 0.478629; gqwt3 = 0.478629; gqwt4 = 0.236927; gqwt5 = 0.236927;
        gqpt = [gqpt1,gqpt2,gqpt3,gqpt4,gqpt5]; gqwt = [gqwt1,gqwt2,gqwt3,gqwt4,gqwt5];
    otherwise
        disp('Incorrect GaussPtOrder')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ====== Initialize Global FEM solver ======
ConvergeOrNot = 0; IterStep = 0;
while ConvergeOrNot < 0.5 && IterStep < 1
    
    % ====== Initialize A matrix and b vector ======
    % clear b; b = sparse(DIM*FEMSize+NodesPerEle*DIM,1); 
    IterStep = IterStep+1;
    
    INDEXAI = []; INDEXAJ = []; INDEXAVAL = []; INDEXBI = []; INDEXBVAL = [];
        
    % ============= Each element, assemble stiffness matrix ============
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
        
        % ------ Find mid pts 5/6/7/8 -------
        % DELETED;
        
        % ------ Calculate coefficients of ksi, eta and zeta ------
        % lMatrix = [ pt1x*pt1y*pt1z, pt1x*pt1y, pt1y*pt1z, pt1z*pt1x, pt1x, pt1y, pt1z, 1;
        %             pt2x*pt2y*pt2z, pt2x*pt2y, pt2y*pt2z, pt2z*pt2x, pt2x, pt2y, pt2z, 1;
        %             pt3x*pt3y*pt3z, pt3x*pt3y, pt3y*pt3z, pt3z*pt3x, pt3x, pt3y, pt3z, 1;
        %             pt4x*pt4y*pt4z, pt4x*pt4y, pt4y*pt4z, pt4z*pt4x, pt4x, pt4y, pt4z, 1;
        %             pt5x*pt5y*pt5z, pt5x*pt5y, pt5y*pt5z, pt5z*pt5x, pt5x, pt5y, pt5z, 1;
        %             pt6x*pt6y*pt6z, pt6x*pt6y, pt6y*pt6z, pt6z*pt6x, pt6x, pt6y, pt6z, 1;
        %             pt7x*pt7y*pt7z, pt7x*pt7y, pt7y*pt7z, pt7z*pt7x, pt7x, pt7y, pt7z, 1;
        %             pt8x*pt8y*pt8z, pt8x*pt8y, pt8y*pt8z, pt8z*pt8x, pt8x, pt8y, pt8z, 1];
        % 
        % lb = [-1;1;1;-1;-1;1;1;-1]; lCoeff = linsolve(lMatrix,lb);
        % mb = [-1;-1;1;1;-1;-1;1;1]; mCoeff = linsolve(lMatrix,mb);
        % nb = [-1;-1;-1;-1;1;1;1;1]; nCoeff = linsolve(lMatrix,nb);
        
        % ------ Input info for CrackPath -----
        % DELETED
        
        % ------ Find the element nodal indices ------
        tp = ones(1,DIM);
        tempIndexU = 3*elementsFEM(j,[tp,2*tp,3*tp,4*tp,5*tp,6*tp,7*tp,8*tp]);
        tempIndexU(1:3:end) = tempIndexU(1:3:end)-2;
        tempIndexU(2:3:end) = tempIndexU(2:3:end)-1; % tempIndexU: 1*24
        
        tp = [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8];
        tempIndexF = [9*elementsFEM(j,tp)'-8, 9*elementsFEM(j,tp)'-5, 9*elementsFEM(j,tp)'-2, ...
            9*elementsFEM(j,tp)'-7, 9*elementsFEM(j,tp)'-4, 9*elementsFEM(j,tp)'-1, ...
            9*elementsFEM(j,tp)'-6, 9*elementsFEM(j,tp)'-3, 9*elementsFEM(j,tp)']; % tempIndexU: 24*9
        
        % Still keep original tempIndexU and rename it as tempIndexUOrig
        % tempIndexUOrig = tempIndexU;
        
        % ------ Find the enriched functions nodal indices -------
        % DELETED part.
        
        % ------ Set Gauss pts ------
        % pt1 = elementsFEM(j,1); pt2 = elementsFEM(j,7);
        % ptOfx = zeros(length(gqwt),1); ptOfy = ptOfx; ptOfz = ptOfx;
        % for tempi = 1:length(gqwt)
        %     ptOfx(tempi) = gqpt(tempi)*0.5*(coordinatesFEM(pt2,1)-coordinatesFEM(pt1,1))+0.5*(coordinatesFEM(pt2,1)+coordinatesFEM(pt1,1));
        %     ptOfy(tempi) = gqpt(tempi)*0.5*(coordinatesFEM(pt2,2)-coordinatesFEM(pt1,2))+0.5*(coordinatesFEM(pt2,2)+coordinatesFEM(pt1,2));
        %     ptOfz(tempi) = gqpt(tempi)*0.5*(coordinatesFEM(pt2,3)-coordinatesFEM(pt1,3))+0.5*(coordinatesFEM(pt2,3)+coordinatesFEM(pt1,3));
        % end
        
        tempA = zeros(24,24);  tempb = zeros(24,1);
        
        for tempk = 1:length(gqwt) % size(ptOfz,1)
            for tempi = 1:length(gqwt) %size(ptOfx,1)
                for tempj = 1:length(gqwt) %size(ptOfy,1)
                    
                    % ------ Calculate ksi, eta and zeta ------
                    % ksi = lCoeff(1)*ptOfx(tempi)*ptOfy(tempj)*ptOfz(tempk) + lCoeff(2)*ptOfx(tempi)*ptOfy(tempj) + lCoeff(3)*ptOfy(tempj)*ptOfz(tempk) + ...
                    %     lCoeff(4)*ptOfz(tempk)*ptOfx(tempi) + lCoeff(5)*ptOfx(tempi) + lCoeff(6)*ptOfy(tempj) + lCoeff(7)*ptOfz(tempk) + lCoeff(8);
                    % eta = mCoeff(1)*ptOfx(tempi)*ptOfy(tempj)*ptOfz(tempk) + mCoeff(2)*ptOfx(tempi)*ptOfy(tempj) + mCoeff(3)*ptOfy(tempj)*ptOfz(tempk) + ...
                    %     mCoeff(4)*ptOfz(tempk)*ptOfx(tempi) + mCoeff(5)*ptOfx(tempi) + mCoeff(6)*ptOfy(tempj) + mCoeff(7)*ptOfz(tempk) + mCoeff(8);
                    % zeta = nCoeff(1)*ptOfx(tempi)*ptOfy(tempj)*ptOfz(tempk) + nCoeff(2)*ptOfx(tempi)*ptOfy(tempj) + nCoeff(3)*ptOfy(tempj)*ptOfz(tempk) + ...
                    %     nCoeff(4)*ptOfz(tempk)*ptOfx(tempi) + nCoeff(5)*ptOfx(tempi) + nCoeff(6)*ptOfy(tempj) + nCoeff(7)*ptOfz(tempk) + nCoeff(8);
                    ksi = gqpt(tempi); eta = gqpt(tempj); zeta = gqpt(tempk);
                    
                    % ------ Calculate N ------
                    N1 = 1/8*(1-ksi)*(1-eta)*(1-zeta); N2 = 1/8*(1+ksi)*(1-eta)*(1-zeta);
                    N3 = 1/8*(1+ksi)*(1+eta)*(1-zeta); N4 = 1/8*(1-ksi)*(1+eta)*(1-zeta);
                    N5 = 1/8*(1-ksi)*(1-eta)*(1+zeta); N6 = 1/8*(1+ksi)*(1-eta)*(1+zeta);
                    N7 = 1/8*(1+ksi)*(1+eta)*(1+zeta); N8 = 1/8*(1-ksi)*(1+eta)*(1+zeta);
                    
                    % ------ Generate [N] shape function matrix ------
                    tpN1 = diag([N1,N1,N1]); tpN2 = diag([N2,N2,N2]); tpN3 = diag([N3,N3,N3]); tpN4 = diag([N4,N4,N4]);
                    tpN5 = diag([N5,N5,N5]); tpN6 = diag([N6,N6,N6]); tpN7 = diag([N7,N7,N7]); tpN8 = diag([N8,N8,N8]);
                    NOrig = [tpN1,tpN2,tpN3,tpN4,tpN5,tpN6,tpN7,tpN8]; N = NOrig;
                    NDiag = diag([N1*ones(1,DIM),N2*ones(1,DIM),N3*ones(1,DIM),N4*ones(1,DIM),N5*ones(1,DIM),N6*ones(1,DIM),N7*ones(1,DIM),N8*ones(1,DIM)]);
                    
                    % ------ Build J matrix ------
                    % Comment: I didn't change Jacobian matrix J when enriched
                    % functions are added.
                    J = [funDN1Dksi(ksi,eta,zeta),funDN2Dksi(ksi,eta,zeta),funDN3Dksi(ksi,eta,zeta),funDN4Dksi(ksi,eta,zeta), ...
                        funDN5Dksi(ksi,eta,zeta),funDN6Dksi(ksi,eta,zeta),funDN7Dksi(ksi,eta,zeta),funDN8Dksi(ksi,eta,zeta);
                        funDN1Deta(ksi,eta,zeta),funDN2Deta(ksi,eta,zeta),funDN3Deta(ksi,eta,zeta),funDN4Deta(ksi,eta,zeta), ...
                        funDN5Deta(ksi,eta,zeta),funDN6Deta(ksi,eta,zeta),funDN7Deta(ksi,eta,zeta),funDN8Deta(ksi,eta,zeta);
                        funDN1Dzeta(ksi,eta,zeta),funDN2Dzeta(ksi,eta,zeta),funDN3Dzeta(ksi,eta,zeta),funDN4Dzeta(ksi,eta,zeta), ...
                        funDN5Dzeta(ksi,eta,zeta),funDN6Dzeta(ksi,eta,zeta),funDN7Dzeta(ksi,eta,zeta),funDN8Dzeta(ksi,eta,zeta)] * ...
                        [pt1x,pt1y,pt1z;pt2x,pt2y,pt2z;pt3x,pt3y,pt3z;pt4x,pt4y,pt4z;pt5x,pt5y,pt5z;pt6x,pt6y,pt6z;pt7x,pt7y,pt7z;pt8x,pt8y,pt8z];
                    
                    Jacobian = det(J);
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
                    
                    % ------- Comment: Calculate DivFMinusW ---------
                    % tempDUDX = DN*FMinusW1(temp);
                    % DivFMinusW1 = tempDUDX(1)+tempDUDX(4);
                    % tempDUDX = DN*FMinusW2(temp);
                    % DivFMinusW2 = tempDUDX(1)+tempDUDX(4);
                    % DivFMinusW = [DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2;
                    %               DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2;
                    %               DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2;
                    %               DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2];
                    % ------ End of comment of calculating DivFMinusW ------
                    
                    % ------ Construct A matrix ------
                    if IterStep == 1
                        % A(tempIndexU,tempIndexU) = A(tempIndexU,tempIndexU) + Jacobian*gqwt(tempi)*gqwt(tempj)*gqwt(tempk) * ( (beta+alpha)*(DN')*(DN) + mu*(N')*N + alpha2*(DNOrig')*DNOrig ) ;
                        tempA = tempA + Jacobian*gqwt(tempi)*gqwt(tempj)*gqwt(tempk) * ( (beta+alpha)*(DN')*(DN) + mu*(N')*N + alpha2*(DNOrig')*DNOrig ) ;
                    end
                    
                    % ------ Construct b vector ------
                    % b(temp) = b(temp) + Jacobian*gqwt(tempi)*gqwt(tempj)* ( -beta*NDiag*(DivFMinusW) + mu*NDiag*(UMinusv(temp)) + alpha*(DN')*DN*U(temp) );
                   
                    % b(tempIndexU) = b(tempIndexU) + Jacobian*gqwt(tempi)*gqwt(tempj)*gqwt(tempk)*( beta*diag((DN')*(FMinusW(tempIndexF)')) + mu*NDiag*(UMinusv(tempIndexU)) +(alpha)*(DN')*DN*U(tempIndexU) );
                    tempb = tempb + Jacobian*gqwt(tempi)*gqwt(tempj)*gqwt(tempk)*( beta*diag((DN')*(FMinusW(tempIndexF)')) + mu*NDiag*(UMinusv(tempIndexU)) +(alpha)*(DN')*DN*U(tempIndexU) );
                end
            end
        end
        
        % Store tempA and tempb
        [IndexAXX,IndexAYY] = ndgrid(tempIndexU,tempIndexU); 
        INDEXAI = [INDEXAI;IndexAXX(:)]; INDEXAJ = [INDEXAJ;IndexAYY(:)]; INDEXAVAL = [INDEXAVAL;tempA(:)];
        INDEXBI = [INDEXBI;tempIndexU(:)]; INDEXBVAL = [INDEXBVAL;tempb(:)];
        
    end
    
    if IterStep == 1
        % A = sparse(DIM*FEMSize+NodesPerEle*DIM, DIM*FEMSize+NodesPerEle*DIM);
        % b = sparse(DIM*FEMSize+NodesPerEle*DIM,1); 
        A = sparse(INDEXAI,INDEXAJ,INDEXAVAL, DIM*FEMSize+NodesPerEle*DIM, DIM*FEMSize+NodesPerEle*DIM);
  
        b = sparse(INDEXBI, ones(length(INDEXBI),1),INDEXBVAL, DIM*FEMSize+NodesPerEle*DIM,1 );
        % AMatrixRegularized = A + 1e-3*max(diag(A))*speye(size(A,1),size(A,2));
        AMatrixRegularized = A;
        % AMatrixRegularized = A; tempVal = 1e-3*max(diag(A));
        % for tempi = 1:size(A,1)
        %     AMatrixRegularized(tempi,tempi) = AMatrixRegularized(tempi,tempi) + tempVal;
        % send
    end
    
    coordsIndexInvolved = unique([0;elementsFEM(:)]); % Need modification for triangle elementsFEM
    %     if CrackOrNot == 1
    %         coordsIndexInvolvedH  = EnrHAndTipEleIndex(:,2);
    %         if CrackTipOrNot == 1
    %             coordsIndexInvolvedC1 = EnrTipEleIndex(:,2); coordsIndexInvolvedC2 = EnrTipEleIndex(:,2)+1;
    %             coordsIndexInvolvedC3 = EnrTipEleIndex(:,2)+2; coordsIndexInvolvedC4 = EnrTipEleIndex(:,2)+3;
    %             coordsIndexInvolved   = [coordsIndexInvolved;coordsIndexInvolvedH;coordsIndexInvolvedC1;coordsIndexInvolvedC2;coordsIndexInvolvedC3;coordsIndexInvolvedC4];
    %         else
    %             coordsIndexInvolved   = [coordsIndexInvolved;coordsIndexInvolvedH];
    %         end
    %     end
     
    % In adaptive mesh, using the following code:
    UIndexInvolved = zeros(DIM*(length(coordsIndexInvolved)-1),1);
    
    % Not including the first 0-th entry
    for tempi = 1:(size(coordsIndexInvolved,1)-1)
        UIndexInvolved(3*tempi-2:3*tempi) = [3*coordsIndexInvolved(tempi+1)-2; ...
                    3*coordsIndexInvolved(tempi+1)-1; 3*coordsIndexInvolved(tempi+1)];
    end
    
    % ========= Set Dirichlet and Neumann boundary conditions =========
    if isempty(dirichlet) ~= 1
        dirichlettemp = [3*dirichlet(:); 3*dirichlet(:)-1; 3*dirichlet(:)-2];
    else
        dirichlettemp = [];
    end
    if isempty(neumann) ~= 1
        neumanntemp = [3*neumann(:,1); 3*neumann(:,1)-1; 3*neumann(:,1)-2; 3*neumann(:,2); 3*neumann(:,2)-1; 3*neumann(:,2)-2];
    else
        neumanntemp = [];
    end
    FreeNodes = setdiff(UIndexInvolved,unique([dirichlettemp]));
    
    % ========= Neumann conditions ===========
    % Last step boundary condition force
    % BCForce = -Global_NodalStrainAvg(coordinatesFEM,elementsFEM,Uhat);
    % for tempj = 1:size(neumann,1)
    %     b(2*neumann(tempj,1:2)-1) = b(2*neumann(tempj,1:2)-1) + norm(coordinatesFEM(neumann(tempj,1),:)-coordinatesFEM(neumann(tempj,2),:)) ...
    %         * ( BCForce(4*neumann(tempj,1:2)-3) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)-1) * neumann(tempj,4) );
    %     b(2*neumann(tempj,1:2))   = b(2*neumann(tempj,1:2)) + norm(coordinatesFEM(neumann(tempj,1),:)-coordinatesFEM(neumann(tempj,2),:)) ...
    %         * ( BCForce(4*neumann(tempj,1:2)-1) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)) * neumann(tempj,4) );
    % end
    
    % ========= Dirichlet conditions ==========
    UhatOld = Uhat;
    UhatNew = sparse(DIM*FEMSize + NodesPerEle*DIM, 1);
    UhatNew(3*unique(dirichlet)) = U(3*unique(dirichlet));
    UhatNew(3*unique(dirichlet)-1) = U(3*unique(dirichlet)-1);
    UhatNew(3*unique(dirichlet)-2) = U(3*unique(dirichlet)-2);
    
    b = b - AMatrixRegularized * UhatNew;
    
    % [tempgridx,tempgridy] = meshgrid(1:FEMSize,1:FEMSize);
    % figure; mesh(tempgridx,tempgridy,AMatrixRegularized(1:FEMSize,1:FEMSize));
    
    % ========= Solve FEM problem ===========
    UhatNew(FreeNodes) = AMatrixRegularized(FreeNodes,FreeNodes) \ b(FreeNodes);
    Uhat = UhatNew(1:end-NodesPerEle*DIM);
     
    
    % if norm(UhatNew-UhatOld)/sqrt(length(UhatOld)) < 1e-3
    %     ConvergeOrNot = 1;
    % end
    
    
end

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





