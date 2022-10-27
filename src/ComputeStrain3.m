% =========================================================
% Compute strain using solved displacement and deformation gradient results
% =========================================================
% 
% Four methods are provided to compute strain fields
%   Method 1: direct solved results;
%   Method 2: central finite difference;
%   Method 3: plane fitting method;
%   Method 4: finite element Gauss points;
%   
% Strain types:
%   
%
% ---------------------------------------------------------
% Author: Jin Yang, Asst.Prof. @UT-Austin; Postdoc @UW-Madison; PhD '19 @Caltech;
% Contact: aldicdvc@gmail.com; jin.yang@austin.utexas.edu
% Date: 2017-2018.02, 2020.06, 2022.09
% =========================================================

switch DVCpara.strainCalculationMethod

    % ======================================
    case 0 % Direct output
        Rad = [0,0,0];
        FStrain = F_accum;   
        FStrain_crop_no_edges = FStrain;
        
    % ======================================
    case 1 % Finite difference method
        Rad = [1,1,1];
        FDOperator3 = funDerivativeOp3(MNL(1),MNL(2),MNL(3),DVCpara.winstepsize); 
        [notNeumannBCInd_F,notNeumannBCInd_U] = funFDNeumannBCInd3(size(DVCmesh.coordinatesFEM,1),[M,N,L],Rad); % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
        F_accum = FDOperator3*U_accum; FStrain_crop_no_edges = F_accum(notNeumannBCInd_F);
        FStrain = F_accum; FStrain(notNeumannBCInd_F) = FStrain_crop_no_edges;
 
    % ======================================
    case 2 % Plane fitting method
        FDOperator3 = funDerivativeOp3(MNL(1),MNL(2),MNL(3),DVCpara.winstepsize); % D = sparse(4*M*N, 2*M*N);
        FStrain = FDOperator3 * reshape(U_accum,length(U_accum),1); % JY!!! change intrisic coords to world coords
        
        % Compute strain method II: Plane Fitting method
        disp('What is the half length of the edge of the used square fitting plane? ')
        prompt = 'Input here (unit: voxel): ';
        Rad = input(prompt); if length(Rad)==1, Rad=ones(1,3)*Rad; end
        
        [UNew,DUDX,DUDY,DUDZ] = funPlaneFit3(reshape(U_accum(1:3:end),M,N,L),DVCpara.winstepsize,Rad,1);
        [VNew,DVDX,DVDY,DVDZ] = funPlaneFit3(reshape(U_accum(2:3:end),M,N,L),DVCpara.winstepsize,Rad,2);
        [WNew,DWDX,DWDY,DWDZ] = funPlaneFit3(reshape(U_accum(3:3:end),M,N,L),DVCpara.winstepsize,Rad,3);
         
        FStrain_crop_no_edges = zeros(9*(M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        xyz_ind_crop = [(Rad(1)+1):M-Rad(1), (Rad(2)+1):N-Rad(2), (Rad(3)+1):L-Rad(3)];
        xyz_pixels_crop = (M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3));

        FStrain_crop_no_edges(1:9:end) = reshape(DUDX(xyz_ind_crop), xyz_pixels_crop, 1);
        FStrain_crop_no_edges(2:9:end) = reshape(DVDX(xyz_ind_crop), xyz_pixels_crop, 1);
        FStrain_crop_no_edges(3:9:end) = reshape(DWDX(xyz_ind_crop), xyz_pixels_crop, 1);
        FStrain_crop_no_edges(4:9:end) = reshape(DUDY(xyz_ind_crop), xyz_pixels_crop, 1);
        FStrain_crop_no_edges(5:9:end) = reshape(DVDY(xyz_ind_crop), xyz_pixels_crop, 1);
        FStrain_crop_no_edges(6:9:end) = reshape(DWDY(xyz_ind_crop), xyz_pixels_crop, 1);
        FStrain_crop_no_edges(7:9:end) = reshape(DUDZ(xyz_ind_crop), xyz_pixels_crop, 1);
        FStrain_crop_no_edges(8:9:end) = reshape(DVDZ(xyz_ind_crop), xyz_pixels_crop, 1);
        FStrain_crop_no_edges(9:9:end) = reshape(DWDZ(xyz_ind_crop), xyz_pixels_crop, 1);
          
        % Find coordinatesFEM that belong to "xyz_ind_crop"
        nodeInd = 1:1:size(coordinatesFEM,1); nodeInd = nodeInd';
        nodeInd = reshape(nodeInd,M,N,L); nodeInd_crop_no_edges = nodeInd(xyz_ind_crop);
        nodeInd_crop_no_edges = reshape(nodeInd_crop_no_edges, xyz_pixels_crop, 1);
        
        % Plotstrain0(FStraintemp,x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad),f,g);
        notNeumannBCInd_F = zeros(9*(M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        for i = 1:(M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3))
            notNeumannBCInd_F(9*i-8:9*i) = 9*nodeInd_crop_no_edges(i)*ones(9,1)+[-8:1:0]';
        end
        FStrain(notNeumannBCInd_F) = FStrain_crop_no_edges;
       
    % ======================================
    case 3 % Finite element method
        
        GaussPtOrder=2; [FStrain,~,~] = funGlobal_NodalStrainAvg3(coordinatesFEM,elementsFEM,U_accum,GaussPtOrder);
     
        Rad = [1,1,1]; %Try to remove results near FE-mesh edges 
        nodeInd = 1:1:size(coordinatesFEM,1); nodeInd = nodeInd';
        nodeInd = reshape(nodeInd,M,N,L); 
        nodeInd_crop_no_edges = nodeInd(Rad(1)+1:M-Rad(1), Rad(2)+1:N-Rad(2), Rad(3)+1:L-Rad(3));
        nodeInd_crop_no_edges = reshape(nodeInd_crop_no_edges, (M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        
        notNeumannBCInd_F = zeros(9*(M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        for i = 1:(M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3))
            notNeumannBCInd_F(9*i-8:9*i) = 9*nodeInd_crop_no_edges(i)*ones(9,1)+[-8:1:0]';
        end
        FStrain_crop_no_edges = FStrain(notNeumannBCInd_F);
     
    % ======================================
    otherwise
        
        disp('Wrong Input to compute strain field!')
        
end
 

%% Update infinitesimal strain to large deformation gradient tensor
FStrainFinite = FStrain;
for tempi = 1:9:length(FStrain)
    
    % Obtain each component of def grad tensor
    dudx = FStrain(tempi);   dvdx = FStrain(tempi+1); dwdx = FStrain(tempi+2);
    dudy = FStrain(tempi+3); dvdy = FStrain(tempi+4); dwdy = FStrain(tempi+5);
    dudz = FStrain(tempi+6); dvdz = FStrain(tempi+7); dwdz = FStrain(tempi+8); 
    
    switch DVCpara.strainType
        % ======================================
        case 0 % Infinitesimal stran
            % Do nothing
            
        % ======================================
        case 1 % Eluerian-Almansi finite strain: e = (I-inv(B))/2; B = F*F';
            tempFMatrix = [1+dudx, dudy, dudz; dvdx, 1+dvdy, dvdz; dwdx, dwdy, 1+dwdz];
            tempBMatrix = tempFMatrix * tempFMatrix';
            tempeMatrix = 0.5*(eye(3)-inv(tempBMatrix));

            FStrainFinite(tempi) = tempeMatrix(1,1);
            FStrainFinite(tempi+1) = tempeMatrix(2,1);
            FStrainFinite(tempi+2) = tempeMatrix(3,1);
            FStrainFinite(tempi+3) = tempeMatrix(1,2);
            FStrainFinite(tempi+4) = tempeMatrix(2,2);
            FStrainFinite(tempi+5) = tempeMatrix(3,2);
            FStrainFinite(tempi+6) = tempeMatrix(1,3);
            FStrainFinite(tempi+7) = tempeMatrix(2,3);
            FStrainFinite(tempi+8) = tempeMatrix(3,3);

            % Don't modify this line
            FStrain_crop_no_edges = FStrainFinite(notNeumannBCInd_F);
            
        % ======================================
        case 2 % Green-Lagrangian finite strain: E=(C-I)/2; C = F'*F;
            tempFMatrix = [1+dudx, dudy, dudz; dvdx, 1+dvdy, dvdz; dwdx, dwdy, 1+dwdz];
            tempCMatrix = tempFMatrix' * tempFMatrix;
            tempEMatrix = 0.5*(tempCMatrix-eye(3));
            
            FStrainFinite(tempi) = tempEMatrix(1,1);
            FStrainFinite(tempi+1) = tempEMatrix(2,1);
            FStrainFinite(tempi+2) = tempEMatrix(3,1);
            FStrainFinite(tempi+3) = tempEMatrix(1,2);
            FStrainFinite(tempi+4) = tempEMatrix(2,2);
            FStrainFinite(tempi+5) = tempEMatrix(3,2);
            FStrainFinite(tempi+6) = tempEMatrix(1,3);
            FStrainFinite(tempi+7) = tempEMatrix(2,3);
            FStrainFinite(tempi+8) = tempEMatrix(3,3);
             
            % Don't modify this line
            FStrain_crop_no_edges = FStrainFinite(notNeumannBCInd_F);
        
        % ======================================
        case 3 % Hencky strain: h = 1/2*ln(B) = sum_{i=1}^{3} (ln lambda_i) n_i \otimes n_i;  B = F*F';
            tempFMatrix = [1+dudx, dudy, dudz; dvdx, 1+dvdy, dvdz; dwdx, dwdy, 1+dwdz];
            tempBMatrix = tempFMatrix * tempFMatrix';
            [tempVMatrix,tempDMatrix] = eig(tempBMatrix);
            
            tempHMatrix = log(tempDMatrix(1,1)) * tempVMatrix(:,1)*tempVMatrix(:,1)' + ...
                          log(tempDMatrix(2,2)) * tempVMatrix(:,2)*tempVMatrix(:,2)' + ...
                          log(tempDMatrix(3,3)) * tempVMatrix(:,3)*tempVMatrix(:,3)';
        
            FStrainFinite(tempi) = tempHMatrix(1,1);
            FStrainFinite(tempi+1) = tempHMatrix(2,1);
            FStrainFinite(tempi+2) = tempHMatrix(3,1);
            FStrainFinite(tempi+3) = tempHMatrix(1,2);
            FStrainFinite(tempi+4) = tempHMatrix(2,2);
            FStrainFinite(tempi+5) = tempHMatrix(3,2);
            FStrainFinite(tempi+6) = tempHMatrix(1,3);
            FStrainFinite(tempi+7) = tempHMatrix(2,3);
            FStrainFinite(tempi+8) = tempHMatrix(3,3);
             
            % Don't modify this line
            FStrain_crop_no_edges = FStrainFinite(notNeumannBCInd_F);

        % ======================================
        case 4
            disp('Press "Ctrl + C" to modify codes by yourself.'); pause;
          
        % ======================================
        otherwise
            disp('Unknown strain type!');
    end

end



    
%% Smoothing strain field
% % Plotstrain(F_accum,x,y,f,g);
% % for tempi = 1:3, FStrain = funSmoothStrain3(FStrain,coordinatesFEM,elementsFEM,winstepsize,0,0); end
%     
% % End of computing strain method I.
% 
% Rad = 2; % FilterSizeInput-2;
% % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
% temp = 1:1:size(coordinatesFEM,1); temp = temp';
% temp = reshape(temp,M,N,L); temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad, Rad+1:L-Rad);
% temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad)*(L-2*Rad),1);
% 
% notNeumannBCInd_F = zeros(9*(M-2*Rad)*(N-2*Rad)*(L-2*Rad),1);
% for i = 1:(M-2*Rad)*(N-2*Rad)*(L-2*Rad)
%     notNeumannBCInd_F(9*i-8:9*i) = 9*temp2(i)*ones(9,1)+[-8:1:0]';
% end
% FStraintemp = FStrain(notNeumannBCInd_F);


