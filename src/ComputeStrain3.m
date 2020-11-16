% ==================================
% Compute strain
% ----------------------------------
% switch MethodToComputeStrain
%   case 0: direct solved results;
%   case 1: central finite difference;
%   case 2: plane fitting method;
%   case 3: finite element Gauss points;
% ==================================

switch DVCpara.MethodToComputeStrain
    case 0 % Direct output
        Rad = [0,0,0];
        FStrain = FLocal;  % JY!!! change intrisic coords to world coords
        FStraintemp = FStrain;
        
    case 1 % Finite difference method
        Rad = [1,1,1];
        %D=funDerivativeOp3((M-2*Rad(1)),(N-2*Rad(2)),(L-2*Rad(3)),DVCpara.winstepsize); % D = sparse(4*(M-2*Rad)*(N-2*Rad), 2*(M-2*Rad)*(N-2*Rad));
        D2=funDerivativeOp3(MNL(1),MNL(2),MNL(3),DVCpara.winstepsize); 
        [temp3,temp4] = funFDNeumannBCInd3(size(DVCmesh.coordinatesFEM,1),[M,N,L],Rad); % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
        temp = D2*ULocal; FStraintemp = FLocal(temp3);
        FStrain = FLocal; FStrain(temp3) = FStraintemp;
 
    case 2 % Plane fitting method
        D = funDerivativeOp3(M,N,L,DVCpara.winstepsize); % D = sparse(4*M*N, 2*M*N);
        FStrain = D*reshape(ULocal,length(ULocal),1); % JY!!! change intrisic coords to world coords
        
        % Compute strain method II: Use Plane Fitting method
        prompt = 'What is your half window size: ';
        Rad = input(prompt); if length(Rad)==1, Rad=ones(1,3)*Rad; end
        
        [UNew,DUDX,DUDY,DUDZ] = PlaneFit3(reshape(ULocal(1:3:end),M,N,L),DVCpara.winstepsize,Rad,1);
        [VNew,DVDX,DVDY,DVDZ] = PlaneFit3(reshape(ULocal(2:3:end),M,N,L),DVCpara.winstepsize,Rad,2);
        [WNew,DWDX,DWDY,DWDZ] = PlaneFit3(reshape(ULocal(3:3:end),M,N,L),DVCpara.winstepsize,Rad,3);
         
        FStraintemp = zeros(9*(M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        FStraintemp(1:9:end) = reshape(DUDX((Rad(1)+1):M-Rad(1),(Rad(2)+1):N-Rad(2),(Rad(3)+1):L-Rad(3)), (M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        FStraintemp(2:9:end) = reshape(DVDX((Rad(1)+1):M-Rad(1),(Rad(2)+1):N-Rad(2),(Rad(3)+1):L-Rad(3)), (M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        FStraintemp(3:9:end) = reshape(DWDX((Rad(1)+1):M-Rad(1),(Rad(2)+1):N-Rad(2),(Rad(3)+1):L-Rad(3)), (M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        FStraintemp(4:9:end) = reshape(DUDY((Rad(1)+1):M-Rad(1),(Rad(2)+1):N-Rad(2),(Rad(3)+1):L-Rad(3)), (M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        FStraintemp(5:9:end) = reshape(DVDY((Rad(1)+1):M-Rad(1),(Rad(2)+1):N-Rad(2),(Rad(3)+1):L-Rad(3)), (M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        FStraintemp(6:9:end) = reshape(DWDY((Rad(1)+1):M-Rad(1),(Rad(2)+1):N-Rad(2),(Rad(3)+1):L-Rad(3)), (M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        FStraintemp(7:9:end) = reshape(DUDZ((Rad(1)+1):M-Rad(1),(Rad(2)+1):N-Rad(2),(Rad(3)+1):L-Rad(3)), (M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        FStraintemp(8:9:end) = reshape(DVDZ((Rad(1)+1):M-Rad(1),(Rad(2)+1):N-Rad(2),(Rad(3)+1):L-Rad(3)), (M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        FStraintemp(9:9:end) = reshape(DWDZ((Rad(1)+1):M-Rad(1),(Rad(2)+1):N-Rad(2),(Rad(3)+1):L-Rad(3)), (M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
          
        % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
        temp = 1:1:size(coordinatesFEM,1); temp = temp';
        temp = reshape(temp,M,N,L); temp2 = temp(Rad(1)+1:M-Rad(1), Rad(2)+1:N-Rad(2), Rad(3)+1:L-Rad(3));
        temp2 = reshape(temp2, (M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        
        % Plotstrain0(FStraintemp,x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad),f,g);
        temp3 = zeros(9*(M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        for i = 1:(M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3))
            temp3(9*i-8:9*i) = 9*temp2(i)*ones(9,1)+[-8:1:0]';
        end
        FStrain(temp3) = FStraintemp;
        
    case 3 % Finite element method
        
        GaussPtOrder=2; [FStrain,~,~] = funGlobal_NodalStrainAvg3(coordinatesFEM,elementsFEM,ULocal,GaussPtOrder);
     
        Rad = [1,1,1];  
        temp = 1:1:size(coordinatesFEM,1); temp = temp';
        temp = reshape(temp,M,N,L); 
        temp2 = temp(Rad(1)+1:M-Rad(1), Rad(2)+1:N-Rad(2), Rad(3)+1:L-Rad(3));
        temp2 = reshape(temp2, (M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        
        temp3 = zeros(9*(M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)),1);
        for i = 1:(M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3))
            temp3(9*i-8:9*i) = 9*temp2(i)*ones(9,1)+[-8:1:0]';
        end
        FStraintemp = FStrain(temp3);
        
    otherwise
        disp('Wrong Input to compute strain field!')
        
end
 

%% Update infinitesimal strain to large deformation gradient tensor
FStrainFinite = FStrain;
for tempi = 1:9:length(FStrain)
    
    % Obtain each component of def grad tensor
    dudx = FStrain(tempi);  dvdx = FStrain(tempi+1); dwdx = FStrain(tempi+2);
    dudy = FStrain(tempi+3); dvdy = FStrain(tempi+4); dwdy = FStrain(tempi+5);
    dudz = FStrain(tempi+6); dvdz = FStrain(tempi+7); dwdz = FStrain(tempi+8); 
    
    switch DVCpara.StrainType
        case 0 % Infinitesimal stran
            % Do nothing
            
        case 1 % Eluerian strain
            FStrainFinite(tempi) = 1/(1-dudx)-1; FStrainFinite(tempi+1) = dvdx/(1-dudx); FStrainFinite(tempi+2) = dwdx/(1-dudx);
            FStrainFinite(tempi+3) = dudy/(1-dvdy); FStrainFinite(tempi+4) = 1/(1-dvdy)-1; FStrainFinite(tempi+5) = dwdy/(1-dvdy);
            FStrainFinite(tempi+6) = dudz/(1-dwdz); FStrainFinite(tempi+7) = dvdz/(1-dwdz); FStrainFinite(tempi+8) = 1/(1-dwdz)-1; 
             
            % Don't modify this line
            FStraintemp = FStrainFinite(temp3);
            
        case 2 % Green-Lagrangian strain: E=(C-I)/2
            tempF = [1+dudx, dudy, dudz; dvdx, 1+dvdy, dvdz; dwdx, dwdy, 1+dwdz];
            tempC = tempF'*tempF;
            tempE = 0.5*(tempC-eye(3));
            FStrainFinite(tempi) = tempE(1,1);
            FStrainFinite(tempi+1) = tempE(2,1);
            FStrainFinite(tempi+2) = tempE(3,1);
            FStrainFinite(tempi+3) = tempE(1,2);
            FStrainFinite(tempi+4) = tempE(2,2);
            FStrainFinite(tempi+5) = tempE(3,2);
            FStrainFinite(tempi+6) = tempE(1,3);
            FStrainFinite(tempi+7) = tempE(2,3);
            FStrainFinite(tempi+8) = tempE(3,3);
             
            % Don't modify this line
            FStraintemp = FStrainFinite(temp3);
            
        case 3 % Principal strains (for infinitesimal strains)
            tempDU = [dudx, dudy, dudz; dvdx, dvdy, dvdz; dwdx, dwdy, dwdz];
            tempPrincipalStrain = eig(tempDU);
            FStrainFinite(tempi) = tempPrincipalStrain(1);
            FStrainFinite(tempi+1) = 0.5*(tempPrincipalStrain(1)-tempPrincipalStrain(2));
            FStrainFinite(tempi+2) = 0.5*(tempPrincipalStrain(1)-tempPrincipalStrain(3));
            FStrainFinite(tempi+3) = 0.5*(tempPrincipalStrain(1)-tempPrincipalStrain(2));
            FStrainFinite(tempi+4) = tempPrincipalStrain(2);
            FStrainFinite(tempi+5) = 0.5*(tempPrincipalStrain(2)-tempPrincipalStrain(3));
            FStrainFinite(tempi+6) = 0.5*(tempPrincipalStrain(1)-tempPrincipalStrain(3));
            FStrainFinite(tempi+7) = 0.5*(tempPrincipalStrain(2)-tempPrincipalStrain(3));
            FStrainFinite(tempi+8) = tempPrincipalStrain(3);
            if tempi == 1 % Only print following lines once.
                disp('Final plotted principal infinitesimal strains are: ');
                disp('  Subplot(1,1): Principal strain e1');
                disp('  Subplot(1,2): Principal strain e2');
                disp('  Subplot(1,3): Principal strain e3');
                disp('  Subplot(2,1): Max shear (e1-e2)/2');
                disp('  Subplot(2,2): Max shear (e1-e3)/2');
                disp('  Subplot(2,3): Max shear (e2-e3)/2');
                disp('  ');
                disp('Final saved principal infinitesimal strains are: ');
                disp('  ResultStrain.Strain(1:9:end) --> e1');
                disp('  ResultStrain.Strain(2:9:end) --> (e1-e2)/2');
                disp('  ResultStrain.Strain(3:9:end) --> (e1-e3)/2');
                disp('  ResultStrain.Strain(4:9:end) --> (e1-e2)/2');
                disp('  ResultStrain.Strain(5:9:end) --> e2');
                disp('  ResultStrain.Strain(6:9:end) --> (e2-e3)/2');
                disp('  ResultStrain.Strain(7:9:end) --> (e1-e3)/2');
                disp('  ResultStrain.Strain(8:9:end) --> (e2-e3)/2');
                disp('  ResultStrain.Strain(9:9:end) --> e3');
                disp('  ');
            end
        case 4
            disp('Press "Ctrl+C" to modify codes by yourself.'); pause;
            
        otherwise
            disp('Wrong strain type!');
    end

end



    
%% Smoothing strain field
% % Plotstrain(FLocal,x,y,f,g);
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
% temp3 = zeros(9*(M-2*Rad)*(N-2*Rad)*(L-2*Rad),1);
% for i = 1:(M-2*Rad)*(N-2*Rad)*(L-2*Rad)
%     temp3(9*i-8:9*i) = 9*temp2(i)*ones(9,1)+[-8:1:0]';
% end
% FStraintemp = FStrain(temp3);

