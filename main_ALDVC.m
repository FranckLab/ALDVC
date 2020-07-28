% ---------------------------------------------------
% Augmented Lagrangian Digital Volume Correlation (ALDVC)
% Author: Jin Yang, Postdoc UW-Madison; PhD '19 Caltech;
% Contact: jyang526@wisc.edu;  aldicdvc@gmail.com
% 2017-2018.02, 2020.06
% ---------------------------------------------------

%% Section 0
% ====== Prepare volumetric image data files ======
% Please go to subfolder code: './DVC_images/GenerateVolMatfile.m' to
% transform your volumetric image stacks to a Matlab matfile for ALDVC code.

%% Section 1
% ====== Clear MATLAB environment & mex set up tricubic interpolation ======
close all; clear all; clc; clearvars -global
fprintf('------------ Section 1 Start ------------ \n')
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')
mex -O ba_interp3.cpp; warning('off'); % dbstop if error 
addpath('./func','./src','./plotFiles','./DVC_images','./plotFiles/export_fig-d966721','./func/regularizeNd');
fprintf('------------ Section 1 Done ------------ \n \n')


%% Section 2 
fprintf('------------ Section 2 Start ------------ \n')
% ====== Read images ======
filename = 'vol_str*.mat'; fileFolder = './DVC_images/';
try if isempty(fileFolder)~=1, cd(fileFolder); end; catch; end
[file_name,Img,DVCpara] = ReadImage3(filename); 
DVCpara.interpmethod='cubic';       % Grayscale interpolation scheme: choose from {'linear','cubic','spline','default'}
DVCpara.displayIterOrNot=0;         % Display Section 4 Subpb1 IC-GN convergence info
DVCpara.Subpb1ICGNMaxIterNum=100;   % Maximum IC-GN iterations in Sectino 4 Subpb1
DVCpara.ICGNtol=1e-2;               % Subproblem 1 IC-GN stopping threshold
try if isempty(fileFolder)~=1, cd('../'); end; catch; end
% ====== Uncomment the behind line and change the value you want ======
% E.g.: gridRange.gridxRange = [352,696]; ... 
% ====== Normalize images ======
[ImgNormalized,DVCpara.gridRange] = funNormalizeImg3(Img,DVCpara.gridRange,'Normalize'); clear Img; % Clear original images to save space
% ====== Initialize variable storage ======
ImgSeqLength=length(ImgNormalized);
ResultDisp = cell(ImgSeqLength-1,1);         ResultDefGrad = cell(ImgSeqLength-1,1);
ResultStrain = cell(ImgSeqLength-1,1);       Resultmubeta = cell(ImgSeqLength-1,1);
ResultConvItPerEle = cell(ImgSeqLength-1,1); ResultcoordinatesFEM = cell(ImgSeqLength-1,1);
ResultelementsFEM = cell(ImgSeqLength-1,1);  ResultSizeOfFFTSearchReg = cell(ImgSeqLength-1,1);
ResultFEMesh = cell(ceil((ImgSeqLength-1)/DVCpara.ImgSeqIncUnit),1); % For incremental DIC mode
fprintf('------------ Section 2 Done ------------ \n \n')
 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start each frame in an image sequence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Img{1} = ImgNormalized{1}; % The first frame is reference 
for ImgSeqNum = 2:length(ImgNormalized)
     
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    
    %% Section 3: Find initial guess
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to find or update initial guess for ALDVC
    % The key idea is to either to use FFT peak fitting, or use last frame
    % results for the next new frame;
    % Particularly in incremental mode, the reference image can also be updated.
    % fNormalized = ImgNormalized{ImgSeqNum-mod(ImgSeqNum-1,ImgSeqIncUnit)};
    Img{2} = ImgNormalized{ImgSeqNum}; NewFFTSearchCheck = 0; DVCpara.NewFFTSearch = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ImgSeqNum==2 || DVCpara.NewFFTSearch==1
        % ====== Integer Search ======
        % DVCpara.InitFFTMethod = 'bigxcorr'; tic % fftmethod \in {'xcorr','phasecorr','bigphasecorr','bigxcorr'}
        tic; [xyz0,uvw0,cc,SizeOfFFTSearchRegion] = IntegerSearch3Mg(Img,DVCpara); DVCpara.SizeOfFFTSearchRegion=SizeOfFFTSearchRegion; toc
        % ======== Find some bad inital guess points ========
        cc.ccThreshold=1.25; % bad cross-correlation threshold (mean - ccThreshold*stdev for q-factor distribution)
        DVCpara.qDICOrNot=0; DVCpara.Thr0=0; [uvw,cc]=RemoveOutliers3(uvw0,cc,DVCpara.qDICOrNot,DVCpara.Thr0); %Last term is threshold value
        % ====== FEM mesh set up ======
        [DVCmesh] = MeshSetUp3(xyz0,DVCpara);
        % ====== Assign initial values ======
        U0 = Init3(uvw,DVCmesh.xyz0); Plotdisp_show3(U0,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM); 
    
        % ====== Deal with incremental mode ======
        Img1NewIndex = ImgSeqNum-mod(ImgSeqNum-2,DVCpara.ImgSeqIncUnit)-1;
        if DVCpara.ImgSeqIncUnit==1, Img1NewIndex = Img1NewIndex-1; end
        ResultFEMesh{1+floor(Img1NewIndex/DVCpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM',DVCmesh.coordinatesFEM,'elementsFEM',DVCmesh.elementsFEM, ...
            'winsize',DVCpara.winsize,'winstepsize',DVCpara.winstepsize,'gridxyROIRange',DVCpara.gridRange );
       
    elseif mod(ImgSeqNum-2,DVCpara.ImgSeqIncUnit)==0 % TO update ref image in incremental mode
        Img1NewIndex = ImgSeqNum-mod(ImgSeqNum-2,DVCpara.ImgSeqIncUnit)-1;
        if DVCpara.ImgSeqIncUnit==1,  Img1NewIndex = Img1NewIndex-1; end
        Img{1} = ImgNormalized{Img1NewIndex}; % Update reference
        [DVCpara,DVCmesh] = ReadImageRefUpdate(file_name,ImgSeqNum,ResultDisp{ImgSeqNum-2}.U,DVCpara,DVCmesh); % Update reference image if needed;
        U0 = zeros(3*size(DVCmesh.coordinatesFEM,1),1); % PlotuvInit;
        ResultFEMesh{1+floor(Img1NewIndex/DVCpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM',DVCmesh.coordinatesFEM,'elementsFEM',DVCmesh.elementsFEM, ...
            'winsize',DVCpara.winsize,'winstepsize',DVCpara.winstepsize,'gridxyROIRange',DVCpara.gridRange );
    else
        U0 = ResultDisp{ImgSeqNum-2}.U;
    end
    % ====== Spline interpolation images ======
    % Df = funImgGradient3(Img{1},'stencil7'); Df.imgSize = size(Img{1}); % Compute image grayscale value gradients
    Df = struct(); Df.imgSize = size(Img{1});
    % ====== Compute f(X)-g(x+u) ======
    % PlotImgDiff(x0,y0,u,v,fNormalized,gNormalized); % Img grayscale value residual
    fprintf('------------ Section 3 Done ------------ \n \n')
 
     
    %% Section 4
    fprintf('------------ Section 4 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve first local step in ALDVC: Subproblem 1
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ====== ALStep 1 Subproblem1: Local Subset DVC ======
    mu=0; beta=0; ALSolveStep=1; ALSub1Time=zeros(6,1); ALSub2Time=zeros(6,1); ConvItPerEle=zeros(size(DVCmesh.coordinatesFEM,1),6);  
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DVCpara.interpmethod='cubic';  % Choose from {'linear','cubic','spline','default'}
    % DVCpara.displayIterOrNot=0; DVCpara.Subpb1ICGNMaxIterNum=100; 
    % ------ Start Local DIC IC-GN iteration ------
    [USubpb1,FSubpb1,HtempPar,ALSub1Timetemp,ConvItPerEletemp] = LocalICGN3( ...
        U0,DVCmesh.coordinatesFEM,Df,Img{1},Img{2},DVCpara,'GaussNewton',DVCpara.ICGNtol);
    ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; toc
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------  Manually find some bad points from Local Subset ICGN step ------
    disp('--- Start to remove bad local points ---')
    MNL = size(DVCmesh.xyz0.x); uvw.u = reshape(USubpb1(1:3:end),MNL); uvw.v = reshape(USubpb1(2:3:end),MNL); uvw.w = reshape(USubpb1(3:3:end),MNL); 
    [uvw,cc,RemoveOutliersList] = RemoveOutliers3(uvw,[],DVCpara.qDICOrNot,DVCpara.Thr0);
    USubpb1 = [uvw.u(:),uvw.v(:),uvw.w(:)]'; USubpb1 = USubpb1(:); FSubpb1 = FSubpb1(:);
    for tempi = 0:8
        FSubpb1(9*RemoveOutliersList-tempi) = nan;
        FSubpb1(9-tempi:9:end) = reshape((inpaint_nans3(reshape(FSubpb1(9-tempi:9:end),MNL),1)),size(DVCmesh.coordinatesFEM,1),1);
    end
    disp('--- Remove bad points done ---')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Plot & save variables ------
    close all; Plotdisp_show3(USubpb1,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM); % Plotuvw(USubpb1,[],xyz0,[]);
    Plotstrain_show3(FSubpb1,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM);
    save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
    ULocalICGN = USubpb1; FLocalICGN = FSubpb1; % These are results from local subset based method 
    fprintf('------------ Section 4 Done ------------ \n \n')
     
    
    %% Section 5
    fprintf('------------ Section 5 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve global step in ALDIC: Subproblem 2
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ======= ALStep 1 Subproblem 2: Global constraint =======
    % ------ Smooth displacements for a little bit better F ------
    DVCpara.DispFilterSize=0; DVCpara.DispFilterStd=0; DVCpara.StrainFilterSize=0; DVCpara.StrainFilterStd=0; LevelNo=1;
    for tempi=1:3, FSubpb1 = funSmoothStrain3(FSubpb1,DVCmesh,DVCpara); end
    % ------ Smooth displacements ------
    % prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    % DoYouWantToSmoothOnceMore = input(prompt); DispFilterSize=0; DispFilterStd=0;
    % if DoYouWantToSmoothOnceMore == 0,USubpb1 = funSmoothDisp(USubpb1,coordinatesFEM,elementsFEM,winstepsize,DispFilterSize,DispFilterStd);end
    % for tempi = 1:SmoothTimesDisp, USubpb1 = funSmoothDisp3(USubpb1,DVCmesh,DVCpara); end
    %Plotdisp_show3(USubpb1,coordinatesFEM,elementsFEM);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ====== Define penalty parameter ======
    mu=1e-3; udual=0*FSubpb1; vdual=0*USubpb1; alpha=0; % Ignore here alpha 
    betaList = [sqrt(1e-5),1e-2,sqrt(1e-3),1e-1,sqrt(1e-1)]*mean(DVCpara.winstepsize).^2.*mu; 
    Err1=zeros(length(betaList),1); Err2=Err1; 
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ====== Check to use FD or FE methods to solve Subpb2 step ======
    if strcmp(DVCpara.Subpb2FDOrFEM,'FD') == 1 % Using finite difference method
        % ====== Build sparse finite difference operator ======
        MNL = 1 + (DVCmesh.coordinatesFEM(end,:)-DVCmesh.coordinatesFEM(1,:))./DVCpara.winstepsize;
        tic; Rad=[1,1,1]; %D=funDerivativeOp3((MNL(1)-2*Rad(1)),(MNL(2)-2*Rad(2)),(MNL(3)-2*Rad(3)),DVCpara.winstepsize); % D = sparse(4*(M-2*Rad)*(N-2*Rad), 2*(M-2*Rad)*(N-2*Rad));
        D2 = funDerivativeOp3(MNL(1),MNL(2),MNL(3),DVCpara.winstepsize); toc
        disp('Finish assembling finite difference operator D');
        
        % ===== Solver using finite difference approximation ======
        tic; a = FSubpb1-udual; b = USubpb1-vdual; a=a(:); b=b(:);
        Rad=[1,1,1]; [temp3,temp4] = funFDNeumannBCInd3(size(DVCmesh.coordinatesFEM,1),MNL,Rad); % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
        atemp = a(temp3); btemp = b(temp4); hbar = waitbar(0,'Please wait for Subproblem 2 global step!'); 
        for tempk = 1:length(betaList)
            beta = betaList(tempk);
            tempAMatrixSub2 = (beta*(D2')*D2) + mu*speye(DVCpara.DIM*(MNL(1))*(MNL(2))*(MNL(3)));
            USubpb2temp = (tempAMatrixSub2) \ (beta*D2'*a + mu*b) ;
            USubpb2 = USubpb1; USubpb2(temp4) = USubpb2temp(temp4);
            FSubpb2 = FSubpb1; temp = D2*USubpb2temp; FSubpb2(temp3) = temp(temp3);
            % FSubpb2 = D2*USubpb2; % Recall: D2 = funDerivativeOp(M,N,winstepsize);
            
            Err1(tempk) = norm(USubpb1(:)-USubpb2(:),2);
            Err2(tempk) = norm(FSubpb1(:)-FSubpb2(:),2);
            waitbar(tempk/(length(betaList)+1));
        end
        ErrSum = Err1+Err2*mean(DVCpara.winstepsize)^2; [~,indexOfbeta] = min(ErrSum);
        try
            [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
            p = coeffvalues(fitobj); beta = 10^(-p(2)/2/p(1));
        catch
            beta = betaList(indexOfbeta);
        end
        % Using optimal beta to solve again
        tempAMatrixSub2 = (beta*(D2')*D2) + mu*speye(DVCpara.DIM*(MNL(1))*(MNL(2))*(MNL(3)));
        USubpb2temp = (tempAMatrixSub2) \ (beta*D2'*a + mu*b) ;
        USubpb2 = USubpb1; USubpb2(temp4) = USubpb2temp(temp4);
        waitbar(1); close(hbar);
        %%%%%%%%%%%%%% End of using finite difference approximation %%%%%%%%%%%%%%
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
	else %Subpb2FDOrFEM: Using finite element method
        MNL = size(DVCmesh.xyz0); GaussPtOrder=2; alpha=0; hbar = waitbar(0,'Please wait for Subproblem 2 global step!'); 
        % ====== Solver using finite element method ======
        for tempk = 1:length(betaList)
            beta = betaList(tempk);
            [USubpb2] = Subpb23(DVCmesh,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,GaussPtOrder);
            [FSubpb2,~,~] = funGlobal_NodalStrainAvg3(DVCmesh.coordinatesFEM,DVCmesh.elementsFEM,USubpb2,GaussPtOrder);
            Err1(tempk) = norm(USubpb1-USubpb2,2);
            Err2(tempk) = norm(FSubpb1-FSubpb2,2);
            waitbar(tempk/(length(betaList)+1));
        end
        Err1Norm = (Err1-mean(Err1))/std(Err1); figure, plot(Err1Norm);
        Err2Norm = (Err2-mean(Err2))/std(Err2); figure, plot(Err2Norm);
        ErrSum = Err1Norm+Err2Norm; figure,plot(ErrSum); [~,indexOfbeta] = min(ErrSum);
        try
            [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
            p = coeffvalues(fitobj); beta = 10^(-p(2)/2/p(1));
        catch
            beta = betaList(indexOfbeta);
        end
        % Using optimal beta to solve again
        [USubpb2] = Subpb23(DVCmesh,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,GaussPtOrder);
        USubpb2 = full(USubpb2); waitbar(1); close(hbar);
    end
	ALSub2Time(ALSolveStep) = toc; toc
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------- Before computing strain, we smooth the displacement field a little bit -------
    % USubpb2 = funSmoothDisp3(USubpb2,coordinatesFEM,elementsFEM,x0,y0,winstepsize,DispFilterSize,DispFilterStd);
    % ------- Compute strain field -------- 
    if strcmp(DVCpara.Subpb2FDOrFEM,'FD') == 1 %Subpb2: Using finite different method
        FSubpb2 = FSubpb1; temp = D2*USubpb2temp; FSubpb2(temp3) = temp(temp3);
    else %Subpb2: Using finite element method
        [FSubpb2,~,~] = funGlobal_NodalStrainAvg3(DVCmesh.coordinatesFEM,DVCmesh.elementsFEM,USubpb2,GaussPtOrder);
    end
    
    % ------- Smooth strain field --------
    %for tempi = 1:3, FSubpb2 = funSmoothStrain3(FSubpb2,DVCmesh,DVCpara); end
     
    % ------- Save data ------
    save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');

    % ------ Plot ------
    Plotdisp_show3(USubpb2,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM);
    Plotstrain_show3(FSubpb2,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM);
    % Or try these codes for a tight view:
    % Plotdisp03(USubpb2,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM,'All');
    % Plotstrain03(full(FSubpb2),xyz0.x,xyz0.y,xyz0.z,size(Img{1}),'All');
    
    % ======= Update dual variables =======
    if strcmp(DVCpara.Subpb2FDOrFEM,'FD') == 1 %Subpb2: Using finite different method
        udualtemp1 = (FSubpb2 - FSubpb1); udualtemp2 = udualtemp1(temp3);
        vdualtemp1 = (USubpb2 - USubpb1); vdualtemp2 = vdualtemp1(temp4);
        udual = zeros(DVCpara.DIM^2*prod(MNL),1); vdual = zeros(DVCpara.DIM*prod(MNL),1);
        udual(temp3) = udualtemp2; vdual(temp4) = vdualtemp2;
    else % FEM or other methods
        udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1;
        for tempi=0:8, udual(9*DVCmesh.dirichlet-tempi) = 0; end
        for tempi=0:2, vdual(3*DVCmesh.dirichlet-tempi) = 0; end
    end
    save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
    fprintf('------------ Section 5 Done ------------ \n \n')


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Section 6
    fprintf('------------ Section 6 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to run ADMM iteration: Subproblem 1 & 2
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ==================== ADMM AL Loop ==========================
    ALSolveStep=1; UpdateY=1e4; tol2=1e-4; % Ignore 1e4, this value>1 will not be used
    HPar=cell(size(HtempPar,2),1); for tempj=1:size(HtempPar,2), HPar{tempj} = HtempPar(:,tempj); end
    while (ALSolveStep < 4)
        ALSolveStep = ALSolveStep + 1;  % Update using the last step
        
        %%%%%%%%%%%%%%%%%%%%%%% Subproblem 1 %%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****']);
        tic;[USubpb1,~,ALSub1Timetemp,ConvItPerEletemp] = Subpb13(USubpb2,FSubpb2,udual,vdual,DVCmesh.coordinatesFEM,...
                Df,Img{1},Img{2},mu,beta,HPar,ALSolveStep,DVCpara,'GaussNewton',DVCpara.ICGNtol); toc
        FSubpb1 = FSubpb2; ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------  Manually find some bad points from Local Subset ICGN step ------
        % Comment START
        % disp('--- Start to manually remove bad points ---')
        % close all; Plotdisp_show3(USubpb1,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM); % Plotuvw(USubpb1,[],xyz0,[]);
        % [row,~] = find(ConvItPerEle(:,ALSolveStep)>DVCpara.Subpb1ICGNMaxIterNum);
        % USubpb1(3*row-2) = NaN; USubpb1(3*row-1) = NaN; USubpb1(3*row) = NaN;
        % MNL = size(DVCmesh.xyz0.x); uvw.u = reshape(USubpb1(1:3:end),MNL); uvw.v = reshape(USubpb1(2:3:end),MNL); uvw.w = reshape(USubpb1(3:3:end),MNL);
        % [uvw,cc,RemoveOutliersList] = RemoveOutliers3(uvw,[],DVCpara.qDICOrNot,DVCpara.Thr0);
        % USubpb1 = [uvw.u(:),uvw.v(:),uvw.w(:)]'; USubpb1 = USubpb1(:); FSubpb1 = FSubpb1(:);
        % for tempi = 0:8,
        %     FSubpb1(9*RemoveOutliersList-tempi) = nan;
        %     FSubpb1(9-tempi:9:end) = reshape((inpaint_nans3(reshape(FSubpb1(9-tempi:9:end),MNL),1)),size(DVCmesh.coordinatesFEM,1),1);
        % end
        % Comment END
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
        %for tempi = 1:3, USubpb1 = funSmoothDisp3(USubpb1,DVCmesh,DVCpara); end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ============== Subproblem 2 ==============
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****'])
        if strcmp(DVCpara.Subpb2FDOrFEM,'FD') == 1 %Subpb2: Using finite different method
            tic; a = FSubpb1-udual; b = USubpb1-vdual;
            USubpb2temp = (tempAMatrixSub2) \ (beta*D2'*a + mu*b) ;
            USubpb2 = USubpb1; USubpb2(temp4) = USubpb2temp(temp4);
        else % FEM or other methods
            tic; [USubpb2] = Subpb23(DVCmesh,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,GaussPtOrder); % [] means I don't apply dirichlet & neumann BC here.
            USubpb2 = full(USubpb2);
        end
        ALSub2Time(ALSolveStep) = toc; toc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for tempi = 0:2, USubpb2(3*dirichlet-tempi) = USubpb1(3*dirichlet-tempi); end
        % ------- Before computing strain, we smooth the displacement field -------
        % USubpb2 = funSmoothDisp3(USubpb2,DVCmesh,DVCpara);
        % ------- Compute strain field --------
        if strcmp(DVCpara.Subpb2FDOrFEM,'FD')==1 %Subpb2: Using finite different method
            FSubpb2 = FSubpb1; temp = D2*USubpb2temp; FSubpb2(temp3) = temp(temp3); 
        else %FEM: Using finite element method
            GaussPtOrder=2; [FSubpb2,~,~] = funGlobal_NodalStrainAvg3(DVCmesh.coordinatesFEM,DVCmesh.elementsFEM,USubpb2,GaussPtOrder);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------- Smooth strain field --------
        % for tempi = 1:3, FSubpb2 = funSmoothStrain3(FSubpb2,DVCmesh,DVCpara); end
        save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute norm of UpdateY
        USubpb2_Old = load(['Subpb2_step',num2str(ALSolveStep-1)],'USubpb2');
        USubpb2_New = load(['Subpb2_step',num2str(ALSolveStep)],'USubpb2');
        USubpb1_Old = load(['Subpb1_step',num2str(ALSolveStep-1)],'USubpb1');
        USubpb1_New = load(['Subpb1_step',num2str(ALSolveStep)],'USubpb1');
        if (mod(ImgSeqNum-2,DVCpara.ImgSeqIncUnit) ~= 0 && (ImgSeqNum>2)) || (ImgSeqNum < DVCpara.ImgSeqIncUnit)
            UpdateY = norm((USubpb2_Old.USubpb2 - USubpb2_New.USubpb2), 2)/sqrt(size(USubpb2_Old.USubpb2,1));
            try
                UpdateY2 = norm((USubpb1_Old.USubpb1 - USubpb1_New.USubpb1), 2)/sqrt(size(USubpb1_Old.USubpb1,1));
            catch
            end
        end
        try
            disp(['Update local step  = ',num2str(UpdateY2)]);
            disp(['Update global step = ',num2str(UpdateY)]);
        catch
        end
        fprintf('*********************************** \n \n');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update dual variables------------------------------
        if strcmp(DVCpara.Subpb2FDOrFEM,'FD') == 1 %%%%%%%%% Using finite difference method %%%%%%%
            udualtemp1 =  (FSubpb2(:) - FSubpb1(:)); udualtemp2 = udualtemp1(temp3);
            vdualtemp1 =  (USubpb2(:) - USubpb1(:)); vdualtemp2 = vdualtemp1(temp4);
            udual(temp3) = udual(temp3)+udualtemp2;
            vdual(temp4) = vdual(temp4)+vdualtemp2;
        else  %%%%%%%%%%%% Using finite element method %%%%%%%%%%% 
            udual = FSubpb2(:) - FSubpb1(:); vdual = USubpb2(:) - USubpb1(:);
            for tempi=0:8, udual(9*DVCmesh.dirichlet-tempi) = 0; end
            for tempi=0:2, vdual(3*DVCmesh.dirichlet-tempi) = 0; end
        end
        
        save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        try if UpdateY<tol2 || UpdateY2<tol2, break; end; catch; end 

    end
    fprintf('------------ Section 6 Done ------------ \n \n')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot and save data
    % ------ Plot ------
    close all; Plotdisp_show3(USubpb2,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM);
    Plotstrain_show3(FSubpb1,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM);
    % Or try these codes for a tight view:
    % Plotdisp03(USubpb2,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM,'All');
    % Plotstrain03(full(FSubpb2),xyz0.x,xyz0.y,xyz0.z,size(Img{1}),'All');
    
    ResultDisp{ImgSeqNum-1}.U = full(USubpb2);
    ResultDisp{ImgSeqNum-1}.ULocalICGN = full(ULocalICGN);
    ResultDisp{ImgSeqNum-1}.U0 = full(U0);
    
    ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb2);
    ResultDefGrad{ImgSeqNum-1}.FLocalICGN = full(FLocalICGN);
    
    Resultmubeta{ImgSeqNum-1}.beta = beta; Resultmubeta{ImgSeqNum-1}.mu = mu;
    ResultConvItPerEle{ImgSeqNum-1}.ConvItPerEle = ConvItPerEle;
    ResultcoordinatesFEM{ImgSeqNum-1}.coordinatesFEM = DVCmesh.coordinatesFEM;
    ResultelementsFEM{ImgSeqNum-1}.elementsFEM = DVCmesh.elementsFEM;
    ResultSizeOfFFTSearchReg{ImgSeqNum-1}.SizeOfFFTSearchReg = SizeOfFFTSearchRegion;
    
end
 

%% Section 7
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to check convergence of ADMM
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('------------ Section 7 Start ------------ \n')
% ====== Check convergence ======
ALSolveStep1 = min(6,ALSolveStep);
disp('====== |u^-u| ======');
for ALSolveStep = 1:ALSolveStep1
    USubpb2 = load(['Subpb2_step',num2str(ALSolveStep )],'USubpb2'); USubpb2=USubpb2.USubpb2(:);
    USubpb1 = load(['Subpb1_step',num2str(ALSolveStep )],'USubpb1'); USubpb1=USubpb1.USubpb1(:);
    UpdateY = norm((USubpb2-USubpb1), 2)/sqrt(length(USubpb2));
    disp(num2str(UpdateY));
end
disp('====== |grad(u^)-F| ======');
for ALSolveStep = 1:ALSolveStep1
    FSubpb1 = load(['Subpb1_step',num2str(ALSolveStep )],'FSubpb1'); FSubpb1=FSubpb1.FSubpb1(:);
    FSubpb2 = load(['Subpb2_step',num2str(ALSolveStep )],'FSubpb2'); FSubpb2=FSubpb2.FSubpb2(:);
    UpdateF = norm((FSubpb1-FSubpb2), 2)/sqrt(length(FSubpb1));
    disp(num2str(UpdateF));
end
disp('====== |delta u^| ======');
for ALSolveStep = 2:ALSolveStep1
    USubpb2_Old = load(['Subpb2_step',num2str(ALSolveStep-1)],'USubpb2'); USubpb2_Old=USubpb2_Old.USubpb2(:);
    USubpb2_New = load(['Subpb2_step',num2str(ALSolveStep)],'USubpb2'); USubpb2_New=USubpb2_New.USubpb2(:);
    UpdateY = norm((USubpb2_Old-USubpb2_New), 2)/sqrt(length(USubpb2_New));
    disp(num2str(UpdateY));
end
disp('====== |delta dual var udual| ======');
for ALSolveStep = 2:ALSolveStep1
    uvdual_Old = load(['uvdual_step',num2str(ALSolveStep-1)],'udual'); uvdual_Old=uvdual_Old.udual(:);
    uvdual_New = load(['uvdual_step',num2str(ALSolveStep)],'udual'); uvdual_New=uvdual_New.udual(:);
    UpdateW = norm((uvdual_Old-uvdual_New), 2)/sqrt(length(uvdual_New));
    disp(num2str(UpdateW));
end
disp('====== |delta dual var vdual| ======');
for ALSolveStep = 2:ALSolveStep1
    uvdual_Old = load(['uvdual_step',num2str(ALSolveStep-1)],'vdual'); uvdual_Old=uvdual_Old.vdual(:);
    uvdual_New = load(['uvdual_step',num2str(ALSolveStep)],'vdual'); uvdual_New=uvdual_New.vdual(:);
    Updatev = norm((uvdual_Old-uvdual_New), 2)/sqrt(length(uvdual_New));
    disp(num2str(Updatev));
end
fprintf('------------ Section 7 Done ------------ \n \n')

% ------ Uncomment these codes to delete temp files ------
% for tempi = 1:ALSolveStep
%     file_name_Subpb1 = ['Subpb1_step',num2str(tempi),'.mat'];
%     file_name_Subpb2 = ['Subpb2_step',num2str(tempi),'.mat'];
%     file_name_dual = ['uvdual_step',num2str(tempi),'.mat'];
%     delete(file_name_Subpb1); delete(file_name_Subpb2); delete(file_name_dual);
% end
% 
% % ------ clear temp variables ------
% clear a ALSub1BadPtNum ALSub1Timetemp atemp b btemp cc ConvItPerEletemp hbar Hbar

 
%% Section 8
fprintf('------------ Section 8 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute strain and plot figures
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Smooth displacements ------
DVCpara.DoYouWantToSmoothOnceMore = funParaInput('SmoothDispOrNot');
% ------ Choose strain computation method ------
DVCpara.MethodToComputeStrain = funParaInput('StrainMethodOp'); 
% ------ Choose strain type (infinitesimal, Eulerian, Green-Lagrangian) ------
DVCpara.StrainType = funParaInput('StrainType');
% ------ Plot displacement & strain components individually or all together ------
DVCpara.PlotComponentEachOrAll = funParaInput('PlotComponentEachOrAll');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Start plotting part -----
for ImgSeqNum = 2:length(ImgNormalized)
    
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DVCpara.ImgSeqIncUnit)-1;
    if DVCpara.ImgSeqIncUnit > 1
        FEMeshIndLast = floor(fNormalizedNewIndex/DVCpara.ImgSeqIncUnit);
    elseif DVCpara.ImgSeqIncUnit == 1
        FEMeshIndLast = floor(fNormalizedNewIndex/DVCpara.ImgSeqIncUnit)-1;
    end
    FEMeshInd = FEMeshIndLast + 1;
    
    if FEMeshInd == 1
        USubpb2 = ResultDisp{ImgSeqNum-1}.U; %+ ResultDisp{10}.U + ResultDisp{20}.U;
        coordinatesFEM = ResultFEMesh{1}.coordinatesFEM; 
        elementsFEM = ResultFEMesh{1}.elementsFEM;
        if (ImgSeqNum-1 == 1) || (DVCpara.ImgSeqIncROIUpdateOrNot==1), UFEMesh = 0*USubpb2; end
    else
        USubpb2 = ResultDisp{ImgSeqNum-1}.U;
        if mod(ImgSeqNum-2,DVCpara.ImgSeqIncUnit) == 0
            coordinatesFEM = ResultFEMesh{FEMeshInd}.coordinatesFEM;
            elementsFEM = ResultFEMesh{FEMeshInd}.elementsFEM;
            coordinatesFEMLast = ResultFEMesh{FEMeshIndLast}.coordinatesFEM;
            UFEMeshLast = ResultDisp{ImgSeqNum-2}.U + UFEMesh;
            xq = coordinatesFEM(:,1); yq = coordinatesFEM(:,2);
            UFEMesh = 0*USubpb2;
            UFEMesh(1:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(1:2:end),xq,yq,'v4');
            UFEMesh(2:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(2:2:end),xq,yq,'v4');
        end
        USubpb2 = USubpb2 + UFEMesh;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %USubpb2 = ResultDisp{ImgSeqNum-1}.U;
    FSubpb2 = ResultDefGrad{ImgSeqNum-1}.F;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xList = min(coordinatesFEM(:,1)):DVCpara.winstepsize(1):max(coordinatesFEM(:,1)); M = length(xList);
    yList = min(coordinatesFEM(:,2)):DVCpara.winstepsize(2):max(coordinatesFEM(:,2)); N = length(yList);
    zList = min(coordinatesFEM(:,3)):DVCpara.winstepsize(3):max(coordinatesFEM(:,3)); L = length(zList);
    [xGrid,yGrid,zGrid] = ndgrid(xList,yList,zList);
    xGrid = xGrid-reshape(UFEMesh(1:3:end),size(xGrid));
    yGrid = yGrid-reshape(UFEMesh(2:3:end),size(yGrid));
    zGrid = zGrid-reshape(UFEMesh(3:3:end),size(zGrid));
    xyz0.x = xGrid; xyz0.y = yGrid; xyz0.z = zGrid;
 
    if size(USubpb2,1)==1
        ULocal = full(USubpb2_New.USubpb2); FLocal = full(FSubpb2.FSubpb2); 
    else
        ULocal = full(USubpb2); FLocal = full(FSubpb2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Smooth displacements ------
    %Plotdisp_show3(full(ULocal),coordinatesFEM,elementsFEM);
    % prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    % DoYouWantToSmoothOnceMore = input(prompt); DispFilterSize=0; DispFilterStd=1;
    SmoothTimes=0;
    try
        while DoYouWantToSmoothOnceMore==0 && SmoothTimes<3
            ULocal = funSmoothDisp3(ULocal,DVCmesh,DVCpara);
            % close all; Plotdisp_show3(full(ULocal),coordinatesFEM,elementsFEM); % Plotuv(ULocal,x0,y0); 
            SmoothTimes = SmoothTimes + 1; %DoYouWantToSmoothOnceMore = input(prompt);
        end
    catch
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ----- Compute strain field ------
    ComputeStrain3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ----- Plot results -----
    close all;
    
    % ------ Plot disp ------
    Plotdisp03(ULocal,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM,DVCpara.PlotComponentEachOrAll);
     
    % ------ Plot strain ------
    %Plotstrain_show3(FLocal,coordinatesFEM,elementsFEM);
    Plotstrain03(full(FStraintemp),xyz0.x(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3)), ...
        xyz0.y(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3)), ...
        xyz0.z(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3)),size(Img{1}),DVCpara.PlotComponentEachOrAll);
    
    % ------ Store strain data ------
    ResultStrain{ImgSeqNum-1}.Strain = FStraintemp;
    tempx = xyz0.x(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3));
    tempy = xyz0.y(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3));
    tempz = xyz0.z(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3));
    ResultStrain{ImgSeqNum-1}.coordinatesFEMStrain = [tempx(:), tempy(:), tempz(:)]; 
    
    % % ------- Add filter and plot strain field -------
    % Plotstrain_Fij;
    % %caxis auto; load('colormap_RdYlBu.mat'); colormap(cMap)
    % Plotstrain03(FStraintemp,xyz0.x(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad),xyz0.y(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad), ...
    %    xyz0.z(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad),size(Img{1}),DVCpara.PlotComponentEachOrAll);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Save figures ------
    % Write down your own codes to save figures! E.g.: % print(['fig_dispu'],'-dpdf');
    if strcmp(DVCpara.PlotComponentEachOrAll,'All')==1
        figure(1); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_disp.fig']);
        figure(2); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_strain.fig']);
    elseif strcmp(DVCpara.PlotComponentEachOrAll,'Individual')==1
        figure(1); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_dispx.fig']);
        figure(2); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_dispy.fig']);
        figure(3); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_dispz.fig']);
        figure(4); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_strainexx.fig']);
        figure(5); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_straineyy.fig']);
        figure(6); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_strainezz.fig']);
        figure(7); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_strainexy.fig']);
        figure(8); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_strainexz.fig']);
        figure(9); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_straineyz.fig']);
    else
        disp('=== Wrong input in DVCpara.PlotComponentEachOrAll! ===')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
fprintf('------------ Section 8 Done ------------ \n \n')


%% Save data again including the computed strain 
results_name = ['results_ws',num2str(DVCpara.winsize(1)),'_st',num2str(DVCpara.winstepsize(1)),'.mat'];
save(results_name, 'file_name','DVCpara','DVCmesh','ResultDisp','ResultDefGrad','ResultStrain','ResultFEMesh',...
    'ALSub1Time','ALSub2Time','ALSolveStep','Resultmubeta','ResultConvItPerEle', ...
    'ResultcoordinatesFEM','ResultelementsFEM','ResultSizeOfFFTSearchReg');

  
%% %%%%%%%%%%%%% Extensions for body and slice plottings %%%%%%%%%%%%%%%%
disp('Extensions for body and slice plottings'); pause;  
plotExt_bodyslice; % Feel free to modify this file (./PlotFiles/plotExt_bodyslice.m) on your purpose.







