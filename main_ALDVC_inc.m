% ---------------------------------------------------
% Augmented Lagrangian Digital Volume Correlation (ALDVC)
% Author: Jin Yang, Postdoc UW-Madison; PhD '19 Caltech;
% Contact: aldicdvc@gmail.com; jin.yang@austin.utexas.edu
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
filename = 'Time*.mat';  fileFolder = './DVC_images/';
try if isempty(fileFolder)~=1, cd(fileFolder); end; catch; end
[file_name_all,Img,DVCpara] = ReadImageLarge3(filename,1); 
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
ImgSeqLength = length(file_name_all);
ResultDisp = cell(ImgSeqLength-1,1);         ResultDefGrad = cell(ImgSeqLength-1,1);
ResultStrain = cell(ImgSeqLength-1,1);       Resultmubeta = cell(ImgSeqLength-1,1);
ResultConvItPerEle = cell(ImgSeqLength-1,1); ResultcoordinatesFEM = cell(ImgSeqLength-1,1);
ResultelementsFEM = cell(ImgSeqLength-1,1);  ResultSizeOfFFTSearchReg = cell(ImgSeqLength-1,1);
Resultxyz0 = cell(ImgSeqLength-1,1);
ResultFEMesh = cell(ceil((ImgSeqLength-1)/DVCpara.ImgSeqIncUnit),1); % For incremental DIC mode
fprintf('------------ Section 2 Done ------------ \n \n')
 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start each frame in an image sequence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ImgSeqNum = 2 : length(file_name_all) % ImgSeqNum is the image stack #
     
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(file_name_all))]);
    
    %% Section 3: Find initial guess
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to find or update initial guess for ALDVC
    % The key idea is to either to use FFT peak fitting, or use last frame
    % results for the next new frame;
    % Particularly in incremental mode, the reference image can also be updated.
    % fNormalized = ImgNormalized{ImgSeqNum-mod(ImgSeqNum-1,ImgSeqIncUnit)};
    % -----------------------------
    try if isempty(fileFolder)~=1, cd(fileFolder); end; catch; end %%%%% Enter image folder 
    % -----------------------------
    load(file_name_all{ImgSeqNum-1}); Img_temp{1} = vol; %%%%% Load image #1
    [ImgNormalized_temp,~] = funNormalizeImg3(Img_temp,DVCpara.gridRange,'Normalize'); %%%%% Normalize image
    Img{1} = ImgNormalized_temp{1}; 
    % -----------------------------
    load(file_name_all{ImgSeqNum}); Img_temp{1} = vol; %%%%% Load image #2
    [ImgNormalized_temp,~] = funNormalizeImg3(Img_temp,DVCpara.gridRange,'Normalize'); %%%%% Normalize image
    Img{2} = ImgNormalized_temp{1}; clear Img_temp ImgNormalized_temp;
    % -----------------------------
    try if isempty(fileFolder)~=1, cd('../'); end; catch; end %%%%% Back to the main path
    % -----------------------------
    NewFFTSearchCheck = 0; DVCpara.NewFFTSearch = 1; %%%%% Perform as NewFFTSearch since it is "incremental" mode
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
        [DVCpara,DVCmesh] = ReadImageRefUpdate(file_name_all,ImgSeqNum,ResultDisp{ImgSeqNum-2}.U,DVCpara,DVCmesh); % Update reference image if needed;
        U0 = zeros(3*size(DVCmesh.coordinatesFEM,1),1); % PlotuvInit;
        ResultFEMesh{1+floor(Img1NewIndex/DVCpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM',DVCmesh.coordinatesFEM,'elementsFEM',DVCmesh.elementsFEM, ...
            'winsize',DVCpara.winsize,'winstepsize',DVCpara.winstepsize,'gridxyROIRange',DVCpara.gridRange );
    else
        U0 = ResultDisp{ImgSeqNum-2}.U;
    end
    
    % ====== Compute image grayscale gradients ======
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
    while (ALSolveStep < 3)
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
    % Plotdisp03(USubpb2,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM,'Individual');
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
    Resultxyz0{ImgSeqNum-1}.xyz0 = xyz0;
    
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


%% ====== Transform "incremental" displacement fields to "cumulative" displacement fields ======
 
tempx = ResultcoordinatesFEM{1}.coordinatesFEM(:,1);
tempy = ResultcoordinatesFEM{1}.coordinatesFEM(:,2);
tempz = ResultcoordinatesFEM{1}.coordinatesFEM(:,3);
coord = [tempx,tempy,tempz]; coordCurr = coord;  

hbar = waitbar(0,'Calculate cumulative disp from incremental disp');

for ImgSeqNum = 2 : length(file_name_all)
     
    waitbar((ImgSeqNum-1)/(size(file_name_all,2)-1));
    
    tempxyz0 = Resultxyz0{ImgSeqNum-1}.xyz0; 
    MNL = size(tempxyz0.x);
    
    tempx = reshape( ResultcoordinatesFEM{ImgSeqNum-1}.coordinatesFEM(:,1), MNL );
    tempy = reshape( ResultcoordinatesFEM{ImgSeqNum-1}.coordinatesFEM(:,2), MNL );
    tempz = reshape( ResultcoordinatesFEM{ImgSeqNum-1}.coordinatesFEM(:,3), MNL );
     
    tempu = reshape( ResultDisp{ImgSeqNum-1}.U(1:3:end), MNL );
    tempv = reshape( ResultDisp{ImgSeqNum-1}.U(2:3:end), MNL );
    tempw = reshape( ResultDisp{ImgSeqNum-1}.U(3:3:end), MNL );
    
    disp_x = interp3(tempy,tempx,tempz,tempu,coordCurr(:,2),coordCurr(:,1),coordCurr(:,3),'makima');
    disp_y = interp3(tempy,tempx,tempz,tempv,coordCurr(:,2),coordCurr(:,1),coordCurr(:,3),'makima');
    disp_z = interp3(tempy,tempx,tempz,tempw,coordCurr(:,2),coordCurr(:,1),coordCurr(:,3),'makima');
    
    disp_x = inpaint_nans3(reshape(disp_x,MNL),0);
    disp_y = inpaint_nans3(reshape(disp_y,MNL),0);
    disp_z = inpaint_nans3(reshape(disp_z,MNL),0);
 
    coordCurr = coordCurr + [disp_x(:), disp_y(:), disp_z(:)];
    U_cum = (coordCurr - coord)'; U_cum = U_cum(:); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Smooth displacements ------
    %Plotdisp_show3(full(ULocal),coordinatesFEM,elementsFEM);
    % prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    % DoYouWantToSmoothOnceMore = input(prompt); DispFilterSize=0; DispFilterStd=1;
    DoYouWantToSmoothOnceMore = 0;
    SmoothTimes = 0;
    try
        while DoYouWantToSmoothOnceMore==0 && SmoothTimes<3
            U_cum = funSmoothDisp3(U_cum,DVCmesh,DVCpara);
            % close all; Plotdisp_show3(full(ULocal),coordinatesFEM,elementsFEM); % Plotuv(ULocal,x0,y0); 
            SmoothTimes = SmoothTimes + 1; %DoYouWantToSmoothOnceMore = input(prompt);
        end
    catch
    end
    
    ResultDisp{ImgSeqNum-1}.U_cum_store = U_cum; % Store cumulative displacement field
    
end

close(hbar);

close all;
for ImgSeqNum = 20
    Plotdisp03(ResultDisp{ImgSeqNum-1}.U_cum_store, ...
        ResultcoordinatesFEM{ImgSeqNum-1}.coordinatesFEM, ...
        ResultelementsFEM{ImgSeqNum-1}.elementsFEM, 0);
end

 
%% Section 8
fprintf('------------ Section 8 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute strain and plot figures
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Convert units from pixels to the physical world ------
DVCpara.um2px = funParaInput('ConvertUnit');
% ------ Smooth displacements ------
DVCpara.DoYouWantToSmoothOnceMore = funParaInput('SmoothDispOrNot');
% ------ Choose strain computation method ------
DVCpara.MethodToComputeStrain = funParaInput('StrainMethodOp'); 
% ------ Choose strain type (infinitesimal, Eulerian, Green-Lagrangian) ------
DVCpara.StrainType = funParaInput('StrainType');
% ------ Plot displacement & strain components individually or all together ------
DVCpara.PlotComponentEachOrAll = funParaInput('PlotComponentEachOrAll');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------ Start plotting part -----
for ImgSeqNum = 2:length(file_name_all)
    
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(file_name_all))]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ====== Old version ======
    % fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DVCpara.ImgSeqIncUnit)-1;
    % if DVCpara.ImgSeqIncUnit > 1
    %     FEMeshIndLast = floor(fNormalizedNewIndex/DVCpara.ImgSeqIncUnit);
    % elseif DVCpara.ImgSeqIncUnit == 1
    %     FEMeshIndLast = floor(fNormalizedNewIndex/DVCpara.ImgSeqIncUnit)-1;
    % end
    % FEMeshInd = FEMeshIndLast + 1;
    % 
    % if FEMeshInd == 1
    %     USubpb2 = ResultDisp{ImgSeqNum-1}.U; %+ ResultDisp{10}.U + ResultDisp{20}.U;
    %     coordinatesFEM = ResultFEMesh{1}.coordinatesFEM; 
    %     elementsFEM = ResultFEMesh{1}.elementsFEM;
    %     if (ImgSeqNum-1 == 1) || (DVCpara.ImgSeqIncROIUpdateOrNot==1), UFEMesh = 0*USubpb2; end
    % else
    %     USubpb2 = ResultDisp{ImgSeqNum-1}.U;
    %     if mod(ImgSeqNum-2,DVCpara.ImgSeqIncUnit) == 0
    %         coordinatesFEM = ResultFEMesh{FEMeshInd}.coordinatesFEM;
    %         elementsFEM = ResultFEMesh{FEMeshInd}.elementsFEM;
    %         coordinatesFEMLast = ResultFEMesh{FEMeshIndLast}.coordinatesFEM;
    %         UFEMeshLast = ResultDisp{ImgSeqNum-2}.U + UFEMesh;
    %         xq = coordinatesFEM(:,1); yq = coordinatesFEM(:,2);
    %         UFEMesh = 0*USubpb2;
    %         UFEMesh(1:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(1:2:end),xq,yq,'v4');
    %         UFEMesh(2:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(2:2:end),xq,yq,'v4');
    %     end
    %     USubpb2 = USubpb2 + UFEMesh;
    % end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %USubpb2 = ResultDisp{ImgSeqNum-1}.U;
    % FSubpb2 = ResultDefGrad{ImgSeqNum-1}.F;
    % ====== Updated code for incremental mode code ======
    USubpb2 = ResultDisp{ImgSeqNum-1}.U_cum_store;
    coordinatesFEM = ResultcoordinatesFEM{1}.coordinatesFEM;
    elementsFEM = ResultelementsFEM{1}.elementsFEM;
    FLocal = 0*ResultDefGrad{1}.F;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xList = min(coordinatesFEM(:,1)):DVCpara.winstepsize(1):max(coordinatesFEM(:,1)); M = length(xList);
    yList = min(coordinatesFEM(:,2)):DVCpara.winstepsize(2):max(coordinatesFEM(:,2)); N = length(yList);
    zList = min(coordinatesFEM(:,3)):DVCpara.winstepsize(3):max(coordinatesFEM(:,3)); L = length(zList);
    [xGrid,yGrid,zGrid] = ndgrid(xList,yList,zList);
    xyz0.x = xGrid; xyz0.y = yGrid; xyz0.z = zGrid;
 
    if size(USubpb2,1)==1
        ULocal = full(USubpb2_New.USubpb2); % FLocal = full(FSubpb2.FSubpb2); 
    else
        ULocal = full(USubpb2); % FLocal = full(FSubpb2);
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
    % ----- Plot results in image (pixel/voxel) unit -----
    %     close all;
    %     
    %     % ------ Plot disp ------
    %     Plotdisp03(ULocal, coordinatesFEM, elementsFEM, DVCpara.PlotComponentEachOrAll);
    %      
    %     % ------ Plot strain ------
    %     % Plotstrain_show3(FLocal,coordinatesFEM,elementsFEM);
    %     Plotstrain03(full(FStraintemp),xyz0.x(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3)), ...
    %         xyz0.y(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3)), ...
    %         xyz0.z(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3)),DVCpara.ImgSize,DVCpara.PlotComponentEachOrAll);
    %     
    %     % ------ Store strain data ------
    %     ResultStrain{ImgSeqNum-1}.Strain = FStraintemp;
    %     tempx = xyz0.x(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3));
    %     tempy = xyz0.y(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3));
    %     tempz = xyz0.z(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3));
    %     ResultStrain{ImgSeqNum-1}.coordinatesFEMStrain = [tempx(:), tempy(:), tempz(:)]; 
    %     
    %     % ------ Calculate "mean + std" value of strain components ------
    %     FStraintemp11 = reshape(FStraintemp(1:9:end), M-2*Rad(1), N-2*Rad(2), L-2*Rad(3));
    %     FStraintemp11_crop = FStraintemp11(4:end-3, 11:end-10, 2:end-1); % TODO
    %     % figure, imagesc3(FStraintemp11_crop);
    %     FStrain11_mean(ImgSeqNum) = mean(FStraintemp11_crop(:))
    %     FStrain11_std(ImgSeqNum) = std(FStraintemp11_crop(:))
    % 
    %     FStraintemp22 = reshape(FStraintemp(5:9:end), M-2*Rad(1), N-2*Rad(2), L-2*Rad(3));
    %     FStraintemp22_crop = FStraintemp22(4:end-3, 6:end-5, 2:end-1); % TODO
    %     % figure, imagesc3(FStraintemp22_crop); caxis auto; colorbar;
    %     FStrain22_mean(ImgSeqNum) = mean(FStraintemp22_crop(:))
    %     FStrain22_std(ImgSeqNum) = std(FStraintemp22_crop(:))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ----- Plot results in physical world unit -----
    close all;
    UWorld = ([ ULocal(1:3:end), ULocal(2:3:end), ULocal(3:3:end) ] * diag(DVCpara.um2px))';   UWorld = UWorld(:);
    coordinatesFEMWorld = coordinatesFEM * diag(DVCpara.um2px);
    
    FStraintempWorld_11 = FStraintemp(1:9:end); % du/dx
    FStraintempWorld_21 = FStraintemp(2:9:end) * DVCpara.um2px(2) / DVCpara.um2px(1); % dv/dx
    FStraintempWorld_31 = FStraintemp(3:9:end) * DVCpara.um2px(3) / DVCpara.um2px(1); % dw/dx
    FStraintempWorld_12 = FStraintemp(4:9:end) * DVCpara.um2px(1) / DVCpara.um2px(2); % du/dy
    FStraintempWorld_22 = FStraintemp(5:9:end); % dv/dy
    FStraintempWorld_32 = FStraintemp(6:9:end) * DVCpara.um2px(3) / DVCpara.um2px(2); % dw/dy
    FStraintempWorld_13 = FStraintemp(7:9:end) * DVCpara.um2px(1) / DVCpara.um2px(3); % du/dz
    FStraintempWorld_23 = FStraintemp(8:9:end) * DVCpara.um2px(2) / DVCpara.um2px(3); % dv/dz
    FStraintempWorld_33 = FStraintemp(9:9:end); % dw/dz

    FStraintempWorld = ([ FStraintempWorld_11, FStraintempWorld_21, FStraintempWorld_31, ...
                         FStraintempWorld_12, FStraintempWorld_22, FStraintempWorld_32, ...
                         FStraintempWorld_13, FStraintempWorld_23, FStraintempWorld_33 ])';
    FStraintempWorld = FStraintempWorld(:);

    xyz0World.x = xyz0.x * DVCpara.um2px(1);
    xyz0World.y = xyz0.y * DVCpara.um2px(2);
    xyz0World.z = xyz0.z * DVCpara.um2px(3);

    % ------ Plot disp ------
    Plotdisp03(UWorld, coordinatesFEMWorld, elementsFEM, DVCpara.PlotComponentEachOrAll);
     
    % ------ Plot strain ------
    % Plotstrain_show3(FLocal,coordinatesFEM,elementsFEM);
    Plotstrain03(full(FStraintempWorld),xyz0World.x(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3)), ...
        xyz0World.y(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3)), ...
        xyz0World.z(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3)),DVCpara.ImgSize,DVCpara.PlotComponentEachOrAll);
    
    % ------ Store strain data ------
    ResultStrain{ImgSeqNum-1}.StrainWorld = FStraintempWorld;
    tempxWorld = xyz0World.x(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3));
    tempyWorld = xyz0World.y(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3));
    tempzWorld = xyz0World.z(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3));
    ResultStrain{ImgSeqNum-1}.coordinatesFEMStrainWorld = [tempxWorld(:), tempyWorld(:), tempzWorld(:)]; 
    
    % ------ Calculate "mean + std" value of strain components ------
    FStraintemp11 = reshape(FStraintempWorld(1:9:end), M-2*Rad(1), N-2*Rad(2), L-2*Rad(3));
    FStraintemp11_crop = FStraintemp11(4:end-3, 21:end-20, 2:end-1); % TODO   % figure, imagesc3(FStraintemp11_crop);
    Strain11_mean(ImgSeqNum) = mean(FStraintemp11_crop(:))
    Strain11_std(ImgSeqNum) = std(FStraintemp11_crop(:))

    FStraintemp22 = reshape(FStraintempWorld(5:9:end), M-2*Rad(1), N-2*Rad(2), L-2*Rad(3));
    FStraintemp22_crop = FStraintemp22(4:end-3, 21:end-20, 2:end-1); % TODO
    Strain22_mean(ImgSeqNum) = mean(FStraintemp22_crop(:))
    Strain22_std(ImgSeqNum) = std(FStraintemp22_crop(:))
    
    FStraintemp33 = reshape(FStraintempWorld(9:9:end), M-2*Rad(1), N-2*Rad(2), L-2*Rad(3));
    FStraintemp33_crop = FStraintemp33(4:end-3, 21:end-20, 2:end-1); % TODO
    Strain33_mean(ImgSeqNum) = mean(FStraintemp33_crop(:))
    Strain33_std(ImgSeqNum) = std(FStraintemp33_crop(:))
 
    FStraintemp21 = reshape(FStraintempWorld(2:9:end), M-2*Rad(1), N-2*Rad(2), L-2*Rad(3));
    FStraintemp21_crop = FStraintemp21(4:end-3, 21:end-20, 2:end-1); 
    FStraintemp12 = reshape(FStraintempWorld(4:9:end), M-2*Rad(1), N-2*Rad(2), L-2*Rad(3));
    FStraintemp12_crop = FStraintemp12(4:end-3, 21:end-20, 2:end-1); 
    Strain12_mean(ImgSeqNum) = mean(0.5*(FStraintemp12_crop(:)+FStraintemp21_crop(:)))
    Strain12_std(ImgSeqNum) = std(0.5*(FStraintemp12_crop(:)+FStraintemp21_crop(:)))

    FStraintemp31 = reshape(FStraintempWorld(3:9:end), M-2*Rad(1), N-2*Rad(2), L-2*Rad(3));
    FStraintemp31_crop = FStraintemp31(4:end-3, 21:end-20, 2:end-1); 
    FStraintemp13 = reshape(FStraintempWorld(7:9:end), M-2*Rad(1), N-2*Rad(2), L-2*Rad(3));
    FStraintemp13_crop = FStraintemp13(4:end-3, 21:end-20, 2:end-1); 
    Strain13_mean(ImgSeqNum) = mean(0.5*(FStraintemp13_crop(:)+FStraintemp31_crop(:)))
    Strain13_std(ImgSeqNum) = std(0.5*(FStraintemp13_crop(:)+FStraintemp31_crop(:)))

    FStraintemp32 = reshape(FStraintempWorld(6:9:end), M-2*Rad(1), N-2*Rad(2), L-2*Rad(3));
    FStraintemp32_crop = FStraintemp32(4:end-3, 21:end-20, 2:end-1); 
    FStraintemp23 = reshape(FStraintempWorld(8:9:end), M-2*Rad(1), N-2*Rad(2), L-2*Rad(3));
    FStraintemp23_crop = FStraintemp23(4:end-3, 21:end-20, 2:end-1); 
    Strain23_mean(ImgSeqNum) = mean(0.5*(FStraintemp32_crop(:)+FStraintemp23_crop(:)))
    Strain23_std(ImgSeqNum) = std(0.5*(FStraintemp32_crop(:)+FStraintemp23_crop(:)))


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
    'ResultcoordinatesFEM','ResultelementsFEM','ResultSizeOfFFTSearchReg','Resultxyz0');

  
%% %%%%%%%%%%%%% Extensions for body and slice plottings %%%%%%%%%%%%%%%%
disp('Extensions for body and slice plottings'); pause;  
plotExt_bodyslice; % Feel free to modify this file (./PlotFiles/plotExt_bodyslice.m) on your purpose.


%%
figure,   errorbar( [2:36], Strain11_mean(2:36), Strain11_std(2:36),'linewidth',1 );
hold on;  errorbar( [2:36], Strain22_mean(2:36), Strain22_std(2:36),'linewidth',1 );
hold on;  errorbar( [2:36], Strain33_mean(2:36), Strain33_std(2:36),'linewidth',1 );
hold on;  errorbar( [2:36], Strain12_mean(2:36), Strain12_std(2:36),'linewidth',1 );
hold on;  errorbar( [2:36], Strain13_mean(2:36), Strain13_std(2:36),'linewidth',1 );
hold on;  errorbar( [2:36], Strain23_mean(2:36), Strain23_std(2:36),'linewidth',1 );

set(gca,'fontsize',20);
xlabel('Frame #');
ylabel('Strain');

lgd = legend('e11','e22','e33','e12','e13','e23');
set(lgd,'location','northeastoutside'); set(lgd,'fontsize',14);









