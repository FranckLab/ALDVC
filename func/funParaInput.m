function paraInput = funParaInput(paraName)
 
switch paraName
    
    case 'initFFTMethod'
        fprintf('\n');
        fprintf('--- Select initial guess method --- \n')
            fprintf("    1: 'bigxcorr' (by default) \n")
            fprintf('       No prior info of the unknown deformation field; \n')
            fprintf('       Multiscale zero normalized cross correlation to guess large deformations; \n')
            fprintf("    2: 'xcorr' \n")
            fprintf('       Zero normalized cross correlation to guess small deformations; \n')
            fprintf("    3: 'phasecorr' \n")
            fprintf('       Phase correlation for small deformation; \n')
            fprintf("    4: 'bigxcorrUni' \n")
            fprintf('       Zero normalized cross correlation to guess large deformations (robust but expensive!); \n')
            fprintf("    5: 'bigphasecorrUni' \n")
            fprintf('       Phase cross correlation to guess large deformations (robust but expensive!); \n')
        prompt = 'Input here (choose from 1-5): '; initFFTMethod = input(prompt);
        while (initFFTMethod ~= 1) && (initFFTMethod ~= 2) && (initFFTMethod ~= 3) && (initFFTMethod ~= 4) && (initFFTMethod ~= 5) 
            disp('****** Wrong input! ******')
            fprintf('--- Select initial guess method --- \n')
            fprintf("    1: 'bigxcorr' (by default) \n")
            fprintf('       No prior info of the unknown deformation field; \n')
            fprintf('       Multiscale zero normalized cross correlation to guess large deformations; \n')
            fprintf("    2: 'xcorr' \n")
            fprintf('       Zero normalized cross correlation to guess small deformations; \n')
            fprintf("    3: 'phasecorr' \n")
            fprintf('       Phase correlation for small deformation; \n')
            fprintf("    4: 'bigxcorrUni' \n")
            fprintf('       Zero normalized cross correlation to guess large deformations (robust but expensive!); \n')
            fprintf("    5: 'bigphasecorrUni' \n")
            fprintf('       Phase cross correlation to guess large deformations (robust but expensive!); \n')
            prompt = 'Input here (choose from 1-5): '; initFFTMethod = input(prompt);
        end
        switch initFFTMethod
            case 1
                paraInput = 'bigxcorr';
            case 2
                paraInput = 'xcorr';
            case 3
                paraInput = 'phasecorr';
            case 4
                paraInput = 'bigxcorrUni';
            case 5
                paraInput = 'bigphasecorrUni';
            otherwise
                paraInput = 'bigxcorr';
        end
        
        
    case 'newFFTSearch'
        fprintf('\n'); 
        fprintf('Since we are dealing with image sequences, for each new frame,   \n');
        fprintf('do we use last frame result as the initial guess or   \n');
        fprintf('Redo FFT initial guess for every new frame? \n    0: Use last frame (by default); \n    1: Redo initial guess.  \n');
        prompt = 'Input here: '; StillFFTSearch = input(prompt); paraInput = StillFFTSearch;
        fprintf('\n');
        

    case 'Subpb2FDOrFEM'
        fprintf('\n'); 
        fprintf('--- Method to solve ALDVC global step Subproblem 2 ---    \n')
        fprintf('    1: Finite difference (recommended) \n    2: Finite element method (might have edge effects) \n');
        prompt = 'Input here (choose 1 or 2): '; Subpb2FDOrFEM = input(prompt); 
        while (Subpb2FDOrFEM ~= 1) && (Subpb2FDOrFEM ~= 2)    
            disp('****** Wrong input! ******')
            fprintf('\n'); 
            fprintf('--- Method to solve ALDVC global step Subproblem 2 ---    \n')
            fprintf('    1: Finite difference (recommended) \n    2: Finite element method (might have edge effects) \n');
            prompt = 'Input here (choose 1 or 2): '; Subpb2FDOrFEM = input(prompt); 
        end
        switch Subpb2FDOrFEM
            case 1
                paraInput = 'finitedifference';
            case 2
                paraInput = 'finiteelement';
            otherwise
        end
        
        
    case 'trackingMode'
        fprintf('--- Choose cumulative or incremental mode ---  \n')
        fprintf('     1: Cumulative (by default);  \n')
        fprintf('     2: Incremental (to track large deformations);  \n')
        prompt = 'Input here (choose 1 or 2): '; trackingMode = input(prompt);
        while (trackingMode ~= 1) && (trackingMode ~= 2)
            disp('****** Wrong input! ******')
            fprintf('\n'); 
            fprintf('--- Choose cumulative or incremental mode ---  \n')
            fprintf('     1: Cumulative (by default);  \n')
            fprintf('     2: Incremental (to track large deformations);  \n')
            prompt = 'Input here (choose 1 or 2): '; trackingMode = input(prompt);
        end
        switch trackingMode
            case 1
                paraInput = 'cumulative';
            case 2
                paraInput = 'incremental';
            otherwise
        end

        
    case 'clusterNo'
        fprintf('\n'); disp('--- Set up Parallel pool ---');
        fprintf('How many parallel pools to open? (Put in 1 if no parallel computing) \n');
        prompt = 'Input here: ';
        clusterNo = input(prompt);
        %%%%%% Codes to delete current cluster assignments and create parallel pool on cluster %%%%%% 
        % if clusterNo > 1
        %     delete(gcp); myCluster = parcluster('local'); delete(myCluster.Jobs);
        %     parpool(clusterNo,'SpmdEnabled',false);
        % end
        paraInput = clusterNo;




 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Section 8: Inputs to plot strain results %%%%% 
    case 'ConvertUnit' % Convert units from pixels to the physical world
        fprintf("Convert units from pixels to the physical world. \n");
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf("Results in Section 8 (ResultStrain) will be converted to the physical world unit instead of the pixel unit.  \n");
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf("If you want to keep the pixel unit, enter '[1,1,1]' or '1'. \n")
        prompt = 'Input here (e.g., mm/px, um/px): ';
        um2px = input(prompt);
        if length(um2px)==1, um2px = um2px*ones(1,3); end
        paraInput = um2px;
        fprintf('------------------------------------- \n');

    case 'SmoothDispOrNot' % Smooth displacements or not
        prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
        DoYouWantToSmoothOnceMore = input(prompt);
        paraInput = DoYouWantToSmoothOnceMore;
    
    case 'StrainMethodOp' % Choose strain computation method  
        fprintf('What method to use to compute strain? \n');
        fprintf('    0: Direct output from ALDVC; \n');
        fprintf('    1: Finite difference(Recommended); \n');
        fprintf('    2: Plane fitting; \n');
        fprintf('    3: Finite element; \n');
        prompt = 'Input here: ';
        MethodToComputeStrain = input(prompt);
        while (MethodToComputeStrain ~= 0) && (MethodToComputeStrain ~= 1) && (MethodToComputeStrain ~= 2) && (MethodToComputeStrain ~= 3)
            disp('****** Wrong input! ******')
            fprintf('What method to use to compute strain? \n');
            fprintf('    0: Direct output from ALDIC; \n');
            fprintf('    1: Finite difference(Recommended); \n');
            fprintf('    2: Plane fitting; \n');
            fprintf('    3: Finite element; \n');
            prompt = 'Input here: ';
            MethodToComputeStrain = input(prompt);
        end
        paraInput = [MethodToComputeStrain];
          
    case 'StrainType' % Choose strain computation method again
        fprintf('Infinitesimal stran or finite strain? \n');
        fprintf('    0: Infinitesimal stran; \n');
        fprintf('    1: Eulerian strain; \n');
        fprintf('    2: Green-Lagrangian strain; \n');
        fprintf('    3: Others: code by yourself; \n');
        prompt = 'Input here: ';
        StrainType = input(prompt);
        while (StrainType ~= 0) && (StrainType ~= 1) && (StrainType ~= 2) && (StrainType ~= 3)
            disp('****** Wrong input! ******')
            fprintf('Infinitesimal stran or finite strain? \n');
            fprintf('    0: Infinitesimal stran; \n');
            fprintf('    1: Eluerian strain; \n');
            fprintf('    2: Green-Lagrangian strain; \n');
            fprintf('    3: Others: code by yourself; \n');
            prompt = 'Input here: ';
            StrainType = input(prompt);
        end
        paraInput = [StrainType];
            
    case 'PlotComponentEachOrAll'
        fprintf('Plot displacement & strain components individually or all together? \n');
        fprintf('    0: Plot each component individually; \n');
        fprintf('    1: Plot all the deformation components together; \n');
        prompt = 'Input here: ';
        PlotComponentEachOrAll = input(prompt);
        while (PlotComponentEachOrAll ~= 0) && (PlotComponentEachOrAll ~= 1)
            disp('****** Wrong input! ******')
            fprintf('Plot displacement & strain components individually or all together? \n');
            fprintf('    0: Plot each component individually; \n');
            fprintf('    1: Plot all the deformation components together; \n');
            prompt = 'Input here: ';
            PlotComponentEachOrAll = input(prompt);
        end
        switch PlotComponentEachOrAll
            case 0
                paraInput = 'Individual';
            case 1
                paraInput = 'All';
            otherwise
                disp('****** Wrong input! ******')
        end
        
    case 'SaveFigFormat'
        fprintf('Save figures into different format: \n');
        fprintf('    1: jpeg(Choose transparency 0~1) \n');
        fprintf('    2: pdf(Choose transparency = 1) \n'); 
        fprintf('    3: Others: Edit codes in ./plotFiles/SaveFigFiles.m \n'); 
        prompt = 'Input here: '; MethodToSaveFig = input(prompt);
        paraInput = MethodToSaveFig;
        
    case 'OrigDICImgTransparency'
        fprintf('Define transparency for overlaying original images: \n')
        fprintf('Input a real number between 0(Only original images) \n')
        fprintf('and 1(Non-transparent deformation results).\n')
        prompt = 'Input here(e.g. 0.5): '; OrigDICImgTransparency = input(prompt);
        paraInput = OrigDICImgTransparency;
        
        
    otherwise
end
