function [fileName,Img,DVCpara] = ReadImageLarge3(varargin)
% [file_name,f,g,winsize,winstepsize,gridRange] = ReadImageLarge3(fileName) is
% the main function that loads 3D volumetric images
%
% INPUTS
% -------------------------------------------------------------------------
%   fileName: string for the fileName prefix for the volumetric images in 
%             the current directory.  
%             Input options:
%             --- If image is not within a cell) ---
%             1) 'fileName*.mat' or 'fileName*' 
% 
%             --- If image is within a cell that contains multichannels ---
%             2) fileName{1} = 'fileName*.mat' or 'fileName*' and
%                fileName{2} = channel number containing images you want to
%                              run IDVC on.
%                (if the channel is not provided, i.e. length(fileName) = 1
%                , then channel = 1
%   define_VOI_or_not: whether to define a new VOI (volume of interest) or not
%                 
% OUTPUTS
% -------------------------------------------------------------------------
%   file_name: Store 3D volumetric images file names
%   Img: reference image;
%       Img{1} = f: reference image
%   
%   DVC paramters:
%       winsize: window size for local subset DVC
%       winstepsize: step size between two neighboring local subsets
%       gridRange: ZOI domain range for 3D volumetric images
%                  gridRange{1} = gridxRange in x-direction
%                  gridRange{2} = gridyRange in y-direction
%                  gridRange{3} = gridzRange in z-direction
%
% NOTES
% -------------------------------------------------------------------------
% none
% 
% For more information please see 
%

% fprintf('Choose method to load images:  \n')
% fprintf('     0: Select images folder;  \n')
% fprintf('     1: Use prefix of image names;  \n')
% fprintf('     2: Manually select images.  \n')
% prompt = 'Input here: ';
% LoadImgMethod = input(prompt);


%% ---- Opening & Reading the First Image into CPU Memory ----
[fileInfo,defineVoiOrNot] = parseInputs(varargin{:});
fileName = fileInfo.fileName;
Img{1} = loadFile(fileInfo,1);

%% ---- Opening and Reading Subsequent Images ---
numImages = length(fileInfo.fileName);
% for i = 2:numImages % Reads Volumes Starting on the Second Volumes
%     Img{i} = loadFile(fileInfo,i);
% end

% ==============================================
% Choose VOI
if defineVoiOrNot == 1
    
    fprintf('--- Open original images to define DVC domain? ---\n');
    prompt = 'Input here (0:Open; 1:Not Open): ';
    openImgOrNot = input(prompt); 
    try 
        openImgOrNot(1); 
    catch
        prompt = 'Input here(0:Open; 1:Not Open): ';
        openImgOrNot = input(prompt); 
    end
    switch openImgOrNot 
        case 0
            disp('--- Define top-left and bottom-right corner points on the xy-plane ---')
            MNL = size(Img{1});
            figure, imshow(squeeze(Img{1}(:,:,round(0.5*(1+MNL(3)))))'); 
            caxis auto; axis on; %Show the middle plane of the volumetric image stack
            title('Click top-left and bottom-right corner points on the xy-plane','fontweight','normal','fontsize',14);

            gridx = zeros(1,2); gridy = zeros(1,2); gridz = zeros(1,2);
            [gridx(1), gridy(1)] = ginput(1); fprintf('xy-coordinates of top-left corner point are (%4.3f,%4.3f) [unit: voxel]\n',gridx(1), gridy(1))
            [gridx(2), gridy(2)] = ginput(1); fprintf('xy-coordinates of bottom-right corner point are (%4.3f,%4.3f) [unit: voxel]\n',gridx(2), gridy(2))

            disp('--- Click top and bottom edges on the yz-plane ---')
            figure, imshow(squeeze(Img{1}(round(0.5*(1+MNL(1))),:,:))'); 
            caxis auto; axis on; %Show the middle plane of the volumetric image stack
            title('Click top and bottom edges on the yz-plane','fontweight','normal','fontsize',14);

            [~,gridz(1)] = ginput(1); fprintf('z-coordinate of the top edge is (%4.3f) [unit: voxel]\n',gridz(1));
            [~,gridz(2)] = ginput(1); fprintf('z-coordinate of the bottom edge is (%4.3f) [unit: voxel]\n',gridz(2));
 
            gridx = round(gridx); gridy = round(gridy); gridz = round(gridz);
            
        otherwise
            gridx = [1,size(Img{1},1)];
            gridy = [1,size(Img{1},2)];
            gridz = [1,size(Img{1},3)];
    end

    gridRange = struct('gridxRange',gridx,'gridyRange',gridy,'gridzRange',gridz);

end


% Choose subset size
fprintf('\n');
fprintf('--- What is subvolume size (unit: voxel)? --- \n');
fprintf("E.g. [20,20,10], or '10' if each subvolume is cubic \n");
prompt = 'Input here: ';
winsize = input(prompt);
if length(winsize) == 1, winsize = winsize*ones(1,3); end

% Choose subset size
fprintf('--- What is subvolume step (unit: voxel)? --- \n');
fprintf("E.g. [10,10,5], or '5' if the step is uniform \n");
prompt = 'Input here: ';
winstepsize = input(prompt);
if length(winstepsize) == 1, winstepsize = winstepsize*ones(1,3); end

% ==============================================
% Initial guess method:
initFFTMethod = funParaInput('initFFTMethod');

% ==============================================
% Subproblem 2 solver: finite difference or finite element
Subpb2FDOrFEM = funParaInput('Subpb2FDOrFEM'); % Subproblem 2 using finite difference or fem?
 
% ==============================================
% Parallel cluster #
clusterNo = funParaInput('clusterNo'); % Assign parpool cluster No



% ==============================================
% Deal with image sequence
if numImages > 2 %More than two frames
    
    % ==============================================
    % Decide DVC as accumulative or incremental mode?
    trackingMode = funParaInput('trackingMode');

    % ==============================================
    % DVC initial guess 
    initFFTMethod = funParaInput('initFFTMethod');

    % ==============================================
    if strcmp(trackingMode,'cumulative')
        newFFTSearch = funParaInput('newFFTSearch'); %0-Use previous frame's results; 1-Calculate a new initial guess
    elseif strcmp(trackingMode,'incremental')
        newFFTSearch = 1; %0-Use previous frame's results; 1-Calculate a new initial guess
    end
     
% ================================    
else %Only two frames
    
    trackingMode = 'incremental';
    newFFTSearch = 1; %0-Use previous frame's results; 1-Calculate a new initial guess
     
end
  
DVCpara.winsize = winsize;
DVCpara.winstepsize = winstepsize;
DVCpara.gridRange = gridRange;
DVCpara.Subpb2FDOrFEM = Subpb2FDOrFEM;
DVCpara.clusterNo = clusterNo;
DVCpara.imgSize = size(Img{1});
DVCpara.trackingMode = trackingMode;
DVCpara.initFFTMethod = initFFTMethod;
DVCpara.newFFTSearch = newFFTSearch;
DVCpara.DIM = 3;
 

close all;


end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supplementary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = loadFile(fileInfo,idx)
I = load(fileInfo.fileName{idx});
fieldName = fieldnames(I);
I = getfield(I,fieldName{1});
if iscell(I)
    if numel(I), I = I{1};
    else
        I = I{fileInfo.dataChannel};
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = parseInputs(varargin) % = parseInputs(fileName, redefine_voi_or_not)
 
% Parse fileNames
fileName = varargin{1};
if iscell(fileName)
    if length(fileName)==1, fileInfo.datachannel = 1;
    else fileInfo.datachannel = fileName{2};
    end
    fileName = fileName{1};
end

[~,fileName,~] = fileparts(fileName);
fileName = dir([fileName,'.mat']);
fileInfo.fileName = {fileName.name};

if isempty(fileInfo), error('File name doesn''t exist'); end

defineVOIOrNot = varargin{2};

% Outputs
varargout{1} = fileInfo;
varargout{2} = defineVOIOrNot;
 
end


