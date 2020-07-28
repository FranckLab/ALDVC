function [file_name,Img,DVCpara] = ReadImage3(varargin)
% [file_name,f,g,winsize,winstepsize,gridRange] = ReadImage3(filename) is
% the main function that loads 3D volumetric images
%
% INPUTS
% -------------------------------------------------------------------------
%   filename: string for the filename prefix for the volumetric images in 
%             the current directory.  
%             Input options:
%             --- If image is not within a cell) ---
%             1) 'filename*.mat' or 'filename*' 
% 
%             --- If image is within a cell that contains multichannels ---
%             2) filename{1} = 'filename*.mat' or 'filename*' and
%                filename{2} = channel number containing images you want to
%                              run IDVC on.
%                (if the channel is not provided, i.e. length(filename) = 1
%                , then channel = 1
%                 
% OUTPUTS
% -------------------------------------------------------------------------
%   file_name: Store 3D volumetric images file names
%   Img: reference and deformed image;
%       Img{1} = f: reference image
%       Img{2} = g: deformed image
%   winsize: window size for local subset DVC
%   winstepsize: step size between two neighboring local subsets
%   gridRange: ZOI domain range for 3D volumetric images
%       gridRange{1} = gridxRange in x-direction
%       gridRange{2} = gridyRange in y-direction
%       gridRange{3} = gridzRange in z-direction
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
[fileInfo,runMode] = parseInputs(varargin{:});
file_name = fileInfo.filename;
Img{1} = loadFile(fileInfo,1);

%% ---- Opening and Reading Subsequent Images ---
numImages = length(fileInfo.filename);
for i = 2:numImages % Reads Volumes Starting on the Second Volumes
    Img{i} = loadFile(fileInfo,i);
end

% ==============================================
% Choose ZOI
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
        figure, imshow(squeeze(Img{1}(:,:,round(0.5*(1+MNL(3)))))'); colormap gray;  
        title('Click top-left and bottom-right corner points on the xy-plane','fontweight','normal','fontsize',16);
        
        gridx = zeros(1,2); gridy = zeros(1,2); gridz = zeros(1,2);
        [gridx(1), gridy(1)] = ginput(1); fprintf('xy-coordinates of top-left corner point are (%4.3f,%4.3f)\n',gridx(1), gridy(1))
        [gridx(2), gridy(2)] = ginput(1); fprintf('xy-coordinates of bottom-right corner point are (%4.3f,%4.3f)\n',gridx(2), gridy(2))
        
        disp('--- Define top and bottom edges on the yz-plane ---')
        figure, imshow(squeeze(Img{1}(round(0.5*(1+MNL(1))),:,:))'); 
        title('Click top and bottom edges on the yz-plane','fontweight','normal','fontsize',16);
        
        [~,gridz(1)] = ginput(1); fprintf('z-coordinate of the top edge is (%4.3f)\n',gridz(1) );
        [~,gridz(2)] = ginput(1); fprintf('z-coordinate of the bottom edge is (%4.3f)\n',gridz(2) );
        

        gridx = round(gridx); gridy = round(gridy); gridz = round(gridz);
    otherwise
        gridx = [1,size(Img{1},1)];
        gridy = [1,size(Img{1},2)];
        gridz = [1,size(Img{1},3)];
end

gridRange = struct('gridxRange',gridx,'gridyRange',gridy,'gridzRange',gridz);

% Choose subset size
fprintf('\n');
fprintf('--- What is the subset size? --- \n');
prompt = 'Input here: ';
winsize = input(prompt);
if length(winsize) == 1, winsize = winsize*ones(1,3); end

% Choose subset size
fprintf('--- What is the subset step? --- \n');
prompt = 'Input here: ';
winstepsize = input(prompt);
if length(winstepsize) == 1, winstepsize = winstepsize*ones(1,3); end

% ==============================================
% Initial guess method:
InitFFTMethod = funParaInput('InitFFTMethod');

% ==============================================
% Subproblem 2 solver: finite difference or finite element
Subpb2FDOrFEM = funParaInput('Subpb2FDOrFEM'); % Subproblem 2 using finite difference or fem?
 
% ==============================================
% Parallel cluster #
ClusterNo = funParaInput('ClusterNo'); % Assign parpool cluster No

% ==============================================
% Deal with image sequence
NewFFTSearch = 1; % By default initialize parameters
if numImages > 2
    
    % ==============================================
    % DVC initial guess 
    % InitFFTMethod = funParaInput('InitFFTMethod');
    NewFFTSearch = funParaInput('NewFFTSearch'); % Use last frame as init guess or not
    
    % ==============================================
    % Decide DVC as accumulative or incremental mode?
    fprintf('--- Choose cumulative or incremental mode ---  \n')
    fprintf('     0: Cumulative(By default);  \n')
    fprintf('     1: Incremental;  \n')
    prompt = 'Input here: ';
    DVCIncOrNot = input(prompt);
    try
        temp=DVCIncOrNot(1)+1;
    catch
        DVCIncOrNot=0;
    end
    
    try
        switch DVCIncOrNot
            case 0
                ImgSeqIncUnit = numImages+1;
                ImgSeqIncROIUpdateOrNot = 1;
            case 1
                fprintf('Incremental mode: How many frames to update reference image once? \n');
                prompt = 'Input here: ';
                ImgSeqIncUnit = input(prompt);
                fprintf('Update ROI at the same time of updating reference image? \n');
                fprintf('    0: Do not update ROI; \n'); 
                fprintf('    1: Manually(Recommended); \n'); 
                fprintf('    2: Automatically; \n'); 
                prompt = 'Input here: ';
                ImgSeqIncROIUpdateOrNot = input(prompt);
            otherwise
                ImgSeqIncUnit = numImages+1;
                ImgSeqIncROIUpdateOrNot = 1;
        end
         
    catch
        ImgSeqIncUnit = numImages+1; 
        ImgSeqIncROIUpdateOrNot = 1;
    end
      
% ================================    
else % Only two frames
    
    ImgSeqIncUnit = numImages+1; 
    ImgSeqIncROIUpdateOrNot = 1;
     
end
  
DVCpara.winsize = winsize;
DVCpara.winstepsize = winstepsize;
DVCpara.gridRange = gridRange;
DVCpara.Subpb2FDOrFEM = Subpb2FDOrFEM;
DVCpara.ClusterNo = ClusterNo;
DVCpara.ImgSize = size(Img{1});
DVCpara.ImgSeqIncUnit = ImgSeqIncUnit;
DVCpara.ImgSeqIncROIUpdateOrNot = ImgSeqIncROIUpdateOrNot;
DVCpara.InitFFTMethod = InitFFTMethod;
DVCpara.NewFFTSearch = NewFFTSearch;
DVCpara.DIM = 3;

% DVCpara.CrackOrNot=0; 
% DVCpara.CrackPath1=[0,0]; 
% DVCpara.CrackPath2=[0,0]; 
% DVCpara.CrackTip=[0,0]; 

close all;


end


%% Supp functions
function I = loadFile(fileInfo,idx)
I = load(fileInfo.filename{idx});
fieldName = fieldnames(I);
I = getfield(I,fieldName{1});
if iscell(I)
    if numel(I), I = I{1};
    else
        I = I{fileInfo.dataChannel};
    end
end
end

function varargout = parseInputs(varargin)
%  = parseInputs(filename, sSize, sSizeMin, runMode)
 
% Parse filenames
filename = varargin{1};
if iscell(filename)
    if length(filename) == 1, fileInfo.datachannel = 1;
    else fileInfo.datachannel = filename{2};
    end
    filename = filename{1};
end


[~,filename,~] = fileparts(filename);
filename = dir([filename,'.mat']);
fileInfo.filename = {filename.name};

if isempty(fileInfo), error('File name doesn''t exist'); end

% % Ensure dimensionality of the subset size
% sSize = varargin{2};
% if numel(sSize) == 1
%     sSize = sSize*[1 1 1];
% elseif numel(sSize) ~=3
%     error('Subset size must be a scalar or a three column array');
% end
% 
% %get minium subset size
% sSizeMin = varargin{3};
%                   
% % Ensure range of subset size
% if min(sSize) < 16 || max(sSize > 128)
%    error('Subset size must be within 16 and 128 pixels');
% end
% 
% % Ensure even subset size
% % if sum(mod(sSize,4)) > 0
% %     error('Subset size must be even');
% % end
% 
% if sum(mod(sSize,32)) ~= 0
%     error('Subset size must be 16, 32, 64, 96, or 128 voxels in each dimension');
% end

% Check run method input
runMode = 'c';
% runMode  = varargin{4};

% if ~(strcmpi(runMode(1),'c') || strcmpi(runMode(1),'i') || strcmpi(runMode(1),'h'))
%     error('Run method must be incremental or cumulative or hybrid');
% end
% 
% % Initial guess of displacement field = [0 0 0];
% u0 = num2cell(zeros(1,3));

% Outputs
varargout{      1} = fileInfo;
% varargout{end + 1} = sSize;
% varargout{end + 1} = sSizeMin;
varargout{end + 1} = runMode;
% varargout{end + 1} = u0;

end


