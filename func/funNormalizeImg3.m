function [ImgNormalized,gridRange] = funNormalizeImg3(Img,gridRange,normalizeOrNot)
% [ImgNormalized,gridRangeCorrected] = funNormalizeImg3(Img,gridRange) is
% to normalize the loaded 3D volumetric images within itself
%
% INPUTS
% -------------------------------------------------------------------------
%   Img: loaded 3D volumetric images
%   gridRange: Input gridROIRange 
%                 
% OUTPUTS
% -------------------------------------------------------------------------
%   ImgNormalized: normalized reference and deformed image;
%       ImgNormalized{1} = fNormalized: normalized reference image
%       ImgNormalized{2} = gNormalized: normalized deformed image
%   gridRangeCorrected: Corrected ZOI domain range for 3D volumetric images
%       gridRangeCorrected{1} = gridRangeCorrected in x-direction
%       gridRangeCorrected{2} = gridRangeCorrected in y-direction
%       gridRangeCorrected{3} = gridRangeCorrected in z-direction
%
% NOTES
% -------------------------------------------------------------------------
% none
% 
% For more information please see 
%


% ======= To garantuee ROI value is within image boundary =======
if gridRange.gridxRange(1) < 1; gridRange.gridxRange(1) = 1; end
if gridRange.gridxRange(2) > size(Img{1},1); gridRange.gridxRange(2) = size(Img{1},1); end

if gridRange.gridyRange(1) < 1; gridRange.gridyRange(1) = 1; end
if gridRange.gridyRange(2) > size(Img{1},2); gridRange.gridyRange(2) = size(Img{1},2); end

if gridRange.gridzRange(1) < 1; gridRange.gridzRange(1) = 1; end
if gridRange.gridzRange(2) > size(Img{1},3); gridRange.gridzRange(2) = size(Img{1},3); end

% temp = gridRange.gridyRange;
% gridRange.gridyRange = gridRange.gridxRange;
% gridRange.gridxRange = temp;

ImgNormalized = cell(size(Img)); 
 
% ============== Normalize or compress images ==============
if (strcmp(normalizeOrNot,'Original')==1) || (strcmp(normalizeOrNot,'original')==1) || (strcmp(normalizeOrNot,'Orig')==1) || (strcmp(normalizeOrNot,'orig')==1)
    for i = 1:length(Img)
        ImgNormalized{i} = double(Img{i});
    end
else
    for i = 1:length(Img)
        tempf = Img{i}(gridRange.gridxRange(1):gridRange.gridxRange(2), ...
                       gridRange.gridyRange(1):gridRange.gridyRange(2), ...
                       gridRange.gridzRange(1):gridRange.gridzRange(2));
        tempf = double(tempf);
        fAvg = mean(tempf(:)); fstd = std(tempf(:));
        ImgNormalized{i} = double((Img{i}-fAvg)/fstd);
    end
    % % ============== Normalize f and g ==============
    % % Normalize f and g
    % Img{1} = double(Img{1});
    % Img{2} = double(Img{2});
    % tempf = Img{1}(gridRange.gridxRange(1):gridRange.gridxRange(2), ...
    %     gridRange.gridyRange(1):gridRange.gridyRange(2), ...
    %     gridRange.gridzRange(1):gridRange.gridzRange(2));
    % fAvg = mean(tempf(:)); fstd = std(tempf(:));
    % 
    % tempg = Img{2}(gridRange.gridxRange(1):gridRange.gridxRange(2), ...
    %     gridRange.gridyRange(1):gridRange.gridyRange(2), ...
    %     gridRange.gridzRange(1):gridRange.gridzRange(2));
    % gAvg = mean(tempg(:)); gstd = std(tempg(:));
    % 
    % ImgNormalized{1} = (Img{1}-fAvg)/fstd;
    % ImgNormalized{2} = (Img{2}-gAvg)/gstd;
end


% ========= Don't Normalized images and use original images =========
% fNormalized = f; gNormalized = g;

% ========= Do wavelets compression =========
% [fNormalized] = func_compress_dw2d(fNormalized0,'sym4',0.1);
% [gNormalized] = func_compress_dw2d(gNormalized0,'sym4',0.1);