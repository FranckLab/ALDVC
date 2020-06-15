% ==============================================
% function funImgGradient 3D
% ==============================================

function [Df] = funImgGradient3(Img,varargin)
%[imgfNormalizedbc,imggNormalizedbc,imgSize,DfAxis] = funImgGradient(fNormalized,gNormalized,varargin)

fNormalized = Img; % gNormalized = Img{2};
imgSize = size(fNormalized);
% disp('--- Start to compute image gradients ---');

if nargin > 1
    method = varargin{1};
else
    method = 'sobel';
end

switch method
    case 'Splines_interp'
        
        %     imgfNormalizedbc = Spline2D('bicubic',[1:1:size(fNormalized,1)],[1:1:size(fNormalized,2)],fNormalized);
        %     imggNormalizedbc = Spline2D('bicubic',[1:1:size(gNormalized,1)],[1:1:size(gNormalized,2)],gNormalized);
        %     % [XX,YY] = ndgrid([1:1:size(fNormalized,1)],[1:1:size(fNormalized,2)]);
        %     % DfDxNormalizedbc = imgfNormalizedbc.eval_Dx(XX,YY);
        %     % DfDyNormalizedbc = imgfNormalizedbc.eval_Dy(XX,YY);
        %     DfAxis = [0,size(fNormalized,1)-1,0,size(fNormalized,2)-1];
        
    case 'stencil7'
        
        % Finite difference operator
        imgradientMatrix = [-1/60 3/20 -3/4 0 3/4 -3/20 1/60]';
        %     for tempi = 1:7
        %     hx(:,:,tempi) = [zeros(3,7); -1/60 3/20 -3/4 0 3/4 -3/20 1/60; zeros(3,7);];
        %     hy(:,:,tempi) = [zeros(7,3),[-1/60 3/20 -3/4 0 3/4 -3/20 1/60]',zeros(7,3)];
        %     end
        %     hz(:,:,:) = zeros(7,7,7); hz(4,4,1) = -1/60; hz(4,4,2) = 3/20; hz(4,4,3) = -3/4;
        %     hz(4,4,4) = 0; hz(4,4,5) = 3/4; hz(4,4,6) = -3/20; hz(4,4,7) = 1/60;
        imgradientMatrixz(1,1,:) = [-1/60 3/20 -3/4 0 3/4 -3/20 1/60]';
        
        DfDxStartx = 4; % max([4,coordinatesFEM(1,1)-0.5*winsize-1]);
        DfDxStarty = 4; % max([4,coordinatesFEM(1,2)-0.5*winsize-1]);
        DfDxStartz = 4;
        DfDxEndx = size(fNormalized,1)-3; % min([coordinatesFEM(M*N,1)+0.5*winsize,size(fNormalized,1)-3]);
        DfDxEndy = size(fNormalized,2)-3; % min([coordinatesFEM(M*N,2)+0.5*winsize,size(fNormalized,2)-3]);
        DfDxEndz = size(fNormalized,3)-3;
        % DfDxStartx = coordinates(1,1)-1; DfDxStarty = coordinates(1,2)-1;
        % I = fNormalized( coordinates(1,1)-3:coordinates((M+1)*(N+1),1)+3, coordinates(1,2)-3:coordinates((M+1)*(N+1),2)+3);
        I = fNormalized( DfDxStartx-3:DfDxEndx+3, DfDxStarty-3:DfDxEndy+3, DfDxStartz-3:DfDxEndz+3 );
        
        DfDxNormalizedtemp = imfilter(I, imgradientMatrix);
        DfDyNormalizedtemp = imfilter(I, imgradientMatrix');
        DfDzNormalizedtemp = imfilter(I, imgradientMatrixz);
        
        
        DfDxNormalized = DfDxNormalizedtemp(4:end-3, 4:end-3, 4:end-3);
        DfDyNormalized = DfDyNormalizedtemp(4:end-3, 4:end-3, 4:end-3);
        DfDzNormalized = DfDzNormalizedtemp(4:end-3, 4:end-3, 4:end-3);
        
        DfAxis = [DfDxStartx,DfDxEndx,DfDxStarty,DfDxEndy,DfDxStartz,DfDxEndz]-1;
        Df.DfAxis = DfAxis; Df.imgSize = imgSize;
        
        Df.DfDx = DfDxNormalized; Df.DfDy = DfDyNormalized; Df.DfDz = DfDzNormalized;
        
    otherwise
        
        DfDxStartx = 1; DfDxEndx = size(fNormalized,1);
        DfDxStarty = 1; DfDxEndy = size(fNormalized,2);
        DfDxStartz = 1; DfDxEndz = size(fNormalized,3);
        
        DfAxis = [DfDxStartx,DfDxEndx,DfDxStarty,DfDxEndy,DfDxStartz,DfDxEndz]-1;
        Df.DfAxis = DfAxis; Df.imgSize = imgSize;
        
        [DfDxNormalized,DfDyNormalized,DfDzNormalized] = imgradientxyz(fNormalized);
        Df.DfDx = DfDxNormalized; Df.DfDy = DfDyNormalized; Df.DfDz = DfDzNormalized;
        
end

% disp('--- Computing image gradients done ---');

end