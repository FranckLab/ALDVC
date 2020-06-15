function [xyz,uvw,cc] = funIntegerSearch3Mg(Img,gridRange,winsize,winstepsize,fftmethod)

MaxVoxelNum = 100^3; MinSearchWS = 30;

% [xyz,uvw,Phi] = funIntegerSearch3(Img,tempSizeOfSearchRegion,gridRange,winsize,winstepsize,method)
% is the main function that loads 3D volumetric images
%
% INPUTS
% -------------------------------------------------------------------------
%   Img: reference and deformed image;
%       Img{1} = f: reference image
%       Img{2} = g: deformed image
%   tempSizeOfSearchRegion: temp size of search region
%   gridRange: ZOI domain range for 3D volumetric images
%       gridRange{1} = gridxRange in x-direction
%       gridRange{2} = gridyRange in y-direction
%       gridRange{3} = gridzRange in z-direction
%   winsize: window size for local subset DVC
%   winstepsize: step size between two neighboring local subsets
%   method: which correlation function are used:
%       method == 0: cross correlation
%       method == 1: normalized cross correlation
%       method == 2: gradient cross correlation (w/ filter |k|^2 in Fourier sp)
%
% OUTPUTS
% -------------------------------------------------------------------------
%   xyz: assigned coordinates of IntegerSearch3 nodes
%   uvw: computed displacements at xyz0 coordinates
%   Phi: correlation function of used correlation function
%
% NOTES
% -------------------------------------------------------------------------
% none
%
% For more information please see
%


% % Only for the same size template search
% % Pad images with zeros so that we don't grad any subset outside of the
% % image domain. This would produce an error
% padSize = (tempSizeOfSearchRegion+1)*[1,1,1];
% Img{1} = padarray(Img{1},padSize,'replicate','both');
% Img{2} = padarray(Img{2},padSize,'replicate','both');

% gridxBackup = gridx; gridyBackup = gridy; gridzBackup = gridz;

 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch fftmethod
    
    case {'xcorr','phasecorr'}
        
        %% use DVC codes
        sSize = winsize; sSpacing = winstepsize;
        I0{1} = single(Img{1});
        I0{2} = single(Img{2});

        % I0{1} = permute(I0{1},[2 1 3]);
        % I0{2} = permute(I0{2},[2 1 3]);

        sizeI0 = size(I0{1});
        sizeI = ceil(sizeI0./sSpacing).*sSpacing;
        prePad1 = ceil((sizeI - sizeI0)/2);
        postPad1 = floor((sizeI - sizeI0)/2);

        I{1} = padarray(I0{1},prePad1,'replicate','pre');
        I{1} = padarray(I{1},postPad1,'replicate','post');

        I{2} = padarray(I0{2},prePad1,'replicate','pre');
        I{2} = padarray(I{2},postPad1,'replicate','post');

        prePad2 = sSize/2;
        postPad2 = sSize/2;

        sizeI = size(I{1});
        I{1} = padarray(I{1},prePad2,0,'pre');
        I{1} = padarray(I{1},postPad2,0,'post');

        I{2} = padarray(I{2},prePad2,0,'pre');
        I{2} = padarray(I{2},postPad2,0,'post');

        % pad images with zeros so that we don't grab any subset outside of the image
        % domain. This would produce an error
        padSize = ceil(0.5*[winsize(1),winsize(2),winsize(3)]);
        Img{1} = padarray(I{1},padSize,'replicate','both');
        Img{2} = padarray(I{2},padSize,'replicate','both');

        %inject small noise to avoid flat subsets
        noise1 = rand(size(Img{1}))/10000;
        noise2 = rand(size(Img{2}))/10000;
        Img{1} = Img{1}+noise1;
        Img{2} = Img{2}+noise2;

        % % run cross-correlation to get an estimate of the displacements
        % [du,cc] = DVC(I,winsize,winstepsize,winstepsize/2,1.25);
        %
        % m = cell(1,3);
        % for i = 1:3, m{i} = (1:sSpacing(i):(sizeI(i) + 1)) + sSize(i)/2 - winstepsize(i)/2; end
        %
        % [x,y,z] = ndgrid(m{1},m{2},m{3});
        %
        % uvw.u = du{1}; uvw.v = du{2}; uvw.w = du{3};
        % xyz.x = x; xyz.y = y; xyz.z = z;
        %
        % end
        
        padSizeTotal = prePad1+postPad1+prePad2+postPad2+2*padSize;
        gridx = gridRange.gridxRange + [0,padSizeTotal(1)];
        gridy = gridRange.gridyRange + [0,padSizeTotal(2)];
        gridz = gridRange.gridzRange + [0,padSizeTotal(3)];

        % Initialize quadratic least squares fitting coefficients
        [mx, my, mz] = meshgrid((-1:1),(-1:1),(-1:1));
        m = [mx(:), my(:), mz(:)];
        
        M{1} = zeros(size(m,1),10);
        for i = 1:size(m,1)
            x = m(i,1); y = m(i,2); z = m(i,3);
            M{1}(i,:) = [1,x,y,z,x^2,x*y,x*z,y^2,y*z,z^2];
        end
        
        M{2} = M{1}'*M{1};
        
        % Generate Moduluar transfer function (see eq. 3)
        [~,~,MTF] = generateMTF((winsize+[1,1,1]));
        
        % may lose boundary regions
        gridxBackup = gridx; gridyBackup = gridy; gridzBackup = gridz;
        
        % may lose boundary regions
        while gridx(1) < 1+3+0.5*winsize(1)+0.5*padSizeTotal(1)
            gridx(1) = gridx(1) + 1;
        end
        while gridy(1) < 1+3+0.5*winsize(2)+0.5*padSizeTotal(2)
            gridy(1) = gridy(1) + 1;
        end
        while gridz(1) < 1+3+0.5*winsize(3)+0.5*padSizeTotal(3)
            gridz(1) = gridz(1) + 1;
        end
        while gridx(end) + 0.5*winsize(1)+0.5*padSizeTotal(1) > size(Img{1},1)-3
            gridx(end) = gridx(end) - 1;
        end
        while gridy(end) + 0.5*winsize(2)+0.5*padSizeTotal(2) > size(Img{1},2)-3
            gridy(end) = gridy(end) - 1;
        end
        while gridz(end) + 0.5*winsize(3)+0.5*padSizeTotal(3) > size(Img{1},3)-3
            gridz(end) = gridz(end) - 1;
        end
        Img1Partemp = Img{1};
        Img2Partemp = Img{2};
        
    otherwise % search within certain domain
        
        Img1Partemp = Img{1};
        Img2Partemp = Img{2};
        gridx = gridRange.gridxRange ;
        gridy = gridRange.gridyRange ;
        gridz = gridRange.gridzRange ;
        
        % may lose boundary regions, backup first
        gridxBackup = gridx; gridyBackup = gridy; gridzBackup = gridz;
        borderGap = 1+round(1*winsize); % 1*3 array
        % Cut border regions a little bit
        if gridx(1)<borderGap(1), gridx(1)=borderGap(1); end
        if gridy(1)<borderGap(2), gridy(1)=borderGap(2); end
        if gridz(1)<borderGap(3), gridz(1)=borderGap(3); end
        if gridx(2)>size(Img1Partemp,1)-borderGap(1), gridx(2)=size(Img1Partemp,1)-borderGap(1); end
        if gridy(2)>size(Img1Partemp,2)-borderGap(2), gridy(2)=size(Img1Partemp,2)-borderGap(2); end
        if gridz(2)>size(Img1Partemp,3)-borderGap(3), gridz(2)=size(Img1Partemp,3)-borderGap(3); end
        
end



% disp('Assemble point position sequence.');
XList = [gridx(1) : winstepsize(1) : gridx(end)];
YList = [gridy(1) : winstepsize(2) : gridy(end)];
ZList = [gridz(1) : winstepsize(3) : gridz(end)];

[XX,YY,ZZ] = ndgrid(XList,YList,ZList);
temparrayLength = length(XList)*length(YList)*length(ZList);
PtPosSeq = zeros(temparrayLength, 3);
PtPosSeq(:,1) = XX(:); PtPosSeq(:,2) = YY(:); PtPosSeq(:,3) = ZZ(:);

u = zeros(length(XList),length(YList),length(ZList));
v = u; w = u; Phi = u; x = u; y = u; z = u;
ck1temp = zeros(temparrayLength,1); cj1temp = ck1temp; ci1temp = ck1temp; Phitemp = ck1temp;
utemp = ck1temp; vtemp = ck1temp; wtemp = ck1temp;
xtemp = ck1temp; ytemp = ck1temp; ztemp = ck1temp;


% prompt = 'How many parallel pools to open? (Put in 0 if no parallel computing)';
% delete(gcp)
% ClusterNo = input(prompt);
% parpool(ClusterNo);
% ------ Choose local fft solver ------
switch fftmethod
    %% --------------------------------------------------------------------------
    case {'xcorr','phasecorr'}
        hbar = waitbar(0,'Integer search, please drink coffee and wait.');
        % hbar = parfor_progressbar(temparrayLength,'Please wait for integer search!');

        for tempi = 1 : temparrayLength
            
            ii = PtPosSeq(tempi,1); jj = PtPosSeq(tempi,2); kk = PtPosSeq(tempi,3);
            
            waitbar(tempi/temparrayLength);
            % hbar.iterate(1);
            
            C = Img1Partemp(ii:ii+winsize(1), jj:jj+winsize(2), kk:kk+winsize(3));
            
            D = Img2Partemp(ii:ii+winsize(1), jj:jj+winsize(2), kk:kk+winsize(3));
            C = MTF.*C; D = MTF.*D;
            if strcmp(fftmethod,'xcorr')
                A = xCorr3(C,D,(winsize+[1,1,1]));
                %close all
                %figure,imagesc3D(A);
            elseif strcmp(fftmethod,'phasecorr')
                A = phaseCorr3(C,D,(winsize+[1,1,1]));
            end
            % find maximum index of the cross-correlaiton
            [max_f, maxIdx] = max(A(:));
            % compute voxel resolution displacements
            [u1, u2, u3] = ind2sub((winsize+[1,1,1]),maxIdx);
            
            %find qfactors
            %cc.A{1} = A;
            %qfactors(tempi,:) = compute_qFactor(cc,tempi);
            qfactors(tempi,:) = compute_qFactorMat(A,tempi);
            
            % gather the 3x3x3 voxel neighborhood around the peak
            try xCorrPeak = reshape(A(u1 + (-1:1), u2 + (-1:1), u3 + (-1:1)),27,1);
                % last squares fitting of the peak to calculate sub-voxel displacements
                du123 = lsqPolyFit3(xCorrPeak, M{1}, M{2});
                % u123(k,:) = [u1 u2 u3] + du123' - sSize/2 - 1;
                utemp(tempi) = -( u1+du123(1) - 0.5*(winsize(1)) - 1 );
                vtemp(tempi) = -( u2+du123(2) - 0.5*(winsize(2)) - 1 );
                wtemp(tempi) = -( u3+du123(3) - 0.5*(winsize(3)) - 1 );
                %--------------------------------------------------------------------------
            catch
                % u123(k,:) = nan;
                utemp(tempi) = nan; vtemp(tempi) = nan; wtemp(tempi) = nan;
            end
            
            xtemp(tempi) = (ii+ii+winsize(1))/2 - padSize(1);
            ytemp(tempi) = (jj+jj+winsize(2))/2 - padSize(2);
            ztemp(tempi) = (kk+kk+winsize(3))/2 - padSize(3);
            
            ck1temp(tempi) = ceil(tempi/(length(YList)*length(XList)));
            cj1temp(tempi) = ceil((tempi-(ck1temp(tempi)-1)*(length(YList)*length(XList)))/length(XList));
            ci1temp(tempi) = tempi - (ck1temp(tempi)-1)*(length(YList)*length(XList)) - (cj1temp(tempi)-1)*length(XList);
            
            Phitemp(tempi) = max_f;
            
        end
        close(hbar); % toc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for k = 1:2
        %     qf_ = (qfactors(:,k)-min(qfactors(:,k)));
        %     cc.qfactors(:,k) = qf_/max(qf_);
        % end
        % hbar = parfor_progressbar(temparrayLegTotal,'Assign results to variables.');
        hbar = waitbar(0,'Assign results to variables.');
        for tempi = 1:temparrayLength
            
            ci1 = ci1temp(tempi); cj1 = cj1temp(tempi); ck1 = ck1temp(tempi);
            
            u(ci1,cj1,ck1) = utemp(tempi);
            v(ci1,cj1,ck1) = vtemp(tempi);
            w(ci1,cj1,ck1) = wtemp(tempi);
            
            Phi(ci1,cj1,ck1) = Phitemp(tempi);
            
            x(ci1,cj1,ck1) = xtemp(tempi);
            y(ci1,cj1,ck1) = ytemp(tempi);
            z(ci1,cj1,ck1) = ztemp(tempi);
            
            waitbar(tempi/temparrayLength);
            % hbar.iterate(1);
        end
        close(hbar);
        
         
        xyz.x = x-floor(0.5*padSizeTotal(1));
        xyz.y = y-floor(0.5*padSizeTotal(2));
        xyz.z = z-floor(0.5*padSizeTotal(3));
        uvw.u = u; uvw.v = v; uvw.w = w;
        cc.max = Phi;
        cc.A = [];

        
        %% --------------------------------------------------------------------------
    case {'bigphasecorr','bigxcorr'}
        
        % Make sure image width/length odd number
        if mod((gridx(2)-gridx(1)),2)==1, gridx(2)=gridx(2)-1;  end
        if mod((gridy(2)-gridy(1)),2)==1, gridy(2)=gridy(2)-1; end
        if mod((gridz(2)-gridz(1)),2)==1, gridz(2)=gridz(2)-1; end
          
        % Level 1: first level cross-correlation
        C = Img1Partemp(gridx(1):gridx(2),gridy(1):gridy(2),gridz(1):gridz(2));
        D = Img2Partemp(gridx(1):gridx(2),gridy(1):gridy(2),gridz(1):gridz(2));
        % Shrink size C & D
        RDTime = 1;
        while size(C,1)*size(C,2)*size(C,3)>MaxVoxelNum
            RDTime = RDTime*2;
            C = imgaussfilt3(C); % C = imgaussfilt3(C); C = imgaussfilt3(C);
            D = imgaussfilt3(D); % D = imgaussfilt3(D); D = imgaussfilt3(D);
            CNew = C(1:2:2*floor(size(C,1)/2),1:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                C(2:2:2*floor(size(C,1)/2), 1:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                C(1:2:2*floor(size(C,1)/2),2:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                C(2:2:2*floor(size(C,1)/2), 2:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                C(1:2:2*floor(size(C,1)/2),1:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)) + ...
                C(2:2:2*floor(size(C,1)/2), 1:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)) + ...
                C(1:2:2*floor(size(C,1)/2),2:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)) + ...
                C(2:2:2*floor(size(C,1)/2), 2:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2));
            DNew = D(1:2:2*floor(size(C,1)/2),1:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                D(2:2:2*floor(size(C,1)/2), 1:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                D(1:2:2*floor(size(C,1)/2),2:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                D(2:2:2*floor(size(C,1)/2), 2:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                D(1:2:2*floor(size(C,1)/2),1:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)) + ...
                D(2:2:2*floor(size(C,1)/2), 1:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)) + ...
                D(1:2:2*floor(size(C,1)/2),2:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)) + ...
                D(2:2:2*floor(size(C,1)/2), 2:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)); 
            C=CNew; D=DNew;
        end
        % figure, imagesc3D(C); figure, imagesc3D(D);
         
        if strcmp(fftmethod,'bigxcorr')
            XCORRF3OfCD0 = normxcorr3(C,D,'full');
        elseif strcmp(fftmethod,'bigphasecorr')
            XCORRF3OfCD0 = bigPhaseCorr3(C,D,'full');
        end
        % figure, imagesc3D(XCORRF3OfCD0);
        
        % find maximum index of the cross-correlaiton
        [v1temp, u1temp, w1temp, max_f] = findpeak3(real(XCORRF3OfCD0),1);
        zero_disp = ceil((size(XCORRF3OfCD0)+[1,1,1])/2);
               
        utemp = RDTime*(u1temp-zero_disp(1));
        vtemp = RDTime*(v1temp-zero_disp(2));
        wtemp = RDTime*(w1temp-zero_disp(3));
        Phitemp = max_f;
        
        % Check: out of image border or not
        gridx_f = gridx; gridy_f = gridy; gridz_f = gridz;
        gridx_g = gridx+repmat(ceil(utemp),1,2); gridy_g = gridy+repmat(ceil(vtemp),1,2);
        gridz_g = gridz+repmat(ceil(wtemp),1,2); %updated ROI in g, could be out of image
        
        if gridx_g(1)<borderGap(1), temp=gridx_g(1); gridx_g(1)=borderGap(1); gridx_f(1)=gridx_f(1)+borderGap(1)-temp; end
        if gridy_g(1)<borderGap(2), temp=gridy_g(1); gridy_g(1)=borderGap(2); gridy_f(1)=gridy_f(1)+borderGap(2)-temp; end
        if gridz_g(1)<borderGap(3), temp=gridz_g(1); gridz_g(1)=borderGap(3); gridz_f(1)=gridz_f(1)+borderGap(3)-temp; end
        
        if gridx_g(2)>size(Img1Partemp,1)-borderGap(1), temp=gridx_g(2); 
            gridx_g(2)=size(Img1Partemp,1)-borderGap(1); gridx_f(2)=gridx_f(2)+((size(Img1Partemp,1)-borderGap(1))-temp); end
        if gridy_g(2)>size(Img1Partemp,2)-borderGap(2), temp=gridy_g(2); 
            gridy_g(2)=size(Img1Partemp,2)-borderGap(2); gridy_f(2)=gridy_f(2)+((size(Img1Partemp,2)-borderGap(2))-temp); end
        if gridz_g(2)>size(Img1Partemp,3)-borderGap(3), temp=gridz_g(2); 
            gridz_g(2)=size(Img1Partemp,3)-borderGap(3); gridz_f(2)=gridz_f(2)+((size(Img1Partemp,3)-borderGap(3))-temp); end
        
        % Make sure image width/length odd number
        if mod((gridx_f(2)-gridx_f(1)),2)==1, gridx_f(2)=gridx_f(2)-1; gridx_g(2)=gridx_g(2)-1; end
        if mod((gridy_f(2)-gridy_f(1)),2)==1, gridy_f(2)=gridy_f(2)-1; gridy_g(2)=gridy_g(2)-1; end
        if mod((gridz_f(2)-gridz_f(1)),2)==1, gridz_f(2)=gridz_f(2)-1; gridz_g(2)=gridz_g(2)-1; end
        
        gridxWidth0 = gridx_f(2)-gridx_f(1); gridyWidth0 = gridy_f(2)-gridy_f(1); gridzWidth0 = gridz_f(2)-gridz_f(1);
        gridx_f0 = gridx_f; gridy_f0 = gridy_f; gridz_f0=gridz_f;
        gridx_g0 = gridx_g; gridy_g0 = gridy_g; gridz_g0=gridz_g;
        
        % ====== Assign values ======
        gridxWidthCurr = gridxWidth0; gridyWidthCurr = gridyWidth0; gridzWidthCurr = gridzWidth0;
        gridxyRatioCurr=gridxWidthCurr/gridyWidthCurr;
        gridxzRatioCurr=gridxWidthCurr/gridzWidthCurr;
        gridyzRatioCurr=gridyWidthCurr/gridzWidthCurr;
        gridx_fCurr = gridx_f0; gridy_fCurr = gridy_f0; gridz_fCurr = gridz_f0;
        gridx_gCurr = gridx_g0; gridy_gCurr = gridy_g0; gridz_gCurr = gridz_g0;
        utempCurr = ceil(utemp); vtempCurr = ceil(vtemp); wtempCurr = ceil(wtemp);
        
        
        %% Level > 1
        levelNoTotal=1; clear gridx_fNew gridx_gNew gridy_fNew gridy_gNew gridz_fNew gridz_gNew utempNew vtempNew wtempNew qfactors
        TotalNo=1; gridxWidthNewtemp=gridxWidthCurr; gridyWidthNewtemp=gridyWidthCurr; gridzWidthNewtemp=gridzWidthCurr;
        xyRatio=gridxWidthNewtemp/gridyWidthNewtemp; xzRatio=gridxWidthNewtemp/gridzWidthNewtemp; yzRatio=gridyWidthNewtemp/gridzWidthNewtemp;
        while gridxWidthNewtemp/gridyWidthNewtemp > 0
            levelNoTotal=levelNoTotal+1;
            if xyRatio>2 && xzRatio>2    % split gridx only
                gridxWidthNewtemp=gridxWidthNewtemp/2; gridyWidthNewtemp=gridyWidthNewtemp; gridzWidthNewtemp=gridzWidthNewtemp;
                TotalNo = TotalNo*2;
            elseif xyRatio<0.5 && yzRatio>2 % split gridy only
                gridxWidthNewtemp=gridxWidthNewtemp; gridyWidthNewtemp=gridyWidthNewtemp/2; gridzWidthNewtemp=gridzWidthNewtemp;
                TotalNo = TotalNo*2;
            elseif xzRatio<0.5 && yzRatio<0.5 % split gridz only
                gridxWidthNewtemp=gridxWidthNewtemp; gridyWidthNewtemp=gridyWidthNewtemp; gridzWidthNewtemp=gridzWidthNewtemp/2;
            elseif xzRatio>2 && yzRatio>2 % split gridx & y
                gridxWidthNewtemp=gridxWidthNewtemp/2; gridyWidthNewtemp=gridyWidthNewtemp/2; gridzWidthNewtemp=gridzWidthNewtemp;
                TotalNo = TotalNo*4;
            elseif xyRatio>2 && yzRatio<0.5 % split gridx & z
                gridxWidthNewtemp=gridxWidthNewtemp/2; gridyWidthNewtemp=gridyWidthNewtemp; gridzWidthNewtemp=gridzWidthNewtemp/2;
                TotalNo = TotalNo*4;
            elseif xyRatio<0.5 && xzRatio<0.5 % split gridy & z
                gridxWidthNewtemp=gridxWidthNewtemp; gridyWidthNewtemp=gridyWidthNewtemp/2; gridzWidthNewtemp=gridzWidthNewtemp/2;
                TotalNo = TotalNo*4;
            else
                gridxWidthNewtemp=gridxWidthNewtemp/2; gridyWidthNewtemp=gridyWidthNewtemp/2; gridzWidthNewtemp=gridzWidthNewtemp/2;
                TotalNo = TotalNo*8;
            end
            xyRatio=gridxWidthNewtemp/gridyWidthNewtemp; xzRatio=gridxWidthNewtemp/gridzWidthNewtemp; yzRatio=gridyWidthNewtemp/gridzWidthNewtemp;
            if ( gridxWidthNewtemp<max([1*winsize(1),MinSearchWS]) || gridyWidthNewtemp<max([1*winsize(2),MinSearchWS]) || gridzWidthNewtemp<max([1*winsize(3),MinSearchWS]) )
                break
            end
        end
        
        
        IterNo=0; levelNo=1; hbar = waitbar(0,'FFT initial guess for large deformations.');
        while 1>0
            levelNo=levelNo+1;
            clear utemp vtemp Phitemp
            if gridxyRatioCurr>2 && gridxzRatioCurr>2    % split gridx only
                gridxWidthNew = gridxWidthCurr/2; gridyWidthNew = gridyWidthCurr; gridzWidthNew = gridzWidthCurr;
                for tempi = 1:size(gridx_fCurr,1)
                    tempj=tempi; tempInd=tempj;
                    utempNew(2*tempInd-1:2*tempInd) = repmat(utempCurr(tempInd),2,1);
                    vtempNew(2*tempInd-1:2*tempInd) = repmat(vtempCurr(tempInd),2,1);
                    wtempNew(2*tempInd-1:2*tempInd) = repmat(wtempCurr(tempInd),2,1);
                    gridx_fNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridx_fCurr(tempInd,1:2), 2,1);
                    gridx_gNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridx_gCurr(tempInd,1:2), 2,1);
                    gridy_fNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridy_fCurr(tempInd,1:2), 2,1);
                    gridy_gNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridy_gCurr(tempInd,1:2), 2,1);
                    gridz_fNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridz_fCurr(tempInd,1:2), 2,1);
                    gridz_gNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridz_gCurr(tempInd,1:2), 2,1);
                    % if ( (gridx_fCurr(tempInd,2)-gridx_fCurr(tempInd,1) > winsize(1))   )
                        gridx_fNew(2*tempInd-1:2*tempInd,1:2) = [ gridx_fCurr(tempInd,1), 0.5*sum(gridx_fCurr(tempInd,:));
                            0.5*sum(gridx_fCurr(tempInd,:)), gridx_fCurr(tempInd,2)];
                        gridx_gNew(2*tempInd-1:2*tempInd,1:2) = [ gridx_gCurr(tempInd,1), 0.5*sum(gridx_gCurr(tempInd,:));
                            0.5*sum(gridx_gCurr(tempInd,:)), gridx_gCurr(tempInd,2)];
                    % end
                end
                
            elseif gridxyRatioCurr<0.5 && gridyzRatioCurr>2 % split gridy only
                gridxWidthNew = gridxWidthCurr; gridyWidthNew = gridyWidthCurr/2; gridzWidthNew = gridzWidthCurr;
                for tempi = 1:size(gridx_fCurr,1)
                    tempj=tempi; tempInd = tempj;
                    utempNew(2*tempInd-1:2*tempInd) = repmat(utempCurr(tempInd),2,1);
                    vtempNew(2*tempInd-1:2*tempInd) = repmat(vtempCurr(tempInd),2,1);
                    wtempNew(2*tempInd-1:2*tempInd) = repmat(wtempCurr(tempInd),2,1);
                    gridx_fNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridx_fCurr(tempInd,1:2), 2,1);
                    gridx_gNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridx_gCurr(tempInd,1:2), 2,1);
                    gridy_fNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridy_fCurr(tempInd,1:2), 2,1);
                    gridy_gNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridy_gCurr(tempInd,1:2), 2,1);
                    gridz_fNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridz_fCurr(tempInd,1:2), 2,1);
                    gridz_gNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridz_gCurr(tempInd,1:2), 2,1);
                    % if (  (gridy_fCurr(tempInd,2)-gridy_fCurr(tempInd,1) > winsize(2))  )
                        gridy_fNew(2*tempInd-1:2*tempInd,1:2) = [ gridy_fCurr(tempInd,1), 0.5*sum(gridy_fCurr(tempInd,:));
                            0.5*sum(gridy_fCurr(tempInd,:)), gridy_fCurr(tempInd,2)];
                        gridy_gNew(2*tempInd-1:2*tempInd,1:2) = [ gridy_gCurr(tempInd,1), 0.5*sum(gridy_gCurr(tempInd,:));
                            0.5*sum(gridy_gCurr(tempInd,:)), gridy_gCurr(tempInd,2)];
                    % end
                end
                
                
            elseif gridxzRatioCurr<0.5 && gridyzRatioCurr<0.5 % split gridz only
                gridxWidthNew = gridxWidthCurr; gridyWidthNew = gridyWidthCurr; gridzWidthNew = gridzWidthCurr/2;
                for tempi = 1:size(gridx_fCurr,1)
                    tempj=tempi; tempInd = tempj;
                    utempNew(2*tempInd-1:2*tempInd) = repmat(utempCurr(tempInd),2,1);
                    vtempNew(2*tempInd-1:2*tempInd) = repmat(vtempCurr(tempInd),2,1);
                    wtempNew(2*tempInd-1:2*tempInd) = repmat(wtempCurr(tempInd),2,1);
                    gridx_fNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridx_fCurr(tempInd,1:2), 2,1);
                    gridx_gNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridx_gCurr(tempInd,1:2), 2,1);
                    gridy_fNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridy_fCurr(tempInd,1:2), 2,1);
                    gridy_gNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridy_gCurr(tempInd,1:2), 2,1);
                    gridz_fNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridz_fCurr(tempInd,1:2), 2,1);
                    gridz_gNew(2*tempInd-1:2*tempInd,1:2) = repmat( gridz_gCurr(tempInd,1:2), 2,1);
                    % if (  (gridz_fCurr(tempInd,2)-gridz_fCurr(tempInd,1) > winsize(3)) )
                        gridz_fNew(2*tempInd-1:2*tempInd,1:2) = [ gridz_fCurr(tempInd,1), 0.5*sum(gridz_fCurr(tempInd,:));
                            0.5*sum(gridz_fCurr(tempInd,:)), gridz_fCurr(tempInd,2)];
                        gridz_gNew(2*tempInd-1:2*tempInd,1:2) = [ gridz_gCurr(tempInd,1), 0.5*sum(gridz_gCurr(tempInd,:));
                            0.5*sum(gridz_gCurr(tempInd,:)), gridz_gCurr(tempInd,2)];
                    % end
                end
                
            elseif gridxzRatioCurr>2 && gridyzRatioCurr>2 % split gridx & y
                gridxWidthNew = gridxWidthCurr/2; gridyWidthNew = gridyWidthCurr/2; gridzWidthNew = gridzWidthCurr;
                for tempi = 1:size(gridx_fCurr,1)
                    tempj=tempi; tempInd = tempj;
                    utempNew(4*tempInd-3:4*tempInd) = repmat(utempCurr(tempInd),4,1);
                    vtempNew(4*tempInd-3:4*tempInd) = repmat(vtempCurr(tempInd),4,1);
                    wtempNew(4*tempInd-3:4*tempInd) = repmat(wtempCurr(tempInd),4,1);
                    gridx_fNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridx_fCurr(tempInd,1:2), 4,1);
                    gridx_gNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridx_gCurr(tempInd,1:2), 4,1);
                    gridy_fNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridy_fCurr(tempInd,1:2), 4,1);
                    gridy_gNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridy_gCurr(tempInd,1:2), 4,1);
                    gridz_fNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridz_fCurr(tempInd,1:2), 4,1);
                    gridz_gNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridz_gCurr(tempInd,1:2), 4,1);
                    % if ( (gridx_fCurr(tempInd,2)-gridx_fCurr(tempInd,1) > winsize(1)) && (gridy_fCurr(tempInd,2)-gridy_fCurr(tempInd,1) > winsize(2))  )
                        gridx_fNew(4*tempInd-3:4*tempInd,1:2) = repmat([ gridx_fCurr(tempInd,1), 0.5*sum(gridx_fCurr(tempInd,:));
                            0.5*sum(gridx_fCurr(tempInd,:)), gridx_fCurr(tempInd,2)], 2,1);
                        gridx_gNew(4*tempInd-3:4*tempInd,1:2) = repmat([ gridx_gCurr(tempInd,1), 0.5*sum(gridx_gCurr(tempInd,:));
                            0.5*sum(gridx_gCurr(tempInd,:)), gridx_gCurr(tempInd,2)], 2,1);
                        gridy_fNew(4*tempInd-3:4*tempInd,1:2) = [ repmat([ gridy_fCurr(tempInd,1), 0.5*sum(gridy_fCurr(tempInd,:))],2,1);
                            repmat([0.5*sum(gridy_fCurr(tempInd,:)), gridy_fCurr(tempInd,2)],2,1) ];
                        gridy_gNew(4*tempInd-3:4*tempInd,1:2) =  [ repmat([ gridy_gCurr(tempInd,1), 0.5*sum(gridy_gCurr(tempInd,:))],2,1);
                            repmat([0.5*sum(gridy_gCurr(tempInd,:)), gridy_gCurr(tempInd,2)],2,1) ];
                    % end
                end
                
            elseif gridxyRatioCurr>2 && gridyzRatioCurr<0.5 % split gridx & z
                gridxWidthNew = gridxWidthCurr/2; gridyWidthNew = gridyWidthCurr; gridzWidthNew = gridzWidthCurr/2;
                for tempi = 1:size(gridx_fCurr,1)
                    tempj=tempi; tempInd = tempj;
                    utempNew(4*tempInd-3:4*tempInd) = repmat(utempCurr(tempInd),4,1);
                    vtempNew(4*tempInd-3:4*tempInd) = repmat(vtempCurr(tempInd),4,1);
                    wtempNew(4*tempInd-3:4*tempInd) = repmat(wtempCurr(tempInd),4,1);
                    gridx_fNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridx_fCurr(tempInd,1:2), 4,1);
                    gridx_gNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridx_gCurr(tempInd,1:2), 4,1);
                    gridy_fNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridy_fCurr(tempInd,1:2), 4,1);
                    gridy_gNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridy_gCurr(tempInd,1:2), 4,1);
                    gridz_fNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridz_fCurr(tempInd,1:2), 4,1);
                    gridz_gNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridz_gCurr(tempInd,1:2), 4,1);
                    % if ( (gridx_fCurr(tempInd,2)-gridx_fCurr(tempInd,1) > winsize(1)) &&  (gridz_fCurr(tempInd,2)-gridz_fCurr(tempInd,1) > winsize(3)) )
                        gridx_fNew(4*tempInd-3:4*tempInd,1:2) = repmat([ gridx_fCurr(tempInd,1), 0.5*sum(gridx_fCurr(tempInd,:));
                            0.5*sum(gridx_fCurr(tempInd,:)), gridx_fCurr(tempInd,2)], 2,1);
                        gridx_gNew(4*tempInd-3:4*tempInd,1:2) = repmat([ gridx_gCurr(tempInd,1), 0.5*sum(gridx_gCurr(tempInd,:));
                            0.5*sum(gridx_gCurr(tempInd,:)), gridx_gCurr(tempInd,2)], 2,1);
                        gridz_fNew(4*tempInd-3:4*tempInd,1:2) = [ repmat([ gridz_fCurr(tempInd,1), 0.5*sum(gridz_fCurr(tempInd,:))],2,1);
                            repmat([0.5*sum(gridz_fCurr(tempInd,:)), gridz_fCurr(tempInd,2)],2,1) ];
                        gridz_gNew(4*tempInd-3:4*tempInd,1:2) =  [ repmat([ gridz_gCurr(tempInd,1), 0.5*sum(gridz_gCurr(tempInd,:))],2,1);
                            repmat([0.5*sum(gridz_gCurr(tempInd,:)), gridz_gCurr(tempInd,2)],2,1) ];
                    % end
                end
                
                
                
            elseif gridxyRatioCurr<0.5 && gridxzRatioCurr<0.5 % split gridy & z
                gridxWidthNew = gridxWidthCurr; gridyWidthNew = gridyWidthCurr/2; gridzWidthNew = gridzWidthCurr/2;
                for tempi = 1:size(gridx_fCurr,1)
                    tempj=tempi; tempInd = tempj;
                    utempNew(4*tempInd-3:4*tempInd) = repmat(utempCurr(tempInd),4,1);
                    vtempNew(4*tempInd-3:4*tempInd) = repmat(vtempCurr(tempInd),4,1);
                    wtempNew(4*tempInd-3:4*tempInd) = repmat(wtempCurr(tempInd),4,1);
                    gridx_fNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridx_fCurr(tempInd,1:2), 4,1);
                    gridx_gNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridx_gCurr(tempInd,1:2), 4,1);
                    gridy_fNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridy_fCurr(tempInd,1:2), 4,1);
                    gridy_gNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridy_gCurr(tempInd,1:2), 4,1);
                    gridz_fNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridz_fCurr(tempInd,1:2), 4,1);
                    gridz_gNew(4*tempInd-3:4*tempInd,1:2) = repmat( gridz_gCurr(tempInd,1:2), 4,1);
                    % if (  (gridy_fCurr(tempInd,2)-gridy_fCurr(tempInd,1) > winsize(2)) && (gridz_fCurr(tempInd,2)-gridz_fCurr(tempInd,1) > winsize(3)) )
                        gridy_fNew(4*tempInd-3:4*tempInd,1:2) = repmat([ gridy_fCurr(tempInd,1), 0.5*sum(gridy_fCurr(tempInd,:));
                            0.5*sum(gridy_fCurr(tempInd,:)), gridy_fCurr(tempInd,2)], 2,1);
                        gridy_gNew(4*tempInd-3:4*tempInd,1:2) = repmat([ gridy_gCurr(tempInd,1), 0.5*sum(gridy_gCurr(tempInd,:));
                            0.5*sum(gridy_gCurr(tempInd,:)), gridy_gCurr(tempInd,2)], 2,1);
                        gridz_fNew(4*tempInd-3:4*tempInd,1:2) = [ repmat([ gridz_fCurr(tempInd,1), 0.5*sum(gridz_fCurr(tempInd,:))],2,1);
                            repmat([0.5*sum(gridz_fCurr(tempInd,:)), gridz_fCurr(tempInd,2)],2,1) ];
                        gridz_gNew(4*tempInd-3:4*tempInd,1:2) =  [ repmat([ gridz_gCurr(tempInd,1), 0.5*sum(gridz_gCurr(tempInd,:))],2,1);
                            repmat([0.5*sum(gridz_gCurr(tempInd,:)), gridz_gCurr(tempInd,2)],2,1) ];
                    % end
                end
                
                
            else % split gridx & y & z
                
                gridxWidthNew = gridxWidthCurr/2; gridyWidthNew = gridyWidthCurr/2; gridzWidthNew = gridzWidthCurr/2;
                for tempi = 1:size(gridx_fCurr,1)
                    tempj=tempi; tempInd = tempj;
                    utempNew(8*tempInd-7:8*tempInd) = repmat(utempCurr(tempInd),8,1);
                    vtempNew(8*tempInd-7:8*tempInd) = repmat(vtempCurr(tempInd),8,1);
                    wtempNew(8*tempInd-7:8*tempInd) = repmat(wtempCurr(tempInd),8,1);
                    gridx_fNew(8*tempInd-7:8*tempInd,1:2) = repmat( gridx_fCurr(tempInd,1:2), 8,1);
                    gridx_gNew(8*tempInd-7:8*tempInd,1:2) = repmat( gridx_gCurr(tempInd,1:2), 8,1);
                    gridy_fNew(8*tempInd-7:8*tempInd,1:2) = repmat( gridy_fCurr(tempInd,1:2), 8,1);
                    gridy_gNew(8*tempInd-7:8*tempInd,1:2) = repmat( gridy_gCurr(tempInd,1:2), 8,1);
                    gridz_fNew(8*tempInd-7:8*tempInd,1:2) = repmat( gridz_fCurr(tempInd,1:2), 8,1);
                    gridz_gNew(8*tempInd-7:8*tempInd,1:2) = repmat( gridz_gCurr(tempInd,1:2), 8,1);
                    % if ( (gridx_fCurr(tempInd,2)-gridx_fCurr(tempInd,1) > winsize(1)) && (gridy_fCurr(tempInd,2)-gridy_fCurr(tempInd,1) > winsize(2)) && (gridz_fCurr(tempInd,2)-gridz_fCurr(tempInd,1) > winsize(3)) )
                        gridx_fNew(8*tempInd-7:8*tempInd,1:2) = repmat([ gridx_fCurr(tempInd,1), 0.5*sum(gridx_fCurr(tempInd,:));
                            0.5*sum(gridx_fCurr(tempInd,:)), gridx_fCurr(tempInd,2)], 4,1);
                        gridx_gNew(8*tempInd-7:8*tempInd,1:2) = repmat([ gridx_gCurr(tempInd,1), 0.5*sum(gridx_gCurr(tempInd,:));
                            0.5*sum(gridx_gCurr(tempInd,:)), gridx_gCurr(tempInd,2)], 4,1);
                        gridy_fNew(8*tempInd-7:8*tempInd,1:2) = [ repmat([ gridy_fCurr(tempInd,1), 0.5*sum(gridy_fCurr(tempInd,:))], 2,1);
                            repmat([ 0.5*sum(gridy_fCurr(tempInd,:)), gridy_fCurr(tempInd,2)], 2,1) ;
                            repmat([ gridy_fCurr(tempInd,1), 0.5*sum(gridy_fCurr(tempInd,:))], 2,1);
                            repmat([ 0.5*sum(gridy_fCurr(tempInd,:)), gridy_fCurr(tempInd,2)], 2,1) ];
                        gridy_gNew(8*tempInd-7:8*tempInd,1:2) = [ repmat([ gridy_gCurr(tempInd,1), 0.5*sum(gridy_gCurr(tempInd,:))], 2,1);
                            repmat([ 0.5*sum(gridy_gCurr(tempInd,:)), gridy_gCurr(tempInd,2)], 2,1) ;
                            repmat([ gridy_gCurr(tempInd,1), 0.5*sum(gridy_gCurr(tempInd,:))], 2,1);
                            repmat([ 0.5*sum(gridy_gCurr(tempInd,:)), gridy_gCurr(tempInd,2)], 2,1) ] ;
                        gridz_fNew(8*tempInd-7:8*tempInd,1:2) = [ repmat([ gridz_fCurr(tempInd,1), 0.5*sum(gridz_fCurr(tempInd,:))],4,1);
                            repmat([ 0.5*sum(gridz_fCurr(tempInd,:)),  gridz_fCurr(tempInd,2)  ],4,1) ];
                        gridz_gNew(8*tempInd-7:8*tempInd,1:2) = [ repmat([ gridz_gCurr(tempInd,1), 0.5*sum(gridz_gCurr(tempInd,:))],4,1);
                            repmat([ 0.5*sum(gridz_gCurr(tempInd,:)),  gridz_gCurr(tempInd,2)  ],4,1) ];
                    % end
                end
            end
            
            
            for tempi = 1:size(gridx_fNew,1)
                IterNo = IterNo+1;
                waitbar((IterNo-1)/TotalNo);
                
                C = Img1Partemp(gridx_fNew(tempi,1):gridx_fNew(tempi,2), ...
                    gridy_fNew(tempi,1):gridy_fNew(tempi,2), ...
                    gridz_fNew(tempi,1):gridz_fNew(tempi,2) );
                if (size(C,1)*size(C,2)*size(C,3)>MinSearchWS^3)
                    D = Img2Partemp(gridx_gNew(tempi,1):gridx_gNew(tempi,2), ...
                        gridy_gNew(tempi,1):gridy_gNew(tempi,2), ...
                        gridz_gNew(tempi,1):gridz_gNew(tempi,2) );
                else
                    D = Img2Partemp(gridx_gNew(tempi,1)-5:gridx_gNew(tempi,2)+5, ...
                    gridy_gNew(tempi,1)-5:gridy_gNew(tempi,2)+5, ...
                    gridz_gNew(tempi,1)-5:gridz_gNew(tempi,2)+5 );
                end
                
                RDTime = 1;
                while size(C,1)*size(C,2)*size(C,3)>MaxVoxelNum
                    C = imgaussfilt3(C); C = imgaussfilt3(C); C = imgaussfilt3(C);
                    D = imgaussfilt3(D); D = imgaussfilt3(D); D = imgaussfilt3(D);
                    RDTime = RDTime*2;
                    CNew = 1/8*( C(1:2:2*floor(size(C,1)/2),1:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                        C(2:2:2*floor(size(C,1)/2), 1:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                        C(1:2:2*floor(size(C,1)/2),2:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                        C(2:2:2*floor(size(C,1)/2), 2:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                        C(1:2:2*floor(size(C,1)/2),1:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)) + ...
                        C(2:2:2*floor(size(C,1)/2), 1:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)) + ...
                        C(1:2:2*floor(size(C,1)/2),2:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)) + ...
                        C(2:2:2*floor(size(C,1)/2), 2:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)) ) ;
                    DNew = 1/8*( D(1:2:2*floor(size(C,1)/2),1:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                        D(2:2:2*floor(size(C,1)/2), 1:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                        D(1:2:2*floor(size(C,1)/2),2:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                        D(2:2:2*floor(size(C,1)/2), 2:2:2*floor(size(C,2)/2),1:2:2*floor(size(C,3)/2)) + ...
                        D(1:2:2*floor(size(C,1)/2),1:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)) + ...
                        D(2:2:2*floor(size(C,1)/2), 1:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)) + ...
                        D(1:2:2*floor(size(C,1)/2),2:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)) + ...
                        D(2:2:2*floor(size(C,1)/2), 2:2:2*floor(size(C,2)/2),2:2:2*floor(size(C,3)/2)) ); 
                    C=CNew; D=DNew;
                end
                % figure, imagesc3D(C); figure, imagesc3D(D);
                % figure, imagesc3D(Img1Partemp); figure, imagesc3D(Img2Partemp);
          
                % cross-correlation
                if strcmp(fftmethod,'bigxcorr') 
                     XCORRF3OfCD0 = normxcorr3(C,D,'full');
                elseif strcmp(fftmethod,'bigphasecorr')
                    XCORRF3OfCD0 = bigPhaseCorr3(C,D,'full');
                end
                % figure, imagesc3D(XCORRF3OfCD0);
                
                % find qfactors
                cc.A{1} = real(XCORRF3OfCD0);
                qfactors(tempi,:) = compute_qFactor(cc,tempi);
                if RDTime==1
                    winsize1 = size(C,1)-1; winsize2 = size(C,2)-1; winsize3 = size(C,3)-1;
                    [v1temp,u1temp,w1temp,max_f] = findpeak3(real(XCORRF3OfCD0(winsize1:end-winsize1+1, ...
                        winsize2:end-winsize2+1, winsize3:end-winsize3+1)),1);
                    zero_disp = ceil((size(XCORRF3OfCD0(winsize1:end-winsize1+1,winsize2:end-winsize2+1,winsize3:end-winsize3+1)))/2);
                else
                    [v1temp, u1temp, w1temp, max_f] = findpeak3(real(XCORRF3OfCD0),1);
                    zero_disp = ceil((size(XCORRF3OfCD0)+[1,1,1])/2);
                end
                ind = tempi;
                
                utempNew(ind) = utempNew(ind) + RDTime*(u1temp-zero_disp(1));
                vtempNew(ind) = vtempNew(ind) + RDTime*(v1temp-zero_disp(2));
                wtempNew(ind) = wtempNew(ind) + RDTime*(w1temp-zero_disp(3));
                Phitemp(ind) = max_f;
                
            end
            
            if ( gridxWidthNew<max([1*winsize(1),MinSearchWS]) || gridyWidthNew<max([1*winsize(2),MinSearchWS]) || gridzWidthNew<max([1*winsize(3),MinSearchWS]) )
                % Finish qfactor computation
                for k=1:2
                    qf_ = (qfactors(:,k)-min(qfactors(:,k)));
                    cc.qfactors(:,k) = qf_/max(qf_);
                end
                break
            else
                %% ====== Check out of image border or not ======
                gridx_f=gridx_fNew; gridy_f=gridy_fNew; gridz_f=gridz_fNew;
                gridx_g=gridx_fNew + ceil(utempNew)'*[1,1];
                gridy_g=gridy_fNew + ceil(vtempNew)'*[1,1];
                gridz_g=gridz_fNew + ceil(wtempNew)'*[1,1];
                for tempi = 1:size(gridx_g,1)
                    if gridx_g(tempi,1)<borderGap(1)
                        temp=gridx_g(tempi,1);
                        gridx_g(tempi,1)=borderGap(1);
                        gridx_f(tempi,1) = gridx_f(tempi,1)+borderGap(1)-temp;
                    end
                    if gridy_g(tempi,1)<borderGap(2)
                        temp=gridy_g(tempi,1);
                        gridy_g(tempi,1)=borderGap(2);
                        gridy_f(tempi,1) = gridy_f(tempi,1)+borderGap(2)-temp;
                    end
                    if gridz_g(tempi,1)<borderGap(3)
                        temp=gridz_g(tempi,1);
                        gridz_g(tempi,1)=borderGap(3);
                        gridz_f(tempi,1) = gridz_f(tempi,1)+borderGap(3)-temp;
                    end
                    if gridx_g(tempi,2)>size(Img1Partemp,1)-borderGap(1)
                        temp=gridx_g(tempi,2);
                        gridx_g(tempi,2)=size(Img1Partemp,1)-borderGap(1);
                        gridx_f(tempi,2) = gridx_f(tempi,2)+(size(Img1Partemp,1)-borderGap(1))-temp;
                    end
                    if gridy_g(tempi,2)>size(Img1Partemp,2)-borderGap(2)
                        temp=gridy_g(tempi,2);
                        gridy_g(tempi,2)=size(Img1Partemp,2)-borderGap(2);
                        gridy_f(tempi,2) = gridy_f(tempi,2)+(size(Img1Partemp,2)-borderGap(2))-temp;
                    end
                    if gridz_g(tempi,2)>size(Img1Partemp,3)-borderGap(3)
                        temp=gridz_g(tempi,2);
                        gridz_g(tempi,2)=size(Img1Partemp,3)-borderGap(3);
                        gridz_f(tempi,2) = gridz_f(tempi,2)+(size(Img1Partemp,3)-borderGap(3))-temp;
                    end
                end
                % Make sure: image width/length odd number
                for tempi = 1:size(gridx_g,1)
                    if mod(gridx_f(tempi,2)-gridx_f(tempi,1),2)==1
                        gridx_f(tempi,2)=gridx_f(tempi,2)-1;
                        gridx_g(tempi,2)=gridx_g(tempi,2)-1;
                    end
                    if mod(gridy_f(tempi,2)-gridy_f(tempi,1),2)==1
                        gridy_f(tempi,2)=gridy_f(tempi,2)-1;
                        gridy_g(tempi,2)=gridy_g(tempi,2)-1;
                    end
                    if mod(gridz_f(tempi,2)-gridz_f(tempi,1),2)==1
                        gridz_f(tempi,2)=gridz_f(tempi,2)-1;
                        gridz_g(tempi,2)=gridz_g(tempi,2)-1;
                    end
                end
                
                %% ====== Assign values ======
                gridxWidthCurr = gridxWidthNew; gridyWidthCurr = gridyWidthNew; gridzWidthCurr = gridzWidthNew;
                gridxyRatioCurr=gridxWidthCurr/gridyWidthCurr;
                gridxzRatioCurr=gridxWidthCurr/gridzWidthCurr;
                gridyzRatioCurr=gridyWidthCurr/gridzWidthCurr;
                gridx_fCurr = gridx_f; gridy_fCurr = gridy_f; gridz_fCurr=gridz_f;
                gridx_gCurr = gridx_g; gridy_gCurr = gridy_g; gridz_gCurr=gridz_g;
                utempCurr =  ceil(utempNew); vtempCurr =  ceil(vtempNew); wtempCurr=ceil(wtempNew);
                
            end
            
        end
        
        close(hbar);
        %%
        tempx=double(0.5*(gridx_fNew(:,1)+gridx_fNew(:,2))); tempx=tempx(:);
        tempy=double(0.5*(gridy_fNew(:,1)+gridy_fNew(:,2))); tempy=tempy(:);
        tempz=0.5*(gridz_fNew(:,1)+gridz_fNew(:,2)); tempz=tempz(:);
        tempu=double(utempNew); tempu=tempu(:);
        tempv=double(vtempNew); tempv=tempv(:);
        tempw=double(wtempNew); tempw=tempw(:);
        tempPhi=double(Phitemp); tempPhi=tempPhi(:);
        qf_1 = double(cc.qfactors(:,1)'); qf_2 = double(cc.qfactors(:,2)');
        qf_1 = qf_1(:); qf_2 = qf_2(:);
        
        %% Update image domain
        [~,indx1]=min(tempx); [~,indx2]=max(tempx);
        [~,indy1]=min(tempy); [~,indy2]=max(tempy);
        [~,indz1]=min(tempz); [~,indz2]=max(tempz);
        borderGapx1=ceil(1+(1.4*winsize(1))+((tempu(indx1)))); borderGapx2=ceil(1+(1.4*winsize(1))+((tempu(indx2)))); 
        borderGapy1=ceil(1+(1.4*winsize(2))+((tempv(indy1)))); borderGapy2=ceil(1+(1.4*winsize(2))+((tempv(indy2))));
        borderGapz1=ceil(1+(1.4*winsize(3))+((tempw(indz1)))); borderGapz2=ceil(1+(1.4*winsize(3))+((tempw(indz2))));
        
        if gridx(1) < borderGapx1, gridx(1) = borderGapx1; end
        if gridy(1) < borderGapy1, gridy(1) = borderGapy1; end
        if gridz(1) < borderGapz1, gridz(1) = borderGapz1; end
        if gridx(2) > size(Img1Partemp,1)-borderGapx2+1, gridx(2) = size(Img1Partemp,1)-borderGapx2+1; end
        if gridy(2) > size(Img1Partemp,2)-borderGapy2+1, gridy(2) = size(Img1Partemp,2)-borderGapy2+1; end
        if gridz(2) > size(Img1Partemp,3)-borderGapz2+1, gridz(2) = size(Img1Partemp,3)-borderGapz2+1; end
  
        %% Interpolate
        % figure, plot3(tempx,tempy,tempz,'.')
        % figure, scatter3(tempx,tempy,tempz,1*ones(length(tempu),1),tempw)
        % figure, quiver3(tempx,tempy,tempz,tempu,tempv,tempw)
        % figure, quiver3(xGrid(:),yGrid(:),zGrid(:),uGrid(:),vGrid(:),wGrid(:))
        % figure, quiver3(xGrid(:),yGrid(:),zGrid(:),uGrid2(:),vGrid2(:),wGrid2(:))
        % figure, surf(squeeze(xGrid(:,:,end)),squeeze(yGrid(:,:,end)),squeeze(vGrid2(:,:,end)),'edgecolor','none');

        try xList=[gridx(1):winstepsize(1):gridx(2)]; catch, xList=[gridx(1):winstepsize:gridx(2)]; end
        try yList=[gridy(1):winstepsize(2):gridy(2)]; catch, yList=[gridy(1):winstepsize:gridy(2)]; end
        try zList=[gridz(1):winstepsize(3):gridz(2)]; catch, zList=[gridz(1):winstepsize:gridz(2)]; end
        [indx]=find(tempx>min(xList) & tempx<max(xList));
        [indy]=find(tempy>min(yList) & tempy<max(yList));
        [indz]=find(tempz>min(zList) & tempz<max(zList));
        [indxyz] = intersect(intersect(indx,indy),indz);
        hbar = waitbar(0,'FFT initial guess: assign values.');
        
        gridPoints = {xList, yList, zList}; smoothness=1e-3;
        uGrid2 = regularizeNd([tempx(indxyz),tempy(indxyz),tempz(indxyz)],tempu(indxyz),gridPoints, smoothness);
        
        [xGrid,yGrid,zGrid]=ndgrid(xList,yList,zList);
        F=scatteredInterpolant(tempx(indxyz),tempy(indxyz),tempz(indxyz),tempu(indxyz),'natural','nearest');
        uGrid=F(xGrid,yGrid,zGrid);
        %uGrid=griddata(tempx(indxyz),tempy(indxyz),tempz(indxyz),tempu(indxyz),xGrid,yGrid,zGrid,'natural');
        waitbar(1/6);
        vGrid2 = regularizeNd([tempx(indxyz),tempy(indxyz),tempz(indxyz)],tempv(indxyz),gridPoints, smoothness);
        
        F=scatteredInterpolant(tempx(indxyz),tempy(indxyz),tempz(indxyz),tempv(indxyz),'natural','nearest');
        vGrid=F(xGrid,yGrid,zGrid);
        %vGrid=griddata(tempx(indxyz),tempy(indxyz),tempz(indxyz),tempv(indxyz),xGrid,yGrid,zGrid,'natural');
        waitbar(2/6);
        wGrid2 = regularizeNd([tempx(indxyz),tempy(indxyz),tempz(indxyz)],tempw(indxyz),gridPoints, smoothness);
        
        F=scatteredInterpolant(tempx(indxyz),tempy(indxyz),tempz(indxyz),tempw(indxyz),'natural','nearest');
        wGrid=F(xGrid,yGrid,zGrid);
        %wGrid=griddata(tempx(indxyz),tempy(indxyz),tempz(indxyz),tempw(indxyz),xGrid,yGrid,zGrid,'natural');
        waitbar(3/6);
        F=scatteredInterpolant(tempx(indxyz),tempy(indxyz),tempz(indxyz),tempPhi(indxyz),'natural','nearest');
        PhiGrid=F(xGrid,yGrid,zGrid);
        %PhiGrid=griddata(tempx(indxyz),tempy(indxyz),tempz(indxyz),tempPhi(indxyz),xGrid,yGrid,zGrid,'natural');
        waitbar(4/6);
        F=scatteredInterpolant(tempx(indxyz),tempy(indxyz),tempz(indxyz),qf_1(indxyz),'natural','nearest');
        qf_1Grid=F(xGrid,yGrid,zGrid);
        %[qf_1Grid]=griddata(tempx(indxyz),tempy(indxyz),tempz(indxyz),qf_1(indxyz),xGrid,yGrid,zGrid,'natural');
        waitbar(5/6);
        F=scatteredInterpolant(tempx(indxyz),tempy(indxyz),tempz(indxyz),qf_2(indxyz),'natural','nearest');
        qf_2Grid=F(xGrid,yGrid,zGrid);
        %[qf_2Grid]=griddata(tempx(indxyz),tempy(indxyz),tempz(indxyz),qf_2(indxyz),xGrid,yGrid,zGrid,'natural');
        waitbar(6/6);
        
        cc.max = PhiGrid; cc.A = []; cc.qfactors =[qf_1Grid(:),qf_2Grid(:)];
        
        xyz.x = xGrid; xyz.y = yGrid; xyz.z = zGrid;
        uvw.u = uGrid; uvw.v = vGrid; uvw.w = wGrid;
        cc.max = PhiGrid;
        cc.A = [];
        
        close(hbar);
        
        %% --------------------------------------------------------------------------
    otherwise
        disp('not finished in funIntegerSearch3.m')
end
% -------- End of Local integer search --------

% -------- Part for Freq integer search --------
% D = g(ii:ii+winsize, jj:jj+winsize);
% D = g(ii-tempSizeOfSearchRegion:ii+winsize+tempSizeOfSearchRegion, ...
%    jj-tempSizeOfSearchRegion:jj+winsize+tempSizeOfSearchRegion);
% tformEstimate0 = imregcorr2(D,C,'purestretch');
% u(cj1,ci1) = tformEstimate0.T(3,2)+tempSizeOfSearchRegion;
% v(cj1,ci1) = tformEstimate0.T(3,1)+tempSizeOfSearchRegion;
% u(cj1,ci1) = tformEstimate0.T(2,2); % But is for stretch
% v(cj1,ci1) = tformEstimate0.T(1,1); % But is for stretch
% Phi(cj1,ci1) = 1; % It's meaningless
% -------- End of Part for Freq integer search --------

% -------- DELETE part --------
% Without interpolating for peak value if correlation function max_f is almost 1.
% I tried following code, unfortunately it doesn't work very well
% [v1temp1, u1temp1, max_f1] = findpeak(XCORRF2OfCD0(winsize:end-winsize+1,winsize:end-winsize+1),1);
% [v1temp2, u1temp2, max_f2] = findpeak(XCORRF2OfCD0(winsize:end-winsize+1,winsize:end-winsize+1),0);
%
% if max_f2 > 0.999
%    v1temp = v1temp2; u1temp = u1temp2; max_f = max_f2;
% else
%    v1temp = v1temp1; u1temp = u1temp1; max_f = max_f1;
% end
% -------- End of DELETE part --------

disp('Finish initial guess search!');
% -------- End of Local integer search --------
 
end

%% =============================================
function A = xCorr3(A,B,sSize)
% performs fft based cross correlation of A and B (see equation 2)

A = fftn(A,sSize);
B = fftn(B,sSize);
B = conj(B);
A = A.*B;
A = ifftn(A);
A = real(A);
A = fftshift(A);
end

%% =============================================
function A = phaseCorr3(A,B,sSize)
% performs fft based cross correlation of A and B (see equation 2)

xFreq = -0.5*(sSize(1)-1):1:0.5*(sSize(1)-1);
yFreq = -0.5*(sSize(2)-1):1:0.5*(sSize(2)-1);
zFreq = -0.5*(sSize(3)-1):1:0.5*(sSize(3)-1);
[XXFreq,YYFreq,ZZFreq] = ndgrid(xFreq,yFreq,zFreq);
kSq = (XXFreq.^2 + YYFreq.^2 + ZZFreq.^2);
kSq = ifftshift(kSq);

A = fftn(A,sSize);
B = fftn(B,sSize);
B = conj(B);
A = kSq.*A.*B;
A = ifftn(A);
A = real(A);
A = fftshift(A);

end

%% ========================================================================
function    duvw = lsqPolyFit3(b, M, trMM)
% LeastSqPoly performs a 3D polynomial fit in the least squares sense
% Solves M*x = b,
% trMM = transpose(M)*M
% trMb = tranpose(M)*b
%
% If you need to generate the coefficients then uncomment the following
% [mx, my, mz] = meshgrid(-1:1,-1:1,-1:1);
% m = [mx(:), my(:), mz(:)];
%
% for i = 1:size(m,1)
%    x = m(i,1); y = m(i,2); z = m(i,3);
%    M1(i,:) = [1,x,y,z,x^2,x*y,x*z,y^2,y*z,z^2];
% end
%
% trMM1 = M'*M;

% b = log(b);
trMb = sum(bsxfun(@times, M, b));

x = trMM\trMb'; %solve for unknown coefficients

A = [x(6), 2*x(5), x(7);
    2*x(8),  x(6) x(9)
    x(9),    x(7), 2*x(10)];

duvw = (A\(-x([2 3 4])));
end



%% ========================================================================
function varargout = generateMTF(sSize)
% MTF functions taken from
% J. Nogueira, A Lecuona, P. A. Rodriguez, J. A. Alfaro, and A. Acosta.
% Limits on the resolution of correlation PIV iterative methods. Practical
% implementation and design of weighting functions. Exp. Fluids,
% 39(2):314{321, July 2005. doi: 10.1007/s00348-005-1017-1

%% equation 4

if prod(single(sSize == 32)) || prod(single(sSize == 64)) || prod(single(sSize == 16))
    sSize = sSize(1);
    
    x = cell(1,3);
    [x{1}, x{2}, x{3}] = meshgrid(1:sSize,1:sSize,1:sSize);
    
    nu{1} = 1;
    for i = 1:3
        x{i} = x{i} - sSize/2 - 0.5;
        x{i} = abs(x{i}/sSize);
        nu{1} = nu{1}.*(3*(4*x{i}.^2-4*x{i}+1));
    end
    
    %% equation 5
    [x{1}, x{2}, x{3}] = meshgrid(1:sSize,1:sSize,1:sSize);
    
    for i = 1:3, x{i} = x{i} - sSize/2 - 0.5; end
    
    r = abs(sqrt(x{1}.^2 + x{2}.^2 + x{3}.^2)/sSize);
    nu{2}  = zeros(size(r));
    nu{2}(r < 0.5) = 24/pi*(4*r(r < 0.5).^2-4*r(r < 0.5)+1);
    
    %% equation 6
    [x{1}, x{2}, x{3}] = meshgrid(1:sSize,1:sSize,1:sSize);
    
    nu{3} = 1;
    for i = 1:3
        x{i} = x{i} - sSize/2 - 0.5;
        x{i} = (x{i}/sSize);
        
        nu{3} = nu{3}.*(12*abs(x{i}).^2 - 12*abs(x{i}) + 3 + ...
            0.15*cos(4*pi*x{i}) + 0.20*cos(6*pi*x{i}) + ...
            0.10*cos(8*pi*x{i}) + 0.05*cos(10*pi*x{i}));
        
    end
    nu{3}(nu{3} < 0) = 0;
    
else
    
    nu{1} = ones(sSize(1),sSize(2),sSize(3));
    nu{2} = nu{1};
    nu{3} = nu{1};
    
end

nu = cellfun(@(x) x/sum(x(:)), nu, 'UniformOutput',0);
nu = cellfun(@sqrt, nu, 'UniformOutput',0);

varargout = nu;

end


%% ========================================================================
function C = normxcorr3(T, A, shape)
% C = normxcorr3(TEMPLATE, IMAGE, SHAPE)
%
%       TEMPLATE - type double, ndims==3, size <= size of image
%       IMAGE    - type double, ndims==3
%       SHAPE    - one of: 'valid', 'same', 'full'. same as conv2 shape parameter
%                  'full' by default
%
%       C        - values in [-1,1]. size depends on SHAPE
%
% the syntax of this function is identical to Matlab's
% normxcorr2, except that it's been extended to 3D matrices,
% and, the SHAPE parameter has been introduced as a convenience
%
% the SHAPE parameter has the same effect as it does for the CONVN function.
% see the documentation for CONVN for a more detailed explanation
%
% caveat emptor: this function does not perform the argument checking that
% normxcorr2 does. for example, it doesn't ensure that std(T(:))~=0
% -------------------------------------------------------------------------
% daniel eaton, 2005, danieljameseaton@gmail.com

if nargin<3, shape = 'full'; end

if ndims(A)~=3 || ndims(T)~=3, error('A and T must be 3 dimensional matrices'); end

szT = size(T); szA = size(A);

if any(szT>szA), error('template must be smaller than image'); end

pSzT = prod(szT);

% =============== Comment =================
% make the running-sum/integral-images of A and A^2, which are
% used to speed up the computation of the NCC denominator
% ---------------------------------------
intImgA = integralImage(A,szT);
intImgA2 = integralImage(A.*A,szT);
% ============ End of comment =============
szOut = [(size(T,1)+size(A,1)-1), (size(T,2)+size(A,2)-1), (size(T,3)+size(A,3)-1)]; % szOut = size(intImgA);

% compute the numerator of the NCC
% emulate 3D correlation by rotating templates dimensions in 3D
% requency-domain correlation is MUCH faster than the spatial-domain variety
rotT = flipdim(flipdim(flipdim(T,1),2),3); % this is rot90 in 3d
fftRotT = fftn(rotT,szOut);
fftA = fftn(A,szOut);

corrTA = real(ifftn(fftA.*fftRotT,'symmetric'));

num = (corrTA - intImgA*sum(T(:))/pSzT ) / (pSzT-1);

% compute the denominator of the NCC
denomA = sqrt( ( intImgA2 - (intImgA.^2)/pSzT ) / (pSzT-1) );
denomT = std(T(:));
denom = denomT*denomA;

% compute the NCC
s = warning('off', 'MATLAB:divideByZero');
C = num ./ denom;
s = warning('on', 'MATLAB:divideByZero');

% replace the NaN (if any) with 0's
zeroInd = find(denomA==0);
C(zeroInd) = 0;

switch( lower(shape) )
    case 'full'
    case 'same'
        szTp = fix((szT-1)/2);
        C = C( szTp(1)+1:szTp(1)+szA(1), szTp(2)+1:szTp(2)+szA(2), szTp(3)+1:szTp(3)+szA(3) );
    case 'valid'
        C = C(szT(1):end-szT(1)+1,szT(2):end-szT(2)+1,szT(3):end-szT(3)+1);
    otherwise
        error(sprintf('unknown SHAPE %s, assuming FULL by default', shape));
end

end

function C = bigPhaseCorr3(T, A, shape)

szOut = [(size(T,1)+size(A,1)-1), (size(T,2)+size(A,2)-1), (size(T,3)+size(A,3)-1)]; % szOut = size(intImgA);

% compute the numerator of the NCC
% emulate 3D correlation by rotating templates dimensions in 3D
% requency-domain correlation is MUCH faster than the spatial-domain variety
rotT = flipdim(flipdim(flipdim(T,1),2),3); % this is rot90 in 3d
fftRotT = fftn(rotT,szOut);
fftA = fftn(A,szOut);

xFreq = -0.5*(szOut(1)-1):1:0.5*(szOut(1)-1);
yFreq = -0.5*(szOut(2)-1):1:0.5*(szOut(2)-1);
zFreq = -0.5*(szOut(3)-1):1:0.5*(szOut(3)-1);
[XXFreq,YYFreq,ZZFreq] = ndgrid(xFreq,yFreq,zFreq);
kSq = (size(T,1)/szOut(1))^2 * (XXFreq.^2 + YYFreq.^2 + ZZFreq.^2);

C = real(ifftn( kSq.*fftA.*fftRotT)); %

end

function integralImageA = integralImage(A,szT)
% this is adapted from Matlab's normxcorr2

szA = size(A);

B = zeros( szA+2*szT-1 );
B( szT(1)+1:szT(1)+szA(1), szT(2)+1:szT(2)+szA(2), szT(3)+1:szT(3)+szA(3) ) = A;

s = cumsum(B,1);
c = s(1+szT(1):end,:,:)-s(1:end-szT(1),:,:);
s = cumsum(c,2);
c = s(:,1+szT(2):end,:)-s(:,1:end-szT(2),:);
s = cumsum(c,3);
integralImageA = s(:,:,1+szT(3):end)-s(:,:,1:end-szT(3));
end

%%
function qfactors = compute_qFactor(cc,qnum)

%get peak locations and cc_min maps (i.e. cc - cc(min))
[peak,cc_min] = cellfun(@(x) cc_max_find(double(x)),cc.A,'UniformOutput',0);

%compute two primary quality metrics, as given in "Xue Z, Particle Image
% Velocimetry Correlation Signal-to-noise Metrics, Particle Image
% Pattern Mutual Information and Measurement uncertainty Quantification.
% MS Thesis, Virginia Tech, 2014.

%peak to corr. energy ratio
pce = cellfun(@(x,y) (abs(y)^2)/(1/numel(x)*(sum(abs(x(:)).^2))),cc_min,peak,'UniformOutput',0);
%min value -> 1 (worst case)

%peak to entropy ratio
ppe = cellfun(@(x) q_entropy(double(x)),cc_min,'UniformOutput',0);%peak to cc (information) entropy
%min value -> 0 (worst case)

qfactors = cell2mat(...
    cellfun(@(x,y) [x(:);y(:)], pce,ppe,'UniformOutput',0))';

end

function qfactors = compute_qFactorMat(A,qnum)

%get peak locations and cc_min maps (i.e. cc - cc(min))
[peak,cc_min] = cc_max_find(A);

%compute two primary quality metrics, as given in "Xue Z, Particle Image
% Velocimetry Correlation Signal-to-noise Metrics, Particle Image
% Pattern Mutual Information and Measurement uncertainty Quantification.
% MS Thesis, Virginia Tech, 2014.

%peak to corr. energy ratio
pce = (abs(peak)^2)/(1/numel(cc_min)*(sum(abs(cc_min(:)).^2)));
%min value -> 1 (worst case)

%peak to entropy ratio
ppe = q_entropy(double(cc_min));%peak to cc (information) entropy
%min value -> 0 (worst case)

qfactors = [pce(:);ppe(:)]';

end


function [peak,cc_min] = cc_max_find(cc)
%find the peak and zero-adjusted cc map

cc_min = cc - min(cc(:));%zero-adjust
% cc_filt = imgaussfilt3(cc_min); %filter to remove noise from peak value

[peak,~] = max(cc_min(:)); %get the index of the peak


end

function [ppe] = q_entropy(cc_min)
%compute entropy q-factor for a given cc map

[cc_hist,~] = histcounts(cc_min,30); %get histogram values

entropy = 0;
p = cc_hist/sum(cc_hist); %compute probablities
for i = 1:length(p)%compute entropy
    if p(i) == 0
        entropy = entropy+p(i);
    else
        entropy = entropy+p(i)*log(1/p(i));
    end
end

ppe = 1/entropy; %peak to cc (information) entropy
%min value -> 0 (worst case)


end


