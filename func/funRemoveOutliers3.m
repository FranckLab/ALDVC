function [uvw,cc,RemoveOutliersList] = funRemoveOutliers3(uvw,cc,qDICOrNot,Thr0)
% =========================================================================
% removes outliers using the universal
% outlier test based on
%
% J. Westerweel and F. Scarano. Universal outlier detection for PIV data.
% Exp. Fluids, 39(6):1096{1100, August 2005. doi: 10.1007/s00348-005-0016-6
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
% needs medFilt3 and John D'Errico's inpaint_nans3 
% (http://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint-nans)function. 
% =========================================================================

u = uvw.u; v = uvw.v; w = uvw.w;
DIM = 3; [M,N,L] = size(u);  % size of the displacement field
mSize = [M*N*L,1];

% ============== qDIC bad points removal ===============
% prompt = 'Input threshold for qDIC bad points removal:';
% ccThreshold = input(prompt); % Thr = 50;
% 
if qDICOrNot == 1
    [cc, ccMask_] = removeBadCorrelations(cc,cc.ccThreshold,1,mSize);
    for ii = 1:2
        ccMask{ii} = reshape(double(ccMask_(:,ii)),mSize);
    end
    qDICpceRmList = ismissing(ccMask{1});
    qDICppeRmList = ismissing(ccMask{2});
else 
    qDICpceRmList = []; 
    qDICppeRmList = [];
end
% figure, imagesc3D(qDICpceRmList); caxis auto; colorbar;
% prompt = 'Input threshold for median test:';
% Thr = input(prompt); % Thr = 50;

% figure, imagesc3D(qDICppeRmList); caxis auto; colorbar;
% prompt = 'Input threshold for median test:';
% Thr = input(prompt); % Thr = 50;
% 
% end

% ============== median test bad points removal ===============
% MedianRes = zeros(M,N,L); % initialize medianresidual
% NormFluct = zeros(M,N,L,DIM); % initialize normalized functuration
% b = 1;          % data-point neighborhood radius (commonly set to 1 or 2)
% DVCeps = 0.1;      % estimated measurement noise level (in pixel units)
% 
% for c = 1:3 % loop over the three velocity components
%     if c == 1; VelComp = u; elseif c == 2; VelComp = v; else VelComp = w; end
%     
%     % loop over all the data points (excluding border points)
%     for tempk = 1+b : L-b
% 	for tempi = 1+b : N-b
% 	for tempj = 1+b : M-b
%         Neigh = VelComp(tempj-b:tempj+b, tempi-b:tempi+b, tempk-b:tempk+b); % data neighborhood with center-point
%         NeighCol = Neigh(:); % in column format
%         NeighCol2 = [NeighCol(1: (2*b+1)^2*b+(2*b+1)*b+b); NeighCol((2*b+1)^2*b+(2*b+1)*b+b+2:end)];
%                             % neighborhood excluding center-point
%         MedianNeigh = median(NeighCol2); % median of the neighborhood
%         Fluct = VelComp(tempj,tempi,tempk) - MedianNeigh; % fluctuation with respect to median
%         Res = NeighCol2 - MedianNeigh;       % residual: neighborhood flucturation w.r.t median
%         MedianRes = median(abs(Res));   % median (absolute) value of residual
%         NormFluct(tempj,tempi,tempk,c) = abs(Fluct/(MedianRes+DVCeps)); % normalized fluctuation w.r.t. neighborhood mesian residual
%     end
% 	end
%     end
%     
% end
% NormFluctAbs = sqrt(NormFluct(:,:,:,1).^2 + NormFluct(:,:,:,2).^2 + NormFluct(:,:,:,3).^2);
medianU = cell(1,3);
normFluct = cell(1,3);
normFluctMag = zeros(size(u));
epsilon = 0.1;
[medianU{1}, normFluct{1}] = funRemoveOutliers(u,epsilon);
[medianU{2}, normFluct{2}] = funRemoveOutliers(v,epsilon);
[medianU{3}, normFluct{3}] = funRemoveOutliers(w,epsilon);
normFluctMag =  normFluct{1}.^2 + normFluct{2}.^2 + normFluct{3}.^2;
normFluctMag = sqrt(normFluctMag);

MedFilterOrNot = 0;
while MedFilterOrNot < 1
    
    figure, imagesc3D(normFluctMag); caxis auto; colorbar;
    if isempty(Thr0) || (Thr0 == 0) 
        prompt = 'Input threshold for median test:';
        Thr = input(prompt); % Thr = 50;
    else
        Thr = Thr0;
    end
    RemoveOutliersList = find( normFluctMag > Thr); % detection criterion
                
    % ============== remove bad points ===============
    u2 = u; u2(qDICpceRmList) = NaN; u2(qDICppeRmList) = NaN; u2(RemoveOutliersList) = NaN; 
    v2 = v; v2(qDICpceRmList) = NaN; u2(qDICppeRmList) = NaN; v2(RemoveOutliersList) = NaN;
    w2 = w; w2(qDICpceRmList) = NaN; u2(qDICppeRmList) = NaN; w2(RemoveOutliersList) = NaN;
    u2 = inpaint_nans3(u2,1); v2 = inpaint_nans3(v2,1); w2 = inpaint_nans3(w2,1);
    % --------------------------------------
    close all; 
    figure, imagesc3D(u2); colorbar; title('Displacement u','fontweight','normal');  
    figure, imagesc3D(v2); colorbar; title('Displacement v','fontweight','normal');  
    figure, imagesc3D(w2); colorbar; title('Displacement w','fontweight','normal');  
    %JY!!!
    %fprintf('Do you want to redo Median test: 0(Yes, redo it!); 1(No, it is good!)  \n')
    %prompt = 'Input here: ';
    MedFilterOrNot = 1; %MedFilterOrNot = input(prompt);   
    
end

% ============== Manual bad points removal ===============
% Find some bad inital guess points
ClearBadInitialPointsOrNot = 1;
% prompt = 'Do you clear bad points by setting upper/lower bounds? (0-yes; 1-no)';
% ClearBadInitialPointsOrNot = input(prompt);

while ClearBadInitialPointsOrNot == 0
    
    prompt = 'What is your upper bound for x-displacement?';
    upperbound = input(prompt);
    [row1,~] = find(u2(:)>upperbound);
    prompt = 'What is your lower bound for x-displacement?';
    lowerbound = input(prompt);
    [row2,~] = find(u2(:)<lowerbound);
    prompt = 'What is your upper bound for y-displacement?';
    upperbound = input(prompt);
    [row3,~] = find(v2(:)>upperbound);
    prompt = 'What is your lower bound for y-displacement?';
    lowerbound = input(prompt);
    [row4,~] = find(v2(:)<lowerbound);
    prompt = 'What is your upper bound for z-displacement?';
    upperbound = input(prompt);
    [row5,~] = find(w2(:)>upperbound);
    prompt = 'What is your lower bound for z-displacement?';
    lowerbound = input(prompt);
    [row6,~] = find(w2(:)<lowerbound);
    
%     prompt = 'What is the bound for correlation function Phi?';
%     lowerbound = input(prompt);
%     [row7,col7,zrow7] = find(Phi<lowerbound);
    row7 = []; col7 = []; zrow7 = [];
    row  = [row1; row2; row3; row4; row5; row6; row7]; 
    [rowsub,colsub,zrowsub] = ind2sub([M,N,L], row);
     RemoveOutliersList = [RemoveOutliersList;rowsub];
    for tempi = 1:length(rowsub)
        u2(rowsub(tempi),colsub(tempi),zrowsub(tempi))=NaN; 
        v2(rowsub(tempi),colsub(tempi),zrowsub(tempi))=NaN;
        w2(rowsub(tempi),colsub(tempi),zrowsub(tempi))=NaN;
    end
     
	% --------------------------------------
    % Take a look at solved u, v, and w
    % --------------------------------------
     close all;
    % Plotuvw([],u2,v2,w2,x,y,z);
    figure,imagesc3D(u2); colorbar; title('Displacement u','fontweight','normal');  
    figure,imagesc3D(v2); colorbar; title('Displacement v','fontweight','normal');  
    figure,imagesc3D(w2); colorbar; title('Displacement w','fontweight','normal');  
    % --------------------------------------
    
    prompt = 'Do you clear bad points by setting upper/lower bounds? (0-yes; 1-no)';
    ClearBadInitialPointsOrNot = input(prompt);
    
end

u = inpaint_nans3(u2,1);
v = inpaint_nans3(v2,1);
w = inpaint_nans3(w2,1);
% u = inpaintn(u2); v = inpaintn(v2); w = inpaintn(w2);

uvw.u = u; uvw.v = v; uvw.w = w;
 

disp('****** Finish removing outliers! ******');

end

%% ========================================================================
function [medianU, normFluct] = funRemoveOutliers(u,epsilon)

nSize = 3*[1 1 1];
skipIdx = ceil(prod(nSize)/2);
padOption = 'replicate';

u =  inpaint_nans3(double(u),0);

medianU = medFilt3(u,nSize,padOption,skipIdx);
fluct = u - medianU;
medianRes = medFilt3(abs(fluct),nSize,padOption,skipIdx);
normFluct = abs(fluct./(medianRes + epsilon));

end

%% ========================================================================
function Vr = medFilt3(V0,nSize, padoption, skipIdx)
% fast median filter for 3D data with extra options.

if nargin < 4, skipIdx = 0; end
if nargin < 3, padoption = 'symmetric'; end
if nargin < 2, nSize = [3 3 3]; end

nLength = prod(nSize);
if mod(nLength,2) == 1, padSize = floor(nSize/2);
elseif mod(nLength,2) == 0, padSize = [nSize(1)/2-1,nSize(2)/2];
end

if strcmpi(padoption,'none')
    V = V0;
else
    V = (padarray(V0,padSize(1)*[1,1,1],padoption,'pre'));
    V = (padarray(V,padSize(2)*[1,1,1],padoption,'post'));
end

S = size(V);
nLength = prod(nSize)-sum(skipIdx>1);
Vn = single(zeros(S(1)-(nSize(1)-1),S(2)-(nSize(2)-1),S(3)-(nSize(3)-1),nLength));  % all the neighbor

%%
% build the neighboor

i = cell(1,nSize(1)); j = cell(1,nSize(2)); k = cell(1,nSize(3));
for m = 1:nSize(1), i{m} = m:(S(1)-(nSize(1)-m)); end
for m = 1:nSize(2), j{m} = m:(S(2)-(nSize(2)-m)); end
for m = 1:nSize(3), k{m} = m:(S(3)-(nSize(3)-m)); end

p = 1;
for m = 1:nSize(1)
    for n = 1:nSize(2)
        for o = 1:nSize(3)
            if p ~= skipIdx || skipIdx == 0
                Vn(:,:,:,p) = V(i{m},j{n},k{o});
            end
            p = p + 1;
        end
    end
end

if skipIdx ~= 0, Vn(:,:,:,skipIdx) = []; end
% perform the processing
Vn = sort(Vn,4);

if mod(nLength,2) == 1 % if odd get the middle element
    Vr = Vn(:,:,:,ceil(nLength/2));
else % if even get the mean of the two middle elements
    Vr = mean(cat(4,Vn(:,:,:,nLength/2),Vn(:,:,:,nLength/2+1)),4);
end

end 

%% ========================================================================
function [cc, ccMask] = ...
    removeBadCorrelations(cc,ccThreshold,sizeChange,mSize)

if sizeChange == 1
    %recompute threshold, only use pce & ppe since these give the best
    %results emprically.
    for ii = 1:2
        
        [qf_para{ii},single_distro] = bimodal_gauss_fit(cc.qfactors(:,ii));
        
        if single_distro == 0%(qf_para{ii}(2) + 2*qf_para{ii}(4)) < (qf_para{ii}(3) - 2*qf_para{ii}(5))
            cc.q_thresh{ii} = qf_para{ii}(3) - ccThreshold*qf_para{ii}(5);
        elseif single_distro == 1
            cc.q_thresh{ii} = qf_para{ii}(3) - ccThreshold*qf_para{ii}(5);
        else
            cc.q_thresh{ii} = qf_para{ii}(3) - ccThreshold*qf_para{ii}(5);
        end
    end
    q_trim = [cc.q_thresh{1};cc.q_thresh{2}];
else
    q_trim = [cc.q_thresh{1};cc.q_thresh{2}];
end

%NaN the qfactor values that are below the threshold
temp = bsxfun(@le,cc.qfactors(:,1:2),q_trim');
qfactors_accept = cc.qfactors(:,1:2);
qfactors_accept(temp) = NaN;

for ii = 1:2
    cc.qfactors_accept{ii} = reshape(double(qfactors_accept(:,ii)),mSize);
end

ccMask = ones(size(qfactors_accept)) + ...
    zeros(size(qfactors_accept)).*qfactors_accept;

end



%% ========================================================================
function [paramEsts,single_distro] = bimodal_gauss_fit(x)
%This function takes a dataset and fits a bimodal Gaussian distro to it.

x = sort(x);

%set function for bimodal Gaussian
pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
    p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);
pdf_single = @(x,mu1,sigma1) ...
    normpdf(x,mu1,sigma1);

%starting params, biased mixture toward "good" values,
%centered at quartiles, equal std dev.
pStart = 0.25;
muStart = quantile(x,[.10 .75]);
sigmaStart(1) = sqrt(var(x(1:round(length(x)/5))));
%- 0.25*diff(quantile(x,[0.01 0.25])).^2);
sigmaStart(2) = sqrt(var(x(ceil(length(x)/10):ceil(3*length(x)/4))));
%... - 0.25*diff(quantile(x,[0.25 0.75])).^2);%1:round(length(x)/2)
start = [pStart muStart sigmaStart];

%set lower and upper bounds
lb = [0 -inf -inf 0.00001 0.00001];
ub = [1 inf inf inf inf];

%do the parameter estimation
options = statset('MaxIter',1800, 'MaxFunEvals',3600);
% options.FunValCheck = 'off';
try
    single_distro = 0;
    paramEsts = mle(x, 'pdf',pdf_normmixture, 'start',start, ...
        'lower',lb, 'upper',ub, 'options',options);%,'optimfun','fmincon'
    
    if paramEsts(2)-paramEsts(4) >= paramEsts(3)+paramEsts(5) || ...
            paramEsts(2)+paramEsts(4) <= paramEsts(3)-paramEsts(5)
        
        single_distro = 1;
        %     disp('Parameters estimated for single peak Gaussian')
        paramEsts = mle(x,'options',options);%,'optimfun','fmincon'
        paramEsts = [0.5,paramEsts(1),paramEsts(1),paramEsts(2),...
            paramEsts(2)];
        
    end
    
catch
    single_distro = 1;
    %     disp('Parameters estimated for single peak Gaussian')
    paramEsts = mle(x,'options',options);%,'optimfun','fmincon'
    paramEsts = [0.5,paramEsts(1),paramEsts(1),paramEsts(2),...
        paramEsts(2)];
end

% % %show the result
% % figure
% % % [~, bins] =
% % histogram(x,100);
% % % bins = -2.5:.5:7.5;
% % % h = bar(bins,histc(x,bins)/(length(x)*0.5),'histc');
% % % histogram(x,100)
% % % h.FaceColor = [0.9 0.9 0.9];
% % xgrid = linspace(1.1*min(x),1.1*max(x),200);
% % pdfgrid = pdf_normmixture(xgrid,paramEsts(1),paramEsts(2),paramEsts(3),...
% %     paramEsts(4),paramEsts(5));
% % hold on
% % plot((paramEsts(3) - 2*paramEsts(5)),pdfgrid,'or')
% % plot((paramEsts(2) + 2*paramEsts(4)),pdfgrid,'*r')
% % plot(xgrid,pdfgrid,'-b')
% % hold off
% % xlabel('x')
% % ylabel('Probability Density')

end
