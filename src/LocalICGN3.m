% ==============================================
% function funICGN
% ==============================================
function [U,F,HtempPar,LocalTime,ConvItPerEle] = LocalICGN3(U0,coordinatesFEM,Df,Img1,Img2,DVCpara,ICGNmethod,tol)

winsize = DVCpara.winsize;
winstepsize = DVCpara.winstepsize;
ClusterNo = DVCpara.ClusterNo;
interpmethod = DVCpara.interpmethod;
displayIterOrNot = DVCpara.displayIterOrNot;
MaxIterNum = DVCpara.Subpb1ICGNMaxIterNum;
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
temp = zeros(size(coordinatesFEM,1),1); 
UtempPar = temp; VtempPar = temp; WtempPar = temp;
F11tempPar = temp; F21tempPar = temp; F31tempPar = temp;
F12tempPar = temp; F22tempPar = temp; F32tempPar = temp;
F13tempPar = temp; F23tempPar = temp; F33tempPar = temp;
HtempPar = zeros(size(coordinatesFEM,1),78);
ConvItPerEle = zeros(size(coordinatesFEM,1),1);
% -----------------------------------------------
% Spline interpolation
% imggNormalizedbc  = Spline2D('bicubic',[1:1:Df.imgSize(1)],[1:1:Df.imgSize(2)],imggNormalizedbc);
%
% disp('--- Set up Parallel pool ---');  
% -------- How to change parallel pools ---------
% myCluster = parcluster('local');
% myCluster.NumWorkers = 4;  % 'Modified' property now TRUE
% saveProfile(myCluster);    % 'local' profile now updated,
%                            % 'Modified' property now FALSE
% -------------- Or we can do this --------------
% Go to the Parallel menu, then select Manage Cluster Profiles.
% Select the "local" profile, and change NumWorkers to 4.
% -----------------------------------------------
switch ClusterNo
    case 0 || 1
        h = waitbar(0,'Please wait for Subproblem 1 IC-GN iterations!'); tic;
        for tempj = 1 : size(coordinatesFEM,1)  % tempj is the element index
            
            try
                
                DfEle = struct(); DfEle.imgSize = Df.imgSize;
                xyz1 = coordinatesFEM(tempj,:)-floor(0.5*winsize); xyz7 = coordinatesFEM(tempj,:)+floor(0.5*winsize);
                
                ImgEle = struct();
                ImgEle.Imgf = Img1(xyz1(1)-3:1:xyz7(1)+3, xyz1(2)-3:1:xyz7(2)+3, xyz1(3)-3:1:xyz7(3)+3); % Don't change these +/- 3
                %ImgEle.Imgg = Img2;
                %ImgEle.ImggAxis = [1,size(Img2,1),1,size(Img2,2),1,size(Img2,3)]-1;
                
                [Utemp, Ftemp, ConvItPerEle(tempj), HtempPar(tempj,:)] = funICGN3(U0(3*tempj-2:3*tempj), ...
                    coordinatesFEM(tempj,:),DfEle,ImgEle,Img2,winsize,tol,ICGNmethod,interpmethod,MaxIterNum);
                
            catch
                
                Utemp = nan(3,1);  Ftemp = nan(9,1); ConvItPerEle(tempj) = 0;
                HtempPar(tempj,:) = zeros(1,78);
                
            end
            
            if displayIterOrNot == 1, disp(['ele ',num2str(tempj),' converge or not is ',num2str(ConvItPerEle(tempj)),' (#>0-converged; 0-unconverged)']); end
            
            % ------ Store solved deformation gradients ------
            UtempPar(tempj) = Utemp(1); VtempPar(tempj) = Utemp(2); WtempPar(tempj) = Utemp(3);
            F11tempPar(tempj) = Ftemp(1); F21tempPar(tempj) = Ftemp(2); F31tempPar(tempj) = Ftemp(3);
            F12tempPar(tempj) = Ftemp(4); F22tempPar(tempj) = Ftemp(5); F32tempPar(tempj) = Ftemp(6);
            F13tempPar(tempj) = Ftemp(7); F23tempPar(tempj) = Ftemp(8); F33tempPar(tempj) = Ftemp(9);
            
            waitbar(tempj/(size(coordinatesFEM,1)));
            
        end
        close(h); LocalTime = toc;
        
    otherwise
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --------------- Start Local ICGN --------------
        
        tic;
        
        % === Old version: to pre-store local subset info before Subpb1 ===
        % imgSubsetPar = cell(1,size(coordinatesFEM,1));
        % 
        % hbar = waitbar(0,'Pre-store subset info before parfor Subproblem 1 IC-GN iterations!');
        % for tempj = 1:size(coordinatesFEM,1)
        % 
        %     % Save RAM and pre-store all the image subset info.
        %     DfEle = struct(); DfEle.imgSize = Df.imgSize; DfElePar{:,tempj} = DfEle;
        % 
        %     xyz1 = coordinatesFEM(tempj,:) - 0.5*winsize; xyz7 = xyz1 + winsize;
        %     ImgEle = struct();
        %     ImgEle.Imgf = Img1(xyz1(1)-3:1:xyz7(1)+3, xyz1(2)-3:1:xyz7(2)+3, xyz1(3)-3:1:xyz7(3)+3);
        %     imgSubsetPar{:,tempj} = ImgEle;
        %     waitbar(tempj/(size(coordinatesFEM,1)));
        % 
        % end
        % close(hbar);
        % === End of old version code ===
        
        Img1Const = parallel.pool.Constant(Img1);
        Img2Const = parallel.pool.Constant(Img2);
        DfEle = struct(); DfEle.imgSize = Df.imgSize;
         
        hbar = parfor_progressbar(size(coordinatesFEM,1),'Subproblem 1: Parallel-computing Subproblem 1 IC-GN iterations!');
        parfor tempj = 1 : size(coordinatesFEM,1)  % tempj is the element index
            
            try 
                xyz1 = coordinatesFEM(tempj,:)-floor(0.5*winsize); xyz7 = coordinatesFEM(tempj,:)+floor(0.5*winsize);
                
                ImgEle = struct();
                ImgEle.Imgf = Img1Const.Value(xyz1(1)-3:1:xyz7(1)+3, xyz1(2)-3:1:xyz7(2)+3, xyz1(3)-3:1:xyz7(3)+3); % Don't change these +/- 3
            
                % === Old version: use pre-stored local subset info to do Subpb1 ===
                % [Utemp, Ftemp, ConvItPerEle(tempj), HtempPar(tempj,:)] = funICGN3(U0(3*tempj-2:3*tempj), ...
                %     coordinatesFEM(tempj,:),DfElePar{:,tempj},imgSubsetPar{:,tempj},Img2Const.Value,winsize,tol,ICGNmethod,interpmethod,MaxIterNum);
                % === End of old version code ===
                [Utemp, Ftemp, ConvItPerEle(tempj), HtempPar(tempj,:)] = funICGN3(U0(3*tempj-2:3*tempj), ...
                    coordinatesFEM(tempj,:),DfEle,ImgEle,Img2Const.Value,winsize,tol,ICGNmethod,interpmethod,MaxIterNum);
                
            catch
                
                Utemp = nan(3,1);  Ftemp = nan(9,1); ConvItPerEle(tempj) = 0;
                HtempPar(tempj,:) = zeros(1,78);
                
            end
            
            if displayIterOrNot == 1, disp(['ele ',num2str(tempj),' converge or not is ',num2str(ConvItPerEle(tempj)),' (#>0-converged; 0-unconverged)']); end
            
            % ------ Store solved deformation gradients ------
            UtempPar(tempj) = Utemp(1); VtempPar(tempj) = Utemp(2); WtempPar(tempj) = Utemp(3);
            F11tempPar(tempj) = Ftemp(1); F21tempPar(tempj) = Ftemp(2); F31tempPar(tempj) = Ftemp(3);
            F12tempPar(tempj) = Ftemp(4); F22tempPar(tempj) = Ftemp(5); F32tempPar(tempj) = Ftemp(6);
            F13tempPar(tempj) = Ftemp(7); F23tempPar(tempj) = Ftemp(8); F33tempPar(tempj) = Ftemp(9);
            
            hbar.iterate(1);
            
        end
         
        LocalTime = toc; close(hbar);  
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = U0; U(1:3:end) = UtempPar ; U(2:3:end) = VtempPar ; U(3:3:end) = WtempPar ;
F = zeros(3^2*size(coordinatesFEM,1),1);
F(1:9:end) = F11tempPar; F(2:9:end) = F21tempPar; F(3:9:end) = F31tempPar; F(4:9:end) = F12tempPar;
F(5:9:end) = F22tempPar; F(6:9:end) = F32tempPar; F(7:9:end) = F13tempPar; F(8:9:end) = F23tempPar; F(9:9:end) = F33tempPar;

% ------ Clear bad points for Local DIC ------
% find bad points after Local Subset ICGN
[row1,~] = find(ConvItPerEle(:,1)>MaxIterNum);
[row2,~] = find(ConvItPerEle(:,1)<1);
row = unique(union(row1,row2)); LocalICGNBadPtNum = length(row);
% find bad points with unusual ICGN steps
% rownans = unique(union (row1,row2) ); rowconv = setdiff(1:1:size(coordinatesFEM,1), rownans);
% meanConvItPerEle = mean(ConvItPerEle(rowconv)); stdConvItPerEle = std(ConvItPerEle(rowconv));
% [row3,~] = find(ConvItPerEle(rowconv)>meanConvItPerEle + 0.1*stdConvItPerEle);
% row = unique(union(rownans,row3));
% row = unique(union(row1,row2));

% disp(['Local step bad subsets total # is: ', num2str(length(row))]);
disp(['Local ICGN bad subsets %: ', num2str(LocalICGNBadPtNum),'/',num2str(size(coordinatesFEM,1)),'=',num2str(100*LocalICGNBadPtNum/size(coordinatesFEM,1)),'%']);

U(3*row-2) = NaN; U(3*row-1) = NaN; U(3*row) = NaN;
F(9*row-8) = NaN; F(9*row-7) = NaN; F(9*row-6) = NaN; F(9*row-5) = NaN;
F(9*row-4) = NaN; F(9*row-3) = NaN; F(9*row-2) = NaN; F(9*row-1) = NaN; F(9*row) = NaN;
% U(3*row-2) = 0; U(3*row-1) = 0; U(3*row) = 0;
% Plotdisp_show3(full(U),coordinatesFEM,elementsFEM );
% Plotuv(full(USubpb1),x0,y0);
% ------ inpaint nans using gridfit ------
Coordxnodes = unique(coordinatesFEM(:,1)); Coordynodes = unique(coordinatesFEM(:,2)); Coordznodes = unique(coordinatesFEM(:,3));
nanindex = find(isnan(U(1:3:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
M = length(Coordxnodes);  N = length(Coordynodes);  L = length(Coordznodes);

% --- Apply only for non-grid meshes ---
%[Xq,Yq,Zq] = ndgrid(Coordxnodes,Coordynodes,Coordznodes);
% u1temp = griddata(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),coordinatesFEM(notnanindex,3),U(3*notnanindex-2),Xq,Yq,Zq,'natural');
% v1temp = griddata(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),coordinatesFEM(notnanindex,3),U(3*notnanindex-1),Xq,Yq,Zq,'natural');
% w1temp = griddata(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),coordinatesFEM(notnanindex,3),U(3*notnanindex-0),Xq,Yq,Zq,'natural');
% F11temp = griddata(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),coordinatesFEM(notnanindex,3),F(9*notnanindex-8),Xq,Yq,Zq,'natural');
% F21temp = griddata(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),coordinatesFEM(notnanindex,3),F(9*notnanindex-7),Xq,Yq,Zq,'natural');
% F31temp = griddata(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),coordinatesFEM(notnanindex,3),F(9*notnanindex-6),Xq,Yq,Zq,'natural');
% F12temp = griddata(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),coordinatesFEM(notnanindex,3),F(9*notnanindex-5),Xq,Yq,Zq,'natural');
% F22temp = griddata(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),coordinatesFEM(notnanindex,3),F(9*notnanindex-4),Xq,Yq,Zq,'natural');
% F32temp = griddata(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),coordinatesFEM(notnanindex,3),F(9*notnanindex-3),Xq,Yq,Zq,'natural');
% F13temp = griddata(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),coordinatesFEM(notnanindex,3),F(9*notnanindex-2),Xq,Yq,Zq,'natural');
% F23temp = griddata(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),coordinatesFEM(notnanindex,3),F(9*notnanindex-1),Xq,Yq,Zq,'natural');
% F33temp = griddata(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),coordinatesFEM(notnanindex,3),F(9*notnanindex-0),Xq,Yq,Zq,'natural');
% --- 2D version codes ---
% [u1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); u1temp = u1temp';
% [v1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex),Coordxnodes,Coordynodes,'regularizer','springs'); v1temp = v1temp';
% [F11temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3),Coordxnodes,Coordynodes,'regularizer','springs'); F11temp = F11temp';
% [F21temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2),Coordxnodes,Coordynodes,'regularizer','springs'); F21temp = F21temp';
% [F12temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); F12temp = F12temp';
% [F22temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0),Coordxnodes,Coordynodes,'regularizer','springs'); F22temp = F22temp';
% --- using inpaint_nans3 ---
u1temp = inpaint_nans3(reshape(U(1:3:end),M,N,L),0);
v1temp = inpaint_nans3(reshape(U(2:3:end),M,N,L),0);
w1temp = inpaint_nans3(reshape(U(3:3:end),M,N,L),0);
F11temp = inpaint_nans3(reshape(F(1:9:end),M,N,L),0);
F21temp = inpaint_nans3(reshape(F(2:9:end),M,N,L),0);
F31temp = inpaint_nans3(reshape(F(3:9:end),M,N,L),0);
F12temp = inpaint_nans3(reshape(F(4:9:end),M,N,L),0);
F22temp = inpaint_nans3(reshape(F(5:9:end),M,N,L),0);
F32temp = inpaint_nans3(reshape(F(6:9:end),M,N,L),0);
F13temp = inpaint_nans3(reshape(F(7:9:end),M,N,L),0);
F23temp = inpaint_nans3(reshape(F(8:9:end),M,N,L),0);
F33temp = inpaint_nans3(reshape(F(9:9:end),M,N,L),0);
% --- using inpaintn ---
% u1temp = inpaintn(reshape(U(1:3:end),M,N,L));
% v1temp = inpaintn(reshape(U(2:3:end),M,N,L));
% w1temp = inpaintn(reshape(U(3:3:end),M,N,L));
% F11temp = inpaintn(reshape(F(1:9:end),M,N,L));
% F21temp = inpaintn(reshape(F(2:9:end),M,N,L));
% F31temp = inpaintn(reshape(F(3:9:end),M,N,L));
% F12temp = inpaintn(reshape(F(4:9:end),M,N,L));
% F22temp = inpaintn(reshape(F(5:9:end),M,N,L));
% F32temp = inpaintn(reshape(F(6:9:end),M,N,L));
% F13temp = inpaintn(reshape(F(7:9:end),M,N,L));
% F23temp = inpaintn(reshape(F(8:9:end),M,N,L));
% F33temp = inpaintn(reshape(F(9:9:end),M,N,L));

% Add remove outliers median-test based on:
% u1=cell(2,1); u1{1}=u1temp; u1{2}=v1temp;
% [u2] = removeOutliersMedian(u1,4); u2temp=u2{1}; v2temp=u2{2};
for tempi = 1:size(coordinatesFEM,1)
    if ismember(tempi,nanindex) == 1
        [row1,~] = find(Coordxnodes==coordinatesFEM(tempi,1));
        [row2,~] = find(Coordynodes==coordinatesFEM(tempi,2));
        [row3,~] = find(Coordznodes==coordinatesFEM(tempi,3));
        U(3*tempi-2) = u1temp(row1,row2,row3);
        U(3*tempi-1) = v1temp(row1,row2,row3);
        U(3*tempi)   = w1temp(row1,row2,row3);
        F(9*tempi-8) = F11temp(row1,row2,row3); F(9*tempi-7) = F21temp(row1,row2,row3); F(9*tempi-6) = F31temp(row1,row2,row3);
        F(9*tempi-5) = F12temp(row1,row2,row3); F(9*tempi-4) = F22temp(row1,row2,row3); F(9*tempi-3) = F32temp(row1,row2,row3);
        F(9*tempi-2) = F13temp(row1,row2,row3); F(9*tempi-1) = F23temp(row1,row2,row3); F(9*tempi-0) = F33temp(row1,row2,row3);
    end
end






