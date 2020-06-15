% ==============================================
% function funSmoothDispCrack
% ==============================================
function U = funSmoothDisp3(U,DVCmesh,DVCpara)
 
coordinatesFEM = DVCmesh.coordinatesFEM;
elementsFEM = DVCmesh.elementsFEM;
winstepsize = DVCpara.winstepsize;
DispFilterSize = DVCpara.DispFilterSize;
DispFilterStd = DVCpara.DispFilterStd; 

FilterStd = DispFilterStd; FilterSizeInput = DispFilterSize; LevelNo = 1;
% switch nargin
%     case 5
%         FilterSizeInput = varargin{1};
%     case 6
%         FilterSizeInput = varargin{1}; FilterStd = varargin{2};
%     case 7
%         FilterSizeInput = varargin{1}; FilterStd = varargin{2}; LevelNo = varargin{3};
%     otherwise
%         disp('Wrong input in funSmoothStrain!');
% end

%%
% close all; Plotdisp_show(ULocal,elementsFEM,coordinatesFEM);

% prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
DoYouWantToSmoothOnceMore = 0; % DoYouWantToSmoothOnceMore = input(prompt);
if DoYouWantToSmoothOnceMore == 0  
    if isempty(FilterStd) == 1
        prompt = 'Choose filter standard deviation(0-default): ';
        FilterStd = input(prompt);
        if FilterStd == 0
            FilterStd = 0.5; 
        end
    else
        if FilterStd == 0
            FilterStd = 0.5;
        end
    end
    if isempty(FilterSizeInput) == 1
        prompt = 'Choose Gaussian filter size(0-default): ';
        FilterSizeInput = input(prompt);
        if FilterSizeInput == 0
            FilterSizeInput = 2*ceil(2*FilterStd)+1; 
        end
    else
        if FilterSizeInput == 0
            FilterSizeInput = 2*ceil(2*FilterStd)+1;
        end
    end
end

SmoothTimes = 1;
while (DoYouWantToSmoothOnceMore==0)
    Coordxnodes = [min(coordinatesFEM(:,1)):winstepsize(1)/(2^(LevelNo-1)):max(coordinatesFEM(:,1))]'; 
    Coordynodes = [min(coordinatesFEM(:,2)):winstepsize(2)/(2^(LevelNo-1)):max(coordinatesFEM(:,2))]';
    Coordznodes = [min(coordinatesFEM(:,3)):winstepsize(3)/(2^(LevelNo-1)):max(coordinatesFEM(:,3))]';
    [Xq,Yq,Zq] = ndgrid(Coordxnodes,Coordynodes,Coordznodes);
    
    Iblur_1 = griddata(coordinatesFEM(:,1),coordinatesFEM(:,2),coordinatesFEM(:,3),U(1:3:end),Xq,Yq,Zq,'natural');
    Iblur_2 = griddata(coordinatesFEM(:,1),coordinatesFEM(:,2),coordinatesFEM(:,3),U(2:3:end),Xq,Yq,Zq,'natural');
    Iblur_3 = griddata(coordinatesFEM(:,1),coordinatesFEM(:,2),coordinatesFEM(:,3),U(3:3:end),Xq,Yq,Zq,'natural');
    
    % Gaussian smooth deformation gradients
%     U1temp = imgaussfilt3(Iblur_1,FilterStd,'FilterSize',FilterSizeInput);
%     U2temp = imgaussfilt3(Iblur_2,FilterStd,'FilterSize',FilterSizeInput);
%     U3temp = imgaussfilt3(Iblur_3,FilterStd,'FilterSize',FilterSizeInput);
    U1temp = medfilt3(Iblur_1,[FilterSizeInput,FilterSizeInput,FilterSizeInput]);
    U2temp = medfilt3(Iblur_2,[FilterSizeInput,FilterSizeInput,FilterSizeInput]);
    U3temp = medfilt3(Iblur_3,[FilterSizeInput,FilterSizeInput,FilterSizeInput]);
    
%     Iblur_Top11 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), FLocal(1:4:end),Coordxnodes,Coordynodes,'regularizer','springs'); 
%     Iblur_Top11=Iblur_Top11';
%     Iblur_Top22 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), FLocal(4:4:end),Coordxnodes,Coordynodes,'regularizer','springs');  
%     Iblur_Top22=Iblur_Top22';
%     Iblur_Top21 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), FLocal(2:4:end),Coordxnodes,Coordynodes,'regularizer','springs'); 
%     Iblur_Top21=Iblur_Top21';
%     Iblur_Top12 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), FLocal(3:4:end),Coordxnodes,Coordynodes,'regularizer','springs');
%     Iblur_Top12=Iblur_Top12';
%     % Iblur_Top = nan(size(ULocal,1),1); Iblur_Top(2*CoordCrackTop-1) = ULocal(2*CoordCrackTop-1); Iblur_Top(2*CoordCrackTop) = ULocal(2*CoordCrackTop); 
%     % Iblur_Top10 = reshape(Iblur_Top(1:2:end),M,N); Iblur_Top20 = reshape(Iblur_Top(2:2:end),M,N);
%     % -------------------------------------------------------
%     % Iblur_Top1 = reshape(imgaussfilt(Iblur_Top10,FilterStd,'FilterSize',FilterSizeInput,'FilterDomain','spatial'), M*N, 1);
%     % Iblur_Top2 = reshape(imgaussfilt(Iblur_Top20,FilterStd,'FilterSize',FilterSizeInput,'FilterDomain','spatial'), M*N, 1);
%     % -------------------------------------------------------
%     imageFilter=fspecial('gaussian',FilterSizeInput,FilterStd);
%     % Iblur_Top1 = reshape(nanconv(Iblur_Top10,imageFilter,'edge','nanout'), M*N, 1);
%     % Iblur_Top2 = reshape(nanconv(Iblur_Top20,imageFilter,'edge','nanout'), M*N, 1);
%     % ULocal(2*CoordCrackTop-1) = Iblur_Top1(CoordCrackTop); ULocal(2*CoordCrackTop) = Iblur_Top2(CoordCrackTop);
%     Iblur_Top1 = nanconv(Iblur_Top11,imageFilter,'edge','nanout');
%     Iblur_Top4 = nanconv(Iblur_Top22,imageFilter,'edge','nanout');
%     Iblur_Top2 = nanconv(Iblur_Top21,imageFilter,'edge','nanout');
%     Iblur_Top3 = nanconv(Iblur_Top12,imageFilter,'edge','nanout');
%     
    for tempi = 1:size(coordinatesFEM,1)
        [row1,~] = find(Coordxnodes==coordinatesFEM((tempi),1));
        [row2,~] = find(Coordynodes==coordinatesFEM((tempi),2));
        [row3,~] = find(Coordznodes==coordinatesFEM(tempi,3));
 
        U(3*tempi-2) = U1temp(row1,row2,row3); U(3*tempi-1) = U2temp(row1,row2,row3); U(3*tempi-0) = U3temp(row1,row2,row3);
    end
    
    % close all; Plotuv(ULocal,x,y); % Plotting u and v
    % close all; Plotdisp_show(ULocal,elementsFEM,coordinatesFEM);
     
    % prompt = 'Do you want to smooth displacement once more? (0-yes; 1-no)';
    % DoYouWantToSmoothOnceMore = input(prompt); 
    SmoothTimes = SmoothTimes+1;
    if SmoothTimes > 2
        DoYouWantToSmoothOnceMore = 1;
    end
    
end


