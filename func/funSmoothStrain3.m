% ==============================================
% function funSmoothDispCrack
% ==============================================
function F = funSmoothStrain3(FLocal,DVCmesh,DVCpara)
 
coordinatesFEM = DVCmesh.coordinatesFEM;
elementsFEM = DVCmesh.elementsFEM;
winstepsize = DVCpara.winstepsize;
StrainFilterSize = DVCpara.StrainFilterSize;
StrainFilterStd = DVCpara.StrainFilterStd; 

FilterStd = StrainFilterStd; FilterSizeInput = StrainFilterSize; LevelNo = 1;
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
    
    Iblur_11 = griddata(coordinatesFEM(:,1),coordinatesFEM(:,2),coordinatesFEM(:,3),FLocal(1:9:end),Xq,Yq,Zq,'natural');
    Iblur_21 = griddata(coordinatesFEM(:,1),coordinatesFEM(:,2),coordinatesFEM(:,3),FLocal(2:9:end),Xq,Yq,Zq,'natural');
    Iblur_31 = griddata(coordinatesFEM(:,1),coordinatesFEM(:,2),coordinatesFEM(:,3),FLocal(3:9:end),Xq,Yq,Zq,'natural');
    Iblur_12 = griddata(coordinatesFEM(:,1),coordinatesFEM(:,2),coordinatesFEM(:,3),FLocal(4:9:end),Xq,Yq,Zq,'natural');
    Iblur_22 = griddata(coordinatesFEM(:,1),coordinatesFEM(:,2),coordinatesFEM(:,3),FLocal(5:9:end),Xq,Yq,Zq,'natural');
    Iblur_32 = griddata(coordinatesFEM(:,1),coordinatesFEM(:,2),coordinatesFEM(:,3),FLocal(6:9:end),Xq,Yq,Zq,'natural');
    Iblur_13 = griddata(coordinatesFEM(:,1),coordinatesFEM(:,2),coordinatesFEM(:,3),FLocal(7:9:end),Xq,Yq,Zq,'natural');
    Iblur_23 = griddata(coordinatesFEM(:,1),coordinatesFEM(:,2),coordinatesFEM(:,3),FLocal(8:9:end),Xq,Yq,Zq,'natural');
    Iblur_33 = griddata(coordinatesFEM(:,1),coordinatesFEM(:,2),coordinatesFEM(:,3),FLocal(9:9:end),Xq,Yq,Zq,'natural');
    
    % Gaussian smooth deformation gradients
%     F11temp = imgaussfilt3(Iblur_11,FilterStd,'FilterSize',FilterSizeInput);
%     F21temp = imgaussfilt3(Iblur_21,FilterStd,'FilterSize',FilterSizeInput);
%     F31temp = imgaussfilt3(Iblur_31,FilterStd,'FilterSize',FilterSizeInput);
%     F12temp = imgaussfilt3(Iblur_12,FilterStd,'FilterSize',FilterSizeInput);
%     F22temp = imgaussfilt3(Iblur_22,FilterStd,'FilterSize',FilterSizeInput);
%     F32temp = imgaussfilt3(Iblur_32,FilterStd,'FilterSize',FilterSizeInput);
%     F13temp = imgaussfilt3(Iblur_13,FilterStd,'FilterSize',FilterSizeInput);
%     F23temp = imgaussfilt3(Iblur_23,FilterStd,'FilterSize',FilterSizeInput);
%     F33temp = imgaussfilt3(Iblur_33,FilterStd,'FilterSize',FilterSizeInput);

    F11temp = medfilt3(Iblur_11,[FilterSizeInput,FilterSizeInput,FilterSizeInput]);
    F21temp = medfilt3(Iblur_21,[FilterSizeInput,FilterSizeInput,FilterSizeInput]);
    F31temp = medfilt3(Iblur_31,[FilterSizeInput,FilterSizeInput,FilterSizeInput]);
    F12temp = medfilt3(Iblur_12,[FilterSizeInput,FilterSizeInput,FilterSizeInput]);
    F22temp = medfilt3(Iblur_22,[FilterSizeInput,FilterSizeInput,FilterSizeInput]);
    F32temp = medfilt3(Iblur_32,[FilterSizeInput,FilterSizeInput,FilterSizeInput]);
    F13temp = medfilt3(Iblur_13,[FilterSizeInput,FilterSizeInput,FilterSizeInput]);
    F23temp = medfilt3(Iblur_23,[FilterSizeInput,FilterSizeInput,FilterSizeInput]);
    F33temp = medfilt3(Iblur_33,[FilterSizeInput,FilterSizeInput,FilterSizeInput]);
    
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
    F=zeros(9*size(coordinatesFEM,1),1);
    for tempi = 1:size(coordinatesFEM,1)
        [row1,~] = find(Coordxnodes==coordinatesFEM((tempi),1));
        [row2,~] = find(Coordynodes==coordinatesFEM((tempi),2));
        [row3,~] = find(Coordznodes==coordinatesFEM(tempi,3));
        F(9*tempi-8) = F11temp(row1,row2,row3); F(9*tempi-7) = F21temp(row1,row2,row3); F(9*tempi-6) = F31temp(row1,row2,row3);
        F(9*tempi-5) = F12temp(row1,row2,row3); F(9*tempi-4) = F22temp(row1,row2,row3); F(9*tempi-3) = F32temp(row1,row2,row3);
        F(9*tempi-2) = F13temp(row1,row2,row3); F(9*tempi-1) = F23temp(row1,row2,row3); F(9*tempi-0) = F33temp(row1,row2,row3);
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