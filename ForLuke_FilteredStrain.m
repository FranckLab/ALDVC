%% Section 8
fprintf('------------ Section 8 Start ------------ \n')
 
% % ------ Delete temp files ------
% for tempi = 1:ALSolveStep-1
%     file_name_Subpb1 = ['Subpb1_step',num2str(tempi),'.mat'];
%     file_name_Subpb2 = ['Subpb2_step',num2str(tempi),'.mat'];
%     file_name_dual = ['uvdual_step',num2str(tempi),'.mat'];
%     delete(file_name_Subpb1); delete(file_name_Subpb2); delete(file_name_dual);
% end

% ------ Plotting and Compute Strain-------
M = 1 + (coordinatesFEM(end,1)-coordinatesFEM(1,1))/winstepsize(1);
N = 1 + (coordinatesFEM(end,2)-coordinatesFEM(1,2))/winstepsize(2);
L = 1 + (coordinatesFEM(end,3)-coordinatesFEM(1,3))/winstepsize(3);
XList = coordinatesFEM(1,1):winstepsize(1):coordinatesFEM(end,1);
YList = coordinatesFEM(1,2):winstepsize(2):coordinatesFEM(end,2);
ZList = coordinatesFEM(1,3):winstepsize(3):coordinatesFEM(end,3);
[XX,YY,ZZ] = ndgrid(XList,YList,ZList);
xyz0.x = XX; xyz0.y = YY; xyz0.z = ZZ;

if size(USubpb2,1) == 1
    ULocal = full(USubpb2_New.USubpb2); USubpb2 = full(ULocal); FLocal = full(FSubpb2.FSubpb2); 
else
    ULocal = full(USubpb2); FLocal = full(FSubpb2);
end

Plotdisp_show3(full(ULocal),coordinatesFEM,elementsFEM);
% ------ Smooth displacements ------
prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
DoYouWantToSmoothOnceMore = input(prompt); DispFilterSize=0; DispFilterStd=1;
while DoYouWantToSmoothOnceMore == 0
    ULocal = funSmoothDisp3(ULocal,coordinatesFEM,elementsFEM,winstepsize,DispFilterSize,DispFilterStd);
    close all; Plotdisp_show3(full(ULocal),coordinatesFEM,elementsFEM); % Plotuv(ULocal,x0,y0); 
    DoYouWantToSmoothOnceMore = input(prompt);
end
% % ------ Choose strain computation method ------
% prompt = 'What method to use to compute strain? (0-None; 1-Finite difference; 2-Plane fitting)';
% MethodToComputeStrain = input(prompt);
% while (MethodToComputeStrain ~= 0) && (MethodToComputeStrain ~= 1) && (MethodToComputeStrain ~= 2)
%     disp('Wrong input!')
%     prompt = 'What method to use to compute strain? (0-None; 1-Finite difference; 2-Plane fitting)';
%     MethodToComputeStrain = input(prompt);
% end
% ----- Compute strain field ------
ComputeStrain3;

close all;
% ------ Plot disp ------
% Plotdisp_show(UWorld,coordinatesFEMWorld,elementsFEM);
Plotdisp03(ULocal,coordinatesFEM,elementsFEM);
Plotdisp03(U0,coordinatesFEM,elementsFEM);
% ------ Plot strain ------
Plotstrain_show3(FLocal,coordinatesFEM,elementsFEM);
Plotstrain03(full(FStraintemp),xyz0.x(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad),xyz0.y(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad), ...
    xyz0.z(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad),size(Img{1})); 
%print(['Stretch',num2str(1),'_WS',num2str(21),'_ST',num2str(8),'_F11'],'-dpdf');

% % ------- Add filter and plot strain field -------
Plotstrain_Fij;
%caxis auto; load('colormap_RdYlBu.mat'); colormap(cMap)
Plotstrain03(FStraintemp,xyz0.x(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad),xyz0.y(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad), ...
    xyz0.z(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad),size(Img{1}),'Individual');

fprintf('------------ Section 8 Done ------------ \n \n')

figure,imagesc3D(reshape(USubpb1(1:3:end),M,N,L));
figure,imagesc3D(reshape(USubpb1(2:3:end),M,N,L));
figure,imagesc3D(reshape(USubpb1(3:3:end),M,N,L));
for tempi = 1:ALSolveStep,figure,imagesc3D(reshape(ConvItPerEle(:,tempi),M,N,L)); end
 figure,imagesc3D(reshape(FStraintemp(1:9:end),M-2*Rad,N-2*Rad,L-2*Rad)); caxis auto;
 figure,imagesc3D(reshape(FStraintemp(3:9:end),M-2*Rad,N-2*Rad,L-2*Rad)); caxis auto;
 figure,imagesc3D(reshape(FStraintemp(9:9:end),M-2*Rad,N-2*Rad,L-2*Rad)); caxis auto;
 
% save results_iden_20190504_cut_ws32_st8_clusterNo1_rmlocbadpt.mat ResultDisp ResultDefGrad Resultbeta ...
% ResultConvItPerEle ResultcoordinatesFEM ResultelementsFEM ResultSizeOfFFTSearchReg
% save results_Serena_HighRes_ws32_st8.mat ResultDisp ResultDefGrad Resultbeta ...
% ResultConvItPerEle ResultcoordinatesFEM ResultelementsFEM ResultSizeOfFFTSearchReg
% save results_S14L1L3L5_ws10_st10_clusterNo1.mat ResultDisp ResultDefGrad Resultbeta ...
% ResultConvItPerEle ResultcoordinatesFEM ResultelementsFEM ResultSizeOfFFTSearchReg
