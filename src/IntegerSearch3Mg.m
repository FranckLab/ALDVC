function [xyz,uvw,cc,SizeOfFFTSearchRegion] = IntegerSearch3Mg(Img,DVCpara)

gridRange = DVCpara.gridRange;
winsize = DVCpara.winsize;
winstepsize = DVCpara.winstepsize;
fftmethod = DVCpara.InitFFTMethod; 
ClusterNo = DVCpara.ClusterNo;

% Input initial integer search zone size
% gridxBackup = gridxROIRange; gridyBackup = gridyROIRange; gridzBackup = gridzROIRange;
% ========= Start FFT to find initial guess with integer accuracy =======
tempSizeOfSearchRegion = 0;
% Start FFT-research
if strcmp(fftmethod,'bigxcorrUni') == 1
    disp('--- What is your initial guess search zone size? ---')
    prompt = 'Input here: ';
    tempSizeOfSearchRegion = input(prompt);
    if length(tempSizeOfSearchRegion) == 1
        tempSizeOfSearchRegion = tempSizeOfSearchRegion*ones(1,3);
    end
    [xyz,uvw,cc] = funIntegerSearch3(Img,tempSizeOfSearchRegion,gridRange,winsize,winstepsize,'bigxcorr',ClusterNo);
elseif strcmp(fftmethod,'bigphasecorrUni') == 1
    disp('--- What is your initial guess search zone size? ---')
    prompt = 'Input here: ';
    tempSizeOfSearchRegion = input(prompt);
    if length(tempSizeOfSearchRegion) == 1
        tempSizeOfSearchRegion = tempSizeOfSearchRegion*ones(1,3);
    end
    [xyz,uvw,cc] = funIntegerSearch3(Img,tempSizeOfSearchRegion,gridRange,winsize,winstepsize,'bigphasecorr',ClusterNo);
else
    [xyz,uvw,cc] = funIntegerSearch3Mg(Img,gridRange,winsize,winstepsize,fftmethod);
end
%[coordinatesFEM,elementsFEM,~,~,xyz] = MeshSetUp3(xyz,winstepsize);
%U0 = full(Init3(uvw,xyz)); Plotdisp_show3(U0,coordinatesFEM,elementsFEM); %PlotuvInit3;

% ======== Output final search region radius ========
SizeOfFFTSearchRegion = tempSizeOfSearchRegion;

% close all
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Have a look at integer search
% figure,imagesc3D(uvw.u);
% figure,imagesc3D(uvw.v);
% figure,imagesc3D(uvw.w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

