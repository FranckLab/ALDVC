function [xyz,uvw,cc,sizeOfFFTSearchRegion] = IntegerSearch3Multigrid(Img,DVCpara)

gridRange = DVCpara.gridRange;
winsize = DVCpara.winsize;
winstepsize = DVCpara.winstepsize;
fftMethod = DVCpara.initFFTMethod; 
clusterNo = DVCpara.clusterNo;

% Input initial integer search zone size
% gridxBackup = gridxROIRange; gridyBackup = gridyROIRange; gridzBackup = gridzROIRange;
% ========= Start FFT to find initial guess with integer accuracy =======
tempSizeOfSearchRegion = 0;
% Start FFT-research
if strcmp(fftMethod,'bigxcorrUni') == 1
    disp('--- What is your initial guess search zone size? ---')
    prompt = 'Input here: ';
    tempSizeOfSearchRegion = input(prompt);
    if length(tempSizeOfSearchRegion) == 1
        tempSizeOfSearchRegion = tempSizeOfSearchRegion*ones(1,3);
    end
    [xyz,uvw,cc] = funIntegerSearch3(Img,tempSizeOfSearchRegion,gridRange,winsize,winstepsize,'bigxcorr',clusterNo);

elseif strcmp(fftMethod,'bigphasecorrUni') == 1
    disp('--- What is your initial guess search zone size? ---')
    prompt = 'Input here: ';
    tempSizeOfSearchRegion = input(prompt);
    if length(tempSizeOfSearchRegion) == 1
        tempSizeOfSearchRegion = tempSizeOfSearchRegion*ones(1,3);
    end
    [xyz,uvw,cc] = funIntegerSearch3(Img,tempSizeOfSearchRegion,gridRange,winsize,winstepsize,'bigphasecorr',clusterNo);

elseif strcmp(fftMethod,'xcorr') == 1
    [xyz,uvw,cc] = funIntegerSearch3(Img,0,gridRange,winsize,winstepsize,'xcorr',clusterNo);

elseif strcmp(fftMethod,'phasecorr') == 1
    [xyz,uvw,cc] = funIntegerSearch3(Img,0,gridRange,winsize,winstepsize,'phasecorr',clusterNo);

elseif strcmp(fftMethod,'bigphasecorr') == 1
    [xyz,uvw,cc] = funIntegerSearch3Multigrid(Img,gridRange,winsize,winstepsize,'phasecorr');
    
else % if strcmp(fftmethod,'bigxcorr') == 1
    [xyz,uvw,cc] = funIntegerSearch3Multigrid(Img,gridRange,winsize,winstepsize,'xcorr');

end
%[coordinatesFEM,elementsFEM,~,~,xyz] = MeshSetUp3(xyz,winstepsize);
%U0 = full(Init3(uvw,xyz)); Plotdisp_show3(U0,coordinatesFEM,elementsFEM); %PlotuvInit3;

% ======== Output final search region radius ========
sizeOfFFTSearchRegion = tempSizeOfSearchRegion;

% close all
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Have a look at integer search
% figure,imagesc3D(uvw.u);
% figure,imagesc3D(uvw.v);
% figure,imagesc3D(uvw.w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

