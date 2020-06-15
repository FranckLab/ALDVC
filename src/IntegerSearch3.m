function [xyz,uvw,cc] = IntegerSearch3(Img,DVCpara)
               
gridRange = DVCpara.gridRange;
winsize = DVCpara.winsize;
winstepsize = DVCpara.winstepsize;
fftmethod = DVCpara.fftmethod; 
 
% Input initial integer search zone size
InitialGuessSatisfied = 1;  % gridxBackup = gridxROIRange; gridyBackup = gridyROIRange; gridzBackup = gridzROIRange;

% ========= Start FFT to find initial guess with integer accuracy =======
while InitialGuessSatisfied == 1
    
    switch fftmethod
        case{'bigphasecorr','bigxcorr'}
            prompt = '--- What is your initial guess search zone size? ---';
            tempSizeOfSearchRegion = input(prompt);
            %tempSizeOfSearchRegion = [10*ImgSeqNum,10*ImgSeqNum,3]; %JY!!!
            if length(tempSizeOfSearchRegion) == 1, tempSizeOfSearchRegion=tempSizeOfSearchRegion*ones(1,3); end
        otherwise
            tempSizeOfSearchRegion = 0;
    end
    
    % Start FFT-research
    [xyz,uvw,cc] = funIntegerSearch3(Img,tempSizeOfSearchRegion,gridRange,winsize,winstepsize,fftmethod);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Have a look at integer search
    figure,imagesc3D(uvw.u);
    figure,imagesc3D(uvw.v);
    figure,imagesc3D(uvw.w);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[coordinatesFEM,elementsFEM,~,~,xyz] = MeshSetUp3(xyz,winstepsize);
    %U0 = full(Init3(uvw,xyz)); Plotdisp_show3(U0,coordinatesFEM,elementsFEM); %PlotuvInit3;
    
    switch fftmethod
        case{'bigphasecorr','bigxcorr'}
           prompt = 'Are you satisfied with initial guess with current search region? (0-yes; 1-no)';
           % InitialGuessSatisfied = input(prompt);
           InitialGuessSatisfied = 0; %JY!!! 
        otherwise
           InitialGuessSatisfied = 0;
    end
    
    
end
 
% ======== Output final search region radius ========
cc.SizeOfFFTSearchRegion = tempSizeOfSearchRegion;

% close all
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Have a look at integer search
% figure,imagesc3D(uvw.u);
% figure,imagesc3D(uvw.v);
% figure,imagesc3D(uvw.w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
