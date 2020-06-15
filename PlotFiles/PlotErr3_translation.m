
close all;
ImgNo = 0:1:10;
%% RMS
DispErrLocalICGN_ws30_st30 = [0.0099845     0.018268     0.021438     0.019801     0.015395     0.012623     0.015115 0.019555     0.021382      0.01806     0.012387];
DispErrAL_ws30_st30 = [0.0090062     0.017508     0.020822     0.019085      0.01437     0.011304     0.014079  0.018836     0.020685     0.017352      0.01116];
DispErrLocalICGN_ws20_st20 = [   0.018643     0.026138     0.028532     0.027402     0.024557     0.022796     0.024289 0.027059     0.028447     0.026038     0.022752];
DispErrAL_ws20_st20 = [0.016564     0.024094     0.026639     0.025319     0.022096     0.020103     0.021845 0.025     0.026498     0.023918     0.020206];
DispErrLocalICGN_ws10_st10 = [0.055704     0.066438      0.06736     0.067137     0.066827     0.065991 0.066547      0.06734     0.067231     0.066535     0.065607];
DispErrAL_ws10_st10 = [0.047324     0.059153     0.060224     0.059665     0.058828     0.058332 0.05876     0.059834      0.06003     0.058999     0.058013];

figure(1); clf; %plot(ImgNo,DispErrFFT_ws10_st10,'d--');hold on
plot(0.1*ImgNo,DispErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,DispErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot(0.1*ImgNo,DispErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,DispErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot(0.1*ImgNo,DispErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,DispErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('$x-$direction translation (voxel)','interpreter','latex'); 
ylabel('Displacement RMS error (voxel)','interpreter','latex');

h = legend({'Local DVC WS=30, ST=30','AL-DVC \quad WS=30, ST=30',...
    'Local DVC WS=20, ST=20','AL-DVC \quad WS=20, ST=20',...
    'Local DVC WS=10, ST=10','AL-DVC \quad WS=10, ST=10'},'FontSize',12,'location','northeastoutside','interpreter','latex');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';
axis([0,1,0.007,0.072]); 
set(gca,'fontsize',20); set(gcf,'color','w');
% set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
% export_fig(gcf,'-png','-r300','painters','Fig_trans_disp_err1rms');

%% mean
DispErrLocalICGN_ws30_st30 = [0.00011545    -0.013244    -0.017614    -0.015473   -0.0091457   -0.0002208    0.0087826 0.01519     0.017562     0.013346   1.9713e-05];
DispErrAL_ws30_st30 = [0.00010739    -0.013318    -0.017674    -0.015513    -0.009164  -0.00019285    0.0087742 0.015225     0.017615     0.013405   1.5414e-05];
 DispErrLocalICGN_ws20_st20 = [0.00012886    -0.013059    -0.017363    -0.015422    -0.009101   -4.962e-05    0.0087836 0.015092     0.017434     0.013203   0.00012868];
 DispErrAL_ws20_st20 = [0.00013909    -0.013256    -0.017523    -0.015547   -0.0090977  -5.1205e-05    0.0088073 0.015189     0.017553     0.013291   5.6865e-05];
DispErrLocalICGN_ws10_st10 = [7.3217e-05    -0.011461    -0.015684    -0.014102   -0.0083997   -9.816e-05 0.0082778     0.013824     0.015867     0.011299   0.00027293];
DispErrAL_ws10_st10 = [0.00011738    -0.013005    -0.017167    -0.015277   -0.0090355  -0.00011944 0.008712      0.01494      0.01727     0.012993   0.00010305];

figure(2); clf; %plot(ImgNo,DispErrFFT_ws10_st10,'d--');hold on
plot(0.1*ImgNo,DispErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,DispErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot(0.1*ImgNo,DispErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,DispErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot(0.1*ImgNo,DispErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,DispErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('$x-$direction translation (voxel)','interpreter','latex'); 
ylabel('Mean of $u$ error (voxel)','interpreter','latex');

h = legend({'Local DVC WS=30, ST=30','AL-DVC \quad WS=30, ST=30',...
    'Local DVC WS=20, ST=20','AL-DVC \quad WS=20, ST=20',...
    'Local DVC WS=10, ST=10','AL-DVC \quad WS=10, ST=10'},'FontSize',12,'location','northeastoutside','interpreter','latex');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';%axis([0,1,0.045,0.07]); 
axis([0,1,-0.022,0.022]); 
set(gca,'fontsize',20); set(gcf,'color','w');
% set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-png','-r300','painters','Fig_trans_disp_err2mean');

%% std
DispErrLocalICGN_ws30_st30 = [0.0058101    0.0068396    0.0066543    0.0063458    0.0062056    0.0062158    0.0060477 0.0061804    0.0063385    0.0066481    0.0072554];
DispErrAL_ws30_st30 = [0.005221    0.0061503    0.0059611    0.0056585    0.0054593    0.0055406    0.0054297 0.0056196    0.0055879    0.0059757    0.0065243];
DispErrLocalICGN_ws20_st20 = [0.010856     0.012323     0.012104     0.011417     0.011386     0.011321     0.011402 0.011532     0.011952      0.01231     0.013113];
 DispErrAL_ws20_st20 = [ 0.0096497     0.010874     0.010667     0.010129    0.0099678      0.00995    0.0099392 0.010188     0.010396     0.010752     0.011717];
DispErrLocalICGN_ws10_st10 = [0.032103     0.036753     0.035577     0.034679     0.034583     0.034135 0.034124     0.034749     0.035591     0.036756     0.037529];
DispErrAL_ws10_st10 = [0.027248     0.032063     0.030936     0.029765     0.029458      0.02932 0.029462     0.029924     0.030635     0.031849     0.033394];

figure(3); clf; %plot(ImgNo,DispErrFFT_ws10_st10,'d--');hold on
plot(0.1*ImgNo,DispErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,DispErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot(0.1*ImgNo,DispErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,DispErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot(0.1*ImgNo,DispErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,DispErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('$x-$direction translation (voxel)','interpreter','latex'); 
ylabel('Std of $u$ error (voxel)','interpreter','latex');

h = legend({'Local DVC WS=30, ST=30','AL-DVC \quad WS=30, ST=30',...
    'Local DVC WS=20, ST=20','AL-DVC \quad WS=20, ST=20',...
    'Local DVC WS=10, ST=10','AL-DVC \quad WS=10, ST=10'},'FontSize',12,'location','northeastoutside','interpreter','latex');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';%axis([0,1,0.045,0.07]); 
axis([0,1,0.003,0.042]); 
set(gca,'fontsize',20); set(gcf,'color','w');
% set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-png','-r300','painters','Fig_trans_disp_err3std');

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Strain RMS
StrainErrLocalICGN_ws30_st30 = [0.0019942    0.0024106    0.0024039    0.0024109    0.0024086    0.0024213     0.002429 0.0023973    0.0024032    0.0023962    0.0024379];
StrainErrAL_ws30_st30 = [0.00028474    0.0003592   0.00034707   0.00035365   0.00034689   0.00036142   0.00034847 0.00034736   0.00034295    0.0003472   0.00035686];
StrainErrLocalICGN_ws20_st20 = [0.0053509    0.0064622    0.0064627     0.006468    0.0065026    0.0065326    0.0065258 0.0065034    0.0064365    0.0064499    0.0065259];
StrainErrAL_ws20_st20 = [0.0007091   0.00086005   0.00086985   0.00085656   0.00086312   0.00086533   0.00086036 0.00085028   0.00085611   0.00085449   0.00087602];
StrainErrLocalICGN_ws10_st10 = [0.031761     0.039675     0.039887     0.040274     0.040465      0.04063 0.040593     0.040222     0.039754     0.039637     0.039746];
StrainErrAL_ws10_st10 = [0.0036584    0.0044674    0.0044594    0.0044586    0.0044906    0.0045041 0.0044851    0.0044607    0.0044547    0.0044447    0.0044864];

figure(4); clf; %plot(ImgNo,StrainErrFFT_ws10_st10,'d--');hold on
plot(0.1*ImgNo,StrainErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,StrainErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot(0.1*ImgNo,StrainErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,StrainErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot(0.1*ImgNo,StrainErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,StrainErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('$x-$direction translation (voxel)','interpreter','latex'); 
ylabel('Strain RMS error','interpreter','latex');

h = legend({'Local DVC WS=30, ST=30','AL-DVC \quad WS=30, ST=30',...
    'Local DVC WS=20, ST=20','AL-DVC \quad WS=20, ST=20',...
    'Local DVC WS=10, ST=10','AL-DVC \quad WS=10, ST=10'},'FontSize',12,'location','northeastoutside','interpreter','latex');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';
axis([0,1,8e-5,1.3e-1]);  
set(gca,'fontsize',20); set(gcf,'color','w');
set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-png','-r300','painters','Fig_trans_strain_err1rms');

%% Strain mean
StrainErrLocalICGN_ws30_st30 = [1.0351e-05  -1.3745e-05   1.4246e-05   5.7605e-06  -1.5883e-05 2.47e-06 2.3388e-06 -4.6795e-06  -1.2682e-05  -1.4255e-05  -1.3741e-05];
StrainErrAL_ws30_st30 = [-2.7135e-06   2.5872e-06  -2.6567e-06   3.1276e-06  -3.7164e-06  -1.9707e-07  -2.2426e-06 8.7619e-07   1.9829e-06  -2.5395e-07   2.5221e-07];
StrainErrLocalICGN_ws20_st20 = [4.9926e-06   3.6595e-05    2.405e-05   5.1926e-05   5.8973e-05   8.4192e-05   4.5266e-05 9.0796e-05   2.7809e-05   2.6734e-05   9.9523e-05];
StrainErrAL_ws20_st20 = [5.0318e-06   4.8992e-07  -2.1921e-07   3.9119e-06   1.2038e-06   3.0681e-06   2.0233e-06 1.2133e-06   1.4215e-06   4.5821e-06  -3.1588e-06];
StrainErrLocalICGN_ws10_st10 = [0.00031911   0.00057262    0.0007886    0.0010564    0.0012434     0.001298 0.0012359    0.0010464   0.00073633   0.00047373    0.0003887];
StrainErrAL_ws10_st10 = [4.9205e-06   3.0969e-06   3.6808e-06   7.5929e-06  -6.0549e-09   6.1017e-06 -5.0587e-08   3.2226e-06  -7.4389e-07   1.0191e-05   2.3172e-06];

figure(5); clf; %plot(ImgNo,StrainErrFFT_ws10_st10,'d--');hold on
plot(0.1*ImgNo,abs(StrainErrLocalICGN_ws30_st30),'o-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,abs(StrainErrAL_ws30_st30),'x-','markersize',7,'linewidth',1.5);
hold on;
plot(0.1*ImgNo,abs(StrainErrLocalICGN_ws20_st20),'s-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,abs(StrainErrAL_ws20_st20),'+-','markersize',7,'linewidth',1.5);
hold on;
plot(0.1*ImgNo,abs(StrainErrLocalICGN_ws10_st10),'^-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,abs(StrainErrAL_ws10_st10),'*-','markersize',7,'linewidth',1.5);

xlabel('$x-$direction translation (voxel)','interpreter','latex'); 
ylabel('Absolute value of mean of $e_{xx}$ error','interpreter','latex');

h = legend({'Local DVC WS=30, ST=30','AL-DVC \quad WS=30, ST=30',...
    'Local DVC WS=20, ST=20','AL-DVC \quad WS=20, ST=20',...
    'Local DVC WS=10, ST=10','AL-DVC \quad WS=10, ST=10'},'FontSize',12,'location','northeastoutside','interpreter','latex');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';%axis([0,1,0.045,0.07]); 
axis([0,1,1e-9,2e-2]); 
set(gca,'fontsize',20); set(gcf,'color','w');
set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-png','-r300','painters','Fig_trans_strain_err2mean');

%% std
StrainErrLocalICGN_ws30_st30 = [0.00066392   0.00074156   0.00072988   0.00070276   0.00068411   0.00068313   0.00068929 0.00072391   0.00072388   0.00075184   0.00083514];
StrainErrAL_ws30_st30 = [8.8388e-05   9.8632e-05   9.8675e-05   9.0059e-05   8.8114e-05   9.0739e-05   8.8803e-05 9.6324e-05    9.075e-05   9.7279e-05   0.00010986];
StrainErrLocalICGN_ws20_st20 = [0.0018178    0.0020599    0.0019964    0.0018958    0.0018734    0.0018518      0.00186 0.0019232    0.0019515    0.0020379    0.0021843];
StrainErrAL_ws20_st20 = [0.00023218   0.00025269    0.0002536   0.00023645    0.0002347   0.00023419   0.00022747  0.00023868   0.00024084   0.00024647   0.00027292];
StrainErrLocalICGN_ws10_st10 = [0.010388     0.012905     0.012583     0.012462     0.012411     0.012386 0.012335     0.012317     0.012447     0.012791     0.013085];
StrainErrAL_ws10_st10 = [ 0.0011765    0.0013782    0.0013335    0.0012879    0.0012837    0.0012697 0.0012694    0.0012917    0.0013271    0.0013638    0.0014419];

figure(6); clf; %plot(ImgNo,StrainErrFFT_ws10_st10,'d--');hold on
plot(0.1*ImgNo,StrainErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,StrainErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot(0.1*ImgNo,StrainErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,StrainErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot(0.1*ImgNo,StrainErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot(0.1*ImgNo,StrainErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('$x-$direction translation (voxel)','interpreter','latex'); 
ylabel('Std of $e_{xx}$ error','interpreter','latex');

h = legend({'Local DVC WS=30, ST=30','AL-DVC \quad WS=30, ST=30',...
    'Local DVC WS=20, ST=20','AL-DVC \quad WS=20, ST=20',...
    'Local DVC WS=10, ST=10','AL-DVC \quad WS=10, ST=10'},'FontSize',12,'location','northeastoutside','interpreter','latex');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';%axis([0,1,0.045,0.07]); 
axis([0,1,5e-5,0.03]); 
set(gca,'fontsize',20); set(gcf,'color','w');
set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-png','-r300','painters','Fig_trans_strain_err3std');

 
 