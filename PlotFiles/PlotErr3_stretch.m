
close all;
ImgNo = 0 :1:6;
%% RMS
% ------ same size search using fft xcorr3 ------
% DispErrLocalICGN_ws30_st30 = [0.01299     0.017077     0.021456      0.02809     0.034026     0.032528];
% DispErrAL_ws30_st30 = [0.011752     0.012946     0.013403     0.013881     0.015072      0.01575];
% DispErrLocalICGN_ws20_st20 = [ 0.022829     0.043601     0.068673     0.065108     0.071871     0.073824];
% DispErrAL_ws20_st20 = [0.021653     0.023386     0.024005     0.024952     0.028388     0.030931];
% DispErrLocalICGN_ws10_st10 = [6.3703       14.029       25.782        30.84       40.993       46.134];
% DispErrAL_ws10_st10 = [0.062745      0.11958     0.068875      0.17289       4.8915       5.0274];
% ------ big search w/ xncorr ------
DispErrLocalICGN_ws30_st30 = [0.0099845 0.012643     0.013003     0.012837      0.01322     0.014477     0.015315];
DispErrAL_ws30_st30 = [0.0090062 0.011729     0.011942     0.011914      0.01199     0.013307     0.013946];
DispErrLocalICGN_ws20_st20 = [0.018643 0.022674     0.023582     0.024147     0.025434     0.026925     0.028888];
DispErrAL_ws20_st20 = [0.016564 0.020593     0.020978     0.021412      0.02246     0.024073     0.025947];
DispErrLocalICGN_ws10_st10 = [0.055704 0.065992 0.068004 0.069813 0.073227 0.077539 0.081599];
DispErrAL_ws10_st10        = [0.047324 0.057962 0.060349 0.063289 0.067777 0.074230 0.080539];

figure(1); clf; %plot(ImgNo,DispErrFFT_ws10_st10,'d--');hold on
plot(1+0.05*ImgNo,DispErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,DispErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot(1+0.05*ImgNo,DispErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,DispErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot(1+0.05*ImgNo,DispErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,DispErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('$x-$direction stretch ratio','interpreter','latex'); 
ylabel('Displacement RMS error (voxel)','interpreter','latex');

h = legend({'Local DVC WS=30, ST=30','AL-DVC \quad WS=30, ST=30',...
    'Local DVC WS=20, ST=20','AL-DVC \quad WS=20, ST=20',...
    'Local DVC WS=10, ST=10','AL-DVC \quad WS=10, ST=10'},'FontSize',12,'location','northeastoutside','interpreter','latex');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';
% axis([1,1.3,0.007,0.082]); 
set(gca,'fontsize',20); set(gcf,'color','w');
% set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-png','-r300','painters','Fig_unistretch_disp_err1rms');

%% mean
% DispErrLocalICGN_ws30_st30 = [-0.00031659   -0.0001971   -0.0022245    0.0030629   -0.0052172    0.0013578];
% DispErrAL_ws30_st30 = [-0.00035228   0.00073567   4.9183e-08    0.0009025  -0.00030585  -0.00036775];
% DispErrLocalICGN_ws20_st20 = [-0.00043607   -0.0073912    0.0025981     0.002483    -0.010216   -0.0024789];
% DispErrAL_ws20_st20 = [-0.00034655    0.0012757  -0.00014506    0.0008963   -0.0013087   -0.0012742];
% DispErrLocalICGN_ws10_st10 = [0.78669        2.229       4.0176      0.28621      0.25789      0.78797];
% DispErrAL_ws10_st10 = [0.00016192    0.0024331   0.00038423   0.00040537      0.92335     -0.25535];
% ------ big search w/ xncorr ------
DispErrLocalICGN_ws30_st30 = [0.00011545  -0.00016807  -0.00055878   9.1883e-05   0.00031357  -0.00030557   4.3151e-06];
DispErrAL_ws30_st30 = [0.00010739  -0.00017984  -0.00057613   0.00012662   0.00030192  -0.00039069   6.0943e-05];
DispErrLocalICGN_ws20_st20 = [0.00012886 0.00021252   0.00048228   0.00038167   0.00054881  -0.00051691  -0.00044425];
DispErrAL_ws20_st20 = [0.00013909  0.00059448   0.00052728   0.00034518   0.00051377   -0.0005066   -0.0005387];
DispErrLocalICGN_ws10_st10 = [7.3217e-05  -7.4260e-5 -3.9745e-4 1.2752e-4 6.6744e-4 -5.7726e-4 -2.0659e-4];
DispErrAL_ws10_st10 = [0.00011738  -2.1872e-5 -1.5294e-3 -1.2621e-4 9.0979e-4 -5.1568e-4 4.0414e-4];

figure(2); clf; %plot(ImgNo,DispErrFFT_ws10_st10,'d--');hold on
plot(1+0.05*ImgNo,DispErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,DispErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot(1+0.05*ImgNo,DispErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,DispErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot(1+0.05*ImgNo,DispErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,DispErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('$x-$direction stretch ratio','interpreter','latex');  
ylabel('Mean of $u$ error (voxel)','interpreter','latex');

h = legend({'Local DVC WS=30, ST=30','AL-DVC \quad WS=30, ST=30',...
    'Local DVC WS=20, ST=20','AL-DVC \quad WS=20, ST=20',...
    'Local DVC WS=10, ST=10','AL-DVC \quad WS=10, ST=10'},'FontSize',12,'location','northeastoutside','interpreter','latex');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';%axis([0,1,0.045,0.07]); 
axis([1,1.3,-0.005,0.005]);
set(gca,'fontsize',20); set(gcf,'color','w');
% set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-png','-r300','painters','Fig_unistretch_disp_err2mean');

%% std
% DispErrLocalICGN_ws30_st30 = [0.0073527    0.0084954    0.0094653     0.017177     0.021744     0.020359];
% DispErrAL_ws30_st30 = [0.0067055    0.0071844    0.0077308    0.0082233    0.0093267     0.010183];
% DispErrLocalICGN_ws20_st20 = [0.012144     0.025184     0.034316     0.037046     0.039422     0.045233];
% DispErrAL_ws20_st20 = [0.011775     0.013271     0.014214     0.015355     0.018022     0.021054];
% DispErrLocalICGN_ws10_st10 = [6.0511        13.11       24.749       30.254       40.598       45.793];
% DispErrAL_ws10_st10 = [0.034663       0.1027     0.041819      0.15848       4.7801       4.9329];
% ------ big search w/ xncorr ------
DispErrLocalICGN_ws30_st30 = [0.0058101  0.0069681    0.0071967    0.0074743    0.0077563    0.0090491    0.0098875];
DispErrAL_ws30_st30 = [0.005221 0.006495    0.0065558    0.0069592     0.006978     0.008338    0.0090409];
DispErrLocalICGN_ws20_st20 = [0.010856 0.012042     0.013378     0.014167      0.01539     0.017001     0.018851];
DispErrAL_ws20_st20 = [0.0096497  0.010887     0.011956     0.012536     0.013503     0.015097     0.016746];
DispErrLocalICGN_ws10_st10 = [0.032103 0.036952 0.039748 0.04251 0.046744 0.051781 0.0568];
DispErrAL_ws10_st10 = [0.027248  0.031989 0.034557 0.038139 0.04277 0.048831 0.054835];

figure(3); clf; %plot(ImgNo,DispErrFFT_ws10_st10,'d--');hold on
plot(1+0.05*ImgNo,DispErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,DispErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot(1+0.05*ImgNo,DispErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,DispErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot(1+0.05*ImgNo,DispErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,DispErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('$x-$direction stretch ratio','interpreter','latex'); 
ylabel('Std of $u$ error (voxel)','interpreter','latex');

h = legend({'Local DVC WS=30, ST=30','AL-DVC \quad WS=30, ST=30',...
    'Local DVC WS=20, ST=20','AL-DVC \quad WS=20, ST=20',...
    'Local DVC WS=10, ST=10','AL-DVC \quad WS=10, ST=10'},'FontSize',12,'location','northeastoutside','interpreter','latex');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';%axis([0,1,0.045,0.07]); 
axis([1,1.3,0.000,0.065]); 
set(gca,'fontsize',20); set(gcf,'color','w');
% set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-png','-r300','painters','Fig_unistretch_disp_err3std');

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Strain RMS
% StrainErrLocalICGN_ws30_st30 = [0.14149      0.28278       0.4242       0.5655      0.70742      0.84874];
% StrainErrAL_ws30_st30 = [0.00036082   0.00038745   0.00038951   0.00040734   0.00044064    0.0004478];
% StrainErrLocalICGN_ws20_st20 = [0.14158      0.28224       0.4259       0.5674      0.70571      0.84791];
% StrainErrAL_ws20_st20 = [0.00088945   0.00093822    0.0009575   0.00098424    0.0010889    0.0011582];
% StrainErrLocalICGN_ws10_st10 = [0.58095       1.4597       1.8638       1.9409       2.4718       2.0528];
% StrainErrAL_ws10_st10 = [0.0046774    0.0081929    0.0049094     0.013064      0.15204      0.16563];
% ------ big search w/ xncorr ------
StrainErrLocalICGN_ws30_st30 = [0.0019942  0.0024123    0.0024945    0.0026499    0.0027304    0.0029278    0.0031146];
StrainErrAL_ws30_st30 = [0.00028474  0.00037291   0.00038987   0.00038798   0.00039515   0.00043566   0.00047581];
StrainErrLocalICGN_ws20_st20 = [0.0053509  0.0067028    0.0068392     0.007172    0.0075879    0.0082861    0.0088324];
StrainErrAL_ws20_st20 = [0.0007091 0.00088462   0.00093153   0.00093736   0.00099327    0.0010647    0.0011381];
StrainErrLocalICGN_ws10_st10 = [0.031761 0.040157 0.042379 0.045775 0.049927 0.055536 0.060792];
StrainErrAL_ws10_st10 = [0.0036584  0.0044781 0.0046472 0.0048539 0.0051752 0.0056623 0.0061382];

figure(4); clf; %plot(ImgNo,StrainErrFFT_ws10_st10,'d--');hold on
plot(1+0.05*ImgNo,StrainErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,StrainErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot(1+0.05*ImgNo,StrainErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,StrainErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot(1+0.05*ImgNo,StrainErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,StrainErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('$x-$direction stretch ratio','interpreter','latex'); 
ylabel('Strain RMS error','interpreter','latex');

h = legend({'Local DVC WS=30, ST=30','AL-DVC \quad WS=30, ST=30',...
    'Local DVC WS=20, ST=20','AL-DVC \quad WS=20, ST=20',...
    'Local DVC WS=10, ST=10','AL-DVC \quad WS=10, ST=10'},'FontSize',12,'location','northeastoutside','interpreter','latex');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';
axis([1,1.3,1e-4,0.2]);  
set(gca,'fontsize',20); set(gcf,'color','w');
set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-png','-r300','painters','Fig_unistretch_strain_err1rms');

%% Strain mean
% StrainErrLocalICGN_ws30_st30 = [4.8402e-05   -0.0001174  -0.00069559   -0.0005067   0.00086679   -0.0010148];
% StrainErrAL_ws30_st30 = [-8.2197e-06   6.5396e-06  -5.6347e-06  -8.2078e-06   -1.492e-06   1.4998e-05];
% StrainErrLocalICGN_ws20_st20 = [0.0010775   0.00038714    0.0025873  -0.00056085  -4.0516e-05  -0.00023556];
% StrainErrAL_ws20_st20 = [1.0137e-06   1.6869e-05   2.0501e-05   1.1543e-05  -9.3138e-06    -2.57e-05];
% StrainErrLocalICGN_ws10_st10 = [-0.0017587    -0.032594    -0.046503     -0.16879     -0.59166      -0.2448];
% StrainErrAL_ws10_st10 = [-0.00012069  -7.7477e-06  -3.6204e-06   2.1798e-05    -0.022663    -0.014077];
% ------ big search w/ xncorr ------
StrainErrLocalICGN_ws30_st30 = [1.0351e-05  -1.1512e-05   8.8192e-05   -0.0001831   3.3174e-05   2.2809e-06  -0.00024189];
StrainErrAL_ws30_st30 = [-2.7135e-06  -3.1562e-08  -3.0368e-06  -5.5087e-06  -6.9692e-06  -1.4787e-06   8.0106e-06];
StrainErrLocalICGN_ws20_st20 = [4.9926e-06 -0.0011968  -0.00026954  -6.7928e-05   4.2119e-05  -0.00014215  -0.00023082];
StrainErrAL_ws20_st20 = [5.0318e-06  8.3906e-06   4.8041e-06   2.6109e-06  -1.3434e-05  -7.3495e-06   1.5935e-05];
StrainErrLocalICGN_ws10_st10 = [0.00031911 -0.0011016 -0.0024241 -0.004334 -0.0057683 -0.0078925 -0.0092213];
StrainErrAL_ws10_st10 = [4.9205e-06 3.993e-5 1.6867e-6 -1.5905e-5 -2.6772e-6 -1.4316e-5 -1.0128e-5];

figure(5); clf; %plot(ImgNo,StrainErrFFT_ws10_st10,'d--');hold on
plot(1+0.05*ImgNo,abs(StrainErrLocalICGN_ws30_st30),'o-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,abs(StrainErrAL_ws30_st30),'x-','markersize',7,'linewidth',1.5);
hold on;
plot(1+0.05*ImgNo,abs(StrainErrLocalICGN_ws20_st20),'s-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,abs(StrainErrAL_ws20_st20),'+-','markersize',7,'linewidth',1.5);
hold on;
plot(1+0.05*ImgNo,abs(StrainErrLocalICGN_ws10_st10),'^-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,abs(StrainErrAL_ws10_st10),'*-','markersize',7,'linewidth',1.5);

xlabel('$x-$direction stretch ratio','interpreter','latex');  
ylabel('Absolute value of mean of $e_{xx}$ error','interpreter','latex');

h = legend({'Local DVC WS=30, ST=30','AL-DVC \quad WS=30, ST=30',...
    'Local DVC WS=20, ST=20','AL-DVC \quad WS=20, ST=20',...
    'Local DVC WS=10, ST=10','AL-DVC \quad WS=10, ST=10'},'FontSize',12,'location','northeastoutside','interpreter','latex');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';%axis([0,1,0.045,0.07]); 
axis([1,1.3,1e-9,0.2]); 
set(gca,'fontsize',20); set(gcf,'color','w');
set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-png','-r300','painters','Fig_unistretch_strain_err2mean');

%% std
% StrainErrLocalICGN_ws30_st30 = [0.00077715    0.0011579     0.001591    0.0016513    0.0025391    0.0023299];
% StrainErrAL_ws30_st30 = [0.00011669   0.00011528   0.00012094   0.00012957   0.00014561   0.00016015];
% StrainErrLocalICGN_ws20_st20 = [0.0020749    0.0044117    0.0083735    0.0068462    0.0090796     0.010305];
% StrainErrAL_ws20_st20 = [0.00025439   0.00028278   0.00029343   0.00031847   0.00036994   0.00044124];
% StrainErrLocalICGN_ws10_st10 = [0.25267      0.58321      0.59294      0.52161      0.83655      0.60912];
% StrainErrAL_ws10_st10 = [0.0014379    0.0039615    0.0016198    0.0068887      0.10508      0.08044];
% ------ big search w/ xncorr ------
StrainErrLocalICGN_ws30_st30 = [0.00066392 0.00078873   0.00082883   0.00097573    0.0010913    0.0011962    0.0013494];
StrainErrAL_ws30_st30 = [8.8388e-05  0.00010573   0.00011988   0.00012258   0.00011949   0.00015073   0.00016925];
StrainErrLocalICGN_ws20_st20 = [0.0018178  0.0022207    0.0023627     0.002691    0.0031602    0.0036567    0.0041098];
StrainErrAL_ws20_st20 = [0.00023218  0.00025342   0.00029016   0.00030241   0.00033481    0.0003746   0.00042896];
StrainErrLocalICGN_ws10_st10 = [0.010388  0.013208 0.015296 0.018212 0.02159 0.02597 0.029859];
StrainErrAL_ws10_st10 = [ 0.0011765  0.0013825 0.0014961 0.001665 0.0018847 0.0021717 0.0024331];

figure(6); clf; %plot(ImgNo,StrainErrFFT_ws10_st10,'d--');hold on
plot(1+0.05*ImgNo,StrainErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,StrainErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot(1+0.05*ImgNo,StrainErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,StrainErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot(1+0.05*ImgNo,StrainErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot(1+0.05*ImgNo,StrainErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('$x-$direction stretch ratio','interpreter','latex'); 
ylabel('Std of $e_{xx}$ error','interpreter','latex');

h = legend({'Local DVC WS=30, ST=30','AL-DVC \quad WS=30, ST=30',...
    'Local DVC WS=20, ST=20','AL-DVC \quad WS=20, ST=20',...
    'Local DVC WS=10, ST=10','AL-DVC \quad WS=10, ST=10'},'FontSize',12,'location','northeastoutside','interpreter','latex');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';%axis([0,1,0.045,0.07]); 
axis([1,1.3,2e-5,1.5e-1]); 
set(gca,'fontsize',20); set(gcf,'color','w');
set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-png','-r300','painters','Fig_unistretch_strain_err3std');

 
 