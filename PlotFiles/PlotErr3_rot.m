
close all;
ImgNo = 0:1:5;
%% RMS
% ------ same size search using fft xcorr3 ------
% DispErrLocalICGN_ws30_st30 = [0.01299     0.017077     0.021456      0.02809     0.034026     0.032528];
% DispErrAL_ws30_st30 = [0.011752     0.012946     0.013403     0.013881     0.015072      0.01575];
% DispErrLocalICGN_ws20_st20 = [ 0.022829     0.043601     0.068673     0.065108     0.071871     0.073824];
% DispErrAL_ws20_st20 = [0.021653     0.023386     0.024005     0.024952     0.028388     0.030931];
% DispErrLocalICGN_ws10_st10 = [6.3703       14.029       25.782        30.84       40.993       46.134];
% DispErrAL_ws10_st10 = [0.062745      0.11958     0.068875      0.17289       4.8915       5.0274];
% ------ big search w/ xncorr ------
DispErrLocalICGN_ws30_st30 = [0.0099845 0.012142     0.012005     0.011725     0.011744     0.012935];
DispErrAL_ws30_st30 = [0.0090062 0.011646     0.010982     0.010827     0.010534     0.011716];
DispErrLocalICGN_ws20_st20 = [0.018643 0.021798      0.02168     0.021781     0.021603     0.021001 ];
DispErrAL_ws20_st20 = [0.016564 0.01922     0.019256     0.019345     0.019658     0.019684];
DispErrLocalICGN_ws10_st10 = [0.055704 0.064258     0.064769       0.0645     0.064206     0.063924];
DispErrAL_ws10_st10        = [0.047324 0.054206     0.053768     0.053733     0.053313     0.054094 ];

figure(1); clf; %plot(ImgNo,DispErrFFT_ws10_st10,'d--');hold on
plot( 5*ImgNo,DispErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,DispErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot( 5*ImgNo,DispErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,DispErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot( 5*ImgNo,DispErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,DispErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('Rotation angle in the $xy$ plane(deg)','interpreter','latex'); 
ylabel('Displacement RMS error (voxel)','interpreter','latex');

h = legend({'Local DIC WS=30, ST=30','AL-DIC     WS=30, ST=30',...
    'Local DIC WS=20, ST=20','AL-DIC     WS=20, ST=20',...
    'Local DIC WS=10, ST=10','AL-DIC     WS=10, ST=10'},'FontSize',12,'location','northeastoutside');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';
axis([0,25,0,0.073]); 
set(gca,'fontsize',20); set(gcf,'color','w');
% set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-pdf','-r300','painters','Fig_rot_disp_err1rms');

%% mean
% DispErrLocalICGN_ws30_st30 = [-0.00031659   -0.0001971   -0.0022245    0.0030629   -0.0052172    0.0013578];
% DispErrAL_ws30_st30 = [-0.00035228   0.00073567   4.9183e-08    0.0009025  -0.00030585  -0.00036775];
% DispErrLocalICGN_ws20_st20 = [-0.00043607   -0.0073912    0.0025981     0.002483    -0.010216   -0.0024789];
% DispErrAL_ws20_st20 = [-0.00034655    0.0012757  -0.00014506    0.0008963   -0.0013087   -0.0012742];
% DispErrLocalICGN_ws10_st10 = [0.78669        2.229       4.0176      0.28621      0.25789      0.78797];
% DispErrAL_ws10_st10 = [0.00016192    0.0024331   0.00038423   0.00040537      0.92335     -0.25535];
% ------ big search w/ xncorr ------
DispErrLocalICGN_ws30_st30 = [.00011545 0.00014177   7.7495e-05  -7.7051e-05   3.7759e-05    0.0015118];
DispErrAL_ws30_st30 = [.00010739 0.00019341   8.8787e-05  -0.00012629   2.4075e-06  -0.00019631 ];
DispErrLocalICGN_ws20_st20 = [.00012886 0.00019471   0.00025798   -0.0002342  -8.1505e-05   0.00030725];
DispErrAL_ws20_st20 = [.00013909 0.0002077   0.00025633  -0.00022427  -6.0294e-05   0.00014204 ];
DispErrLocalICGN_ws10_st10 = [7.3217e-05 0.00017266   0.00041092  -0.00046676  -0.00016369   8.7772e-05];
DispErrAL_ws10_st10 = [0.00011738 0.00017144   0.00032466  -0.00015176  -5.9103e-05    1.418e-05];

figure(2); clf; %plot(ImgNo,DispErrFFT_ws10_st10,'d--');hold on
plot(5*ImgNo,DispErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot(5*ImgNo,DispErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot(5*ImgNo,DispErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot(5*ImgNo,DispErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot(5*ImgNo,DispErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot(5*ImgNo,DispErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('Rotation angle in the $xy$ plane(deg)','interpreter','latex');
ylabel('Mean of $u$ error (voxel)','interpreter','latex');

h = legend({'Local DIC WS=30, ST=30','AL-DIC     WS=30, ST=30',...
    'Local DIC WS=20, ST=20','AL-DIC     WS=20, ST=20',...
    'Local DIC WS=10, ST=10','AL-DIC     WS=10, ST=10'},'FontSize',12,'location','northeastoutside');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';%axis([0,1,0.045,0.07]); 
axis([0,20,-6e-4,6.3e-4]);
set(gca,'fontsize',20); set(gcf,'color','w');
% set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-pdf','-r300','painters','Fig_rot_disp_err2mean');

%% std
% DispErrLocalICGN_ws30_st30 = [0.0073527    0.0084954    0.0094653     0.017177     0.021744     0.020359];
% DispErrAL_ws30_st30 = [0.0067055    0.0071844    0.0077308    0.0082233    0.0093267     0.010183];
% DispErrLocalICGN_ws20_st20 = [0.012144     0.025184     0.034316     0.037046     0.039422     0.045233];
% DispErrAL_ws20_st20 = [0.011775     0.013271     0.014214     0.015355     0.018022     0.021054];
% DispErrLocalICGN_ws10_st10 = [6.0511        13.11       24.749       30.254       40.598       45.793];
% DispErrAL_ws10_st10 = [0.034663       0.1027     0.041819      0.15848       4.7801       4.9329];
% ------ big search w/ xncorr ------
DispErrLocalICGN_ws30_st30 = [.0058101 0.0068549    0.0065498    0.0064418    0.0064699    0.0078064];
DispErrAL_ws30_st30 = [.005221 0.0065315    0.0059901    0.0059058    0.0058108    0.0065262];
DispErrLocalICGN_ws20_st20 = [.010856 0.011839     0.011902     0.011659     0.011994     0.011737];
DispErrAL_ws20_st20 = [.0096497 0.010409     0.010538     0.010433     0.010927     0.011014];
DispErrLocalICGN_ws10_st10 = [.032103 0.036015     0.036125     0.036379     0.036341     0.036764];
DispErrAL_ws10_st10 = [.027248 0.029937     0.029791     0.029619     0.029802     0.030431];

figure(3); clf; %plot(ImgNo,DispErrFFT_ws10_st10,'d--');hold on
plot( 5*ImgNo,DispErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,DispErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot( 5*ImgNo,DispErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,DispErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot( 5*ImgNo,DispErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,DispErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('Rotation angle in the $xy$ plane(deg)','interpreter','latex');
ylabel('Std of $u$ error (voxel)','interpreter','latex');

h = legend({'Local DIC WS=30, ST=30','AL-DIC     WS=30, ST=30',...
    'Local DIC WS=20, ST=20','AL-DIC     WS=20, ST=20',...
    'Local DIC WS=10, ST=10','AL-DIC     WS=10, ST=10'},'FontSize',12,'location','northeastoutside');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';%axis([0,1,0.045,0.07]); 
axis([0,25,0,0.042]); 
set(gca,'fontsize',20); set(gcf,'color','w');
% set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-pdf','-r300','painters','Fig_rot_disp_err3std');

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Strain RMS
% StrainErrLocalICGN_ws30_st30 = [0.14149      0.28278       0.4242       0.5655      0.70742      0.84874];
% StrainErrAL_ws30_st30 = [0.00036082   0.00038745   0.00038951   0.00040734   0.00044064    0.0004478];
% StrainErrLocalICGN_ws20_st20 = [0.14158      0.28224       0.4259       0.5674      0.70571      0.84791];
% StrainErrAL_ws20_st20 = [0.00088945   0.00093822    0.0009575   0.00098424    0.0010889    0.0011582];
% StrainErrLocalICGN_ws10_st10 = [0.58095       1.4597       1.8638       1.9409       2.4718       2.0528];
% StrainErrAL_ws10_st10 = [0.0046774    0.0081929    0.0049094     0.013064      0.15204      0.16563];
% ------ big search w/ xncorr ------
StrainErrLocalICGN_ws30_st30 = [.0019942 0.0023324    0.0022878      0.00229    0.0022621    0.0023129];
StrainErrAL_ws30_st30 = [.00028474 0.00037849   0.00036837   0.00036466   0.00036096   0.00042431];
StrainErrLocalICGN_ws20_st20 = [.0053509 0.0063334    0.0062707    0.0062819    0.0062608    0.0061894];
StrainErrAL_ws20_st20 = [.0007091 0.00085314   0.00084795      0.00086    0.0008854   0.00089093];
StrainErrLocalICGN_ws10_st10 = [.031761 0.039231     0.039702     0.040023     0.040821     0.051132];
StrainErrAL_ws10_st10 = [.0036584 0.0042371    0.0041932      0.00423    0.0041989    0.0043016];

figure(4); clf; %plot(ImgNo,StrainErrFFT_ws10_st10,'d--');hold on
plot( 5*ImgNo,StrainErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,StrainErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot( 5*ImgNo,StrainErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,StrainErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot( 5*ImgNo,StrainErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,StrainErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('Rotation angle in the $xy$ plane(deg)','interpreter','latex');
ylabel('Strain RMS error','interpreter','latex');

h = legend({'Local DIC WS=30, ST=30','AL-DIC     WS=30, ST=30',...
    'Local DIC WS=20, ST=20','AL-DIC     WS=20, ST=20',...
    'Local DIC WS=10, ST=10','AL-DIC     WS=10, ST=10'},'FontSize',12,'location','northeastoutside');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';
axis([0,25,1e-4,0.13]);  
set(gca,'fontsize',20); set(gcf,'color','w');
set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-pdf','-r300','painters','Fig_rot_strain_err1rms');

%% Strain mean
% StrainErrLocalICGN_ws30_st30 = [4.8402e-05   -0.0001174  -0.00069559   -0.0005067   0.00086679   -0.0010148];
% StrainErrAL_ws30_st30 = [-8.2197e-06   6.5396e-06  -5.6347e-06  -8.2078e-06   -1.492e-06   1.4998e-05];
% StrainErrLocalICGN_ws20_st20 = [0.0010775   0.00038714    0.0025873  -0.00056085  -4.0516e-05  -0.00023556];
% StrainErrAL_ws20_st20 = [1.0137e-06   1.6869e-05   2.0501e-05   1.1543e-05  -9.3138e-06    -2.57e-05];
% StrainErrLocalICGN_ws10_st10 = [-0.0017587    -0.032594    -0.046503     -0.16879     -0.59166      -0.2448];
% StrainErrAL_ws10_st10 = [-0.00012069  -7.7477e-06  -3.6204e-06   2.1798e-05    -0.022663    -0.014077];
% ------ big search w/ xncorr ------
StrainErrLocalICGN_ws30_st30 = [1.0351e-05 -0.0038053    -0.015202     -0.03411    -0.060293    -0.093684];
StrainErrAL_ws30_st30 = [-2.7135e-06 -0.0038063    -0.015195    -0.034066     -0.06031    -0.093688];
StrainErrLocalICGN_ws20_st20 = [4.9926e-06 -0.0037315     -0.01508    -0.033903    -0.060223    -0.093602];
StrainErrAL_ws20_st20 = [5.0318e-06 -0.0038104    -0.015195    -0.034074     -0.06029    -0.093681];
StrainErrLocalICGN_ws10_st10 = [0.00031911 -0.0026562    -0.013717    -0.032462    -0.058609    -0.091453 ];
StrainErrAL_ws10_st10 = [4.9205e-06 -0.0038053    -0.015171    -0.034039    -0.060237    -0.093614];

figure(5); clf; %plot(ImgNo,StrainErrFFT_ws10_st10,'d--');hold on
plot( 5*ImgNo,abs(StrainErrLocalICGN_ws30_st30),'o-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,abs(StrainErrAL_ws30_st30),'x-','markersize',7,'linewidth',1.5);
hold on;
plot( 5*ImgNo,abs(StrainErrLocalICGN_ws20_st20),'s-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,abs(StrainErrAL_ws20_st20),'+-','markersize',7,'linewidth',1.5);
hold on;
plot( 5*ImgNo,abs(StrainErrLocalICGN_ws10_st10),'^-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,abs(StrainErrAL_ws10_st10),'*-','markersize',7,'linewidth',1.5);

xlabel('Rotation angle in the $xy$ plane(deg)','interpreter','latex');
ylabel('Absolute value of mean of $e_{xx}$ error','interpreter','latex');

h = legend({'Local DIC WS=30, ST=30','AL-DIC     WS=30, ST=30',...
    'Local DIC WS=20, ST=20','AL-DIC     WS=20, ST=20',...
    'Local DIC WS=10, ST=10','AL-DIC     WS=10, ST=10'},'FontSize',12,'location','northeastoutside');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';%axis([0,1,0.045,0.07]); 
axis([0,25,1e-6,0.13]);  
set(gca,'fontsize',20); set(gcf,'color','w');
set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-pdf','-r300','painters','Fig_rot_strain_err2mean');

%% std
% StrainErrLocalICGN_ws30_st30 = [0.00077715    0.0011579     0.001591    0.0016513    0.0025391    0.0023299];
% StrainErrAL_ws30_st30 = [0.00011669   0.00011528   0.00012094   0.00012957   0.00014561   0.00016015];
% StrainErrLocalICGN_ws20_st20 = [0.0020749    0.0044117    0.0083735    0.0068462    0.0090796     0.010305];
% StrainErrAL_ws20_st20 = [0.00025439   0.00028278   0.00029343   0.00031847   0.00036994   0.00044124];
% StrainErrLocalICGN_ws10_st10 = [0.25267      0.58321      0.59294      0.52161      0.83655      0.60912];
% StrainErrAL_ws10_st10 = [0.0014379    0.0039615    0.0016198    0.0068887      0.10508      0.08044];
% ------ big search w/ xncorr ------
StrainErrLocalICGN_ws30_st30 = [.00066392 0.00073029   0.00075031   0.00076806   0.00069644   0.00072287];
StrainErrAL_ws30_st30 = [8.8388e-05 0.00011132   0.00010818   0.00010623   0.00011383   0.00011362 ];
StrainErrLocalICGN_ws20_st20 = [.0018178 0.0019965     0.002019    0.0020258    0.0019823    0.0020113   ];
StrainErrAL_ws20_st20 = [.00023218 0.00024597   0.00025434   0.00025162   0.00027103   0.00028581];
StrainErrLocalICGN_ws10_st10 = [.010388 0.012594     0.013102     0.013196     0.013493     0.017359];
StrainErrAL_ws10_st10 = [.0011765 0.0013023    0.0012945    0.0013232    0.0013331     0.001383];

figure(6); clf; %plot(ImgNo,StrainErrFFT_ws10_st10,'d--');hold on
plot(5*ImgNo,StrainErrLocalICGN_ws30_st30,'o-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,StrainErrAL_ws30_st30,'x-','markersize',7,'linewidth',1.5);
hold on;
plot( 5*ImgNo,StrainErrLocalICGN_ws20_st20,'s-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,StrainErrAL_ws20_st20,'+-','markersize',7,'linewidth',1.5);
hold on;
plot( 5*ImgNo,StrainErrLocalICGN_ws10_st10,'^-.','markersize',7,'linewidth',1.5);hold on;
plot( 5*ImgNo,StrainErrAL_ws10_st10,'*-','markersize',7,'linewidth',1.5);

xlabel('Rotation angle in the $xy$ plane(deg)','interpreter','latex'); 
ylabel('Std of $e_{xx}$ error','interpreter','latex');

h = legend({'Local DIC WS=30, ST=30','AL-DIC     WS=30, ST=30',...
    'Local DIC WS=20, ST=20','AL-DIC     WS=20, ST=20',...
    'Local DIC WS=10, ST=10','AL-DIC     WS=10, ST=10'},'FontSize',12,'location','northeastoutside');
axis 'auto x' ;  a = gca; a.TickLabelInterpreter = 'latex';%axis([0,1,0.045,0.07]); 
axis([0,25,5e-5,0.04]); 
set(gca,'fontsize',20); set(gcf,'color','w');
set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')
export_fig(gcf,'-pdf','-r300','painters','Fig_rot_strain_err3std');

 
 