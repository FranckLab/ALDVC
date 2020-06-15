
close all;
ImgNo = 0:1:10;

DispErrFFT = [7.2465e-3
0.027257
0.050454
0.07039
0.085725
0.093704
0.092143
0.079653
0.076149
0.079483
0.087878
];
DispErrLocalICGN = [9.9846e-3
0.018268
0.021438
0.019801
0.015395
0.012623
0.015115
0.019555
0.021382
0.01806
0.012387
];
DispErrAL = [9.8634e-3
0.018206
0.021388
0.019732
0.015239
0.012402
0.014919
0.019483
0.021313
0.018009
0.012209
];

figure(1);
plot(ImgNo,DispErrFFT,'d--');hold on
plot(ImgNo,DispErrLocalICGN,'o-');hold on
plot(ImgNo,DispErrAL,'s:');

xlabel('Image No','interpreter','latex'); ylabel('Displacement RMS error (voxel)','interpreter','latex');
set(gca,'fontsize',20);
h = legend('FFT cross correlation','Local Subset DIC','AL-DIC','location','southeast');
axis 'auto x' 
% set(gca,'yscale','log');
% print('Fig_Sample1_disp','-dpng','-r600')


% StrainErrLocal = [2.5353e-10;5.3581e-4;6.1097e-4;6.0613e-4;5.7968e-4;5.9336e-4;6.6720e-4;7.8526e-4;9.0707e-4;1.1443e-3;0.0015197];
% StrainErrLocal = [3.4192e-12; 6.5791e-6; 7.0069e-6; 7.6150e-6; 6.2491e-6; 5.1284e-6; 5.5628e-6; 6.8593e-6; 6.9496e-6; 8.6889e-6; 1.1862e-5];
% figure(2);
% plot(ImgNo,StrainErrLocal,'ko-')
% 
% hold on
% % StrainErrAL = [7.4806e-7;1.4470e-5;7.2775e-6; 7.9931e-6; 7.3891e-6; 6.0268e-6; 1.0383e-5; 1.2120e-5; 1.3979e-5; 1.70399e-5; 2.0884e-5;];
% StrainErrAL = [2.8411e-9;5.8436e-6;6.2465e-6;6.2724e-6;5.2046e-6;4.2727e-6;4.9722e-6;5.2532e-6;5.4865e-6;7.0949e-6;8.6825e-6];
% 
% plot(ImgNo,StrainErrAL,'r+--')
% xlabel('Image No'); ylabel('Strain deviation error');
% set(gca,'fontsize',20);
% h = legend('Local DIC','Aug-Lag DIC','location','southeast');
% % set(gca,'yscale','log');
% print('Fig_Sample1_strain','-dpng','-r600')
 