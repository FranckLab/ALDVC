
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessing for Lauren's ind exp
% convergence in Subprblem 1 ADMM scheme
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 52; N = 52; L = 19; winstepsize = [8,8,8];
ImgSeqNum = 2; %:length(ResultDisp)+1
USubpb2 = ResultDisp{ImgSeqNum-1}.U;
ULocalICGN = ResultDisp{ImgSeqNum-1}.ULocalICGN;
U0 = ResultDisp{ImgSeqNum-1}.U0;

FSubpb2 = ResultDefGrad{ImgSeqNum-1}.F;
FLocalICGN = ResultDefGrad{ImgSeqNum-1}.FLocalICGN;

coordinatesFEM = ResultcoordinatesFEM{ImgSeqNum-1}.coordinatesFEM;
elementsFEM = ResultelementsFEM{ImgSeqNum-1}.elementsFEM;
ConvItPerEle = ResultConvItPerEle{ImgSeqNum-1}.ConvItPerEle;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ConvIterNumtemp = reshape(ConvItPerEle(:,1),M,N,L);
% slice 1
ConvIterNumSlice = ConvIterNumtemp(:,:,1); figure,  imagesc(ConvIterNumSlice); %surf(ConvIterNumSlice,'edgecolor','none'); 
view(2); axis equal; axis tight; colormap hot;  box on; set(gca,'fontsize',18); caxis([0,100])
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar;  b.TickLabelInterpreter = 'latex';
% slice 10
ConvIterNumSlice = ConvIterNumtemp(:,:,10); figure, imagesc(ConvIterNumSlice); %surf(ConvIterNumSlice,'edgecolor','none'); 
view(2); axis equal; axis tight; colormap hot;  box on; set(gca,'fontsize',18); caxis([0,100])
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar;  b.TickLabelInterpreter = 'latex';
% slice 19
ConvIterNumSlice = ConvIterNumtemp(:,:,19); figure,  imagesc(ConvIterNumSlice); %surf(ConvIterNumSlice,'edgecolor','none'); 
view(2); axis equal; axis tight; colormap hot;  box on; set(gca,'fontsize',18); caxis([0,100])
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar;  b.TickLabelInterpreter = 'latex';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ConvIterNumtemp = reshape(ConvItPerEle(:,6),M,N,L);
% slice 1
ConvIterNumSlice = ConvIterNumtemp(:,:,1); figure,  imagesc(ConvIterNumSlice); %surf(ConvIterNumSlice,'edgecolor','none'); 
view(2); axis equal; axis tight; colormap hot;  box on; set(gca,'fontsize',18); caxis([0,100])
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar;  b.TickLabelInterpreter = 'latex';
% slice 10
ConvIterNumSlice = ConvIterNumtemp(:,:,10); figure, imagesc(ConvIterNumSlice); %surf(ConvIterNumSlice,'edgecolor','none'); 
view(2); axis equal; axis tight; colormap hot;  box on; set(gca,'fontsize',18); caxis([0,100])
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar;  b.TickLabelInterpreter = 'latex';
% slice 19
ConvIterNumSlice = ConvIterNumtemp(:,:,19); figure,  imagesc(ConvIterNumSlice); %surf(ConvIterNumSlice,'edgecolor','none'); 
view(2); axis equal; axis tight; colormap hot;  box on; set(gca,'fontsize',18); caxis([0,100])
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar;  b.TickLabelInterpreter = 'latex';

 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============================================= 
utemp = reshape(USubpb2(1:3:end),M,N,L); vtemp = reshape(USubpb2(2:3:end),M,N,L); wtemp = reshape(USubpb2(3:3:end),M,N,L);

% [XXGrid,ZZGrid] = ndgrid(320:4:728,20:4:164); 
[XXGrid,ZZGrid] = ndgrid(coordinatesFEM(1,1):winstepsize(1):coordinatesFEM(end,1), ...
    coordinatesFEM(1,3):winstepsize(3):coordinatesFEM(end,3)); 
 XX2Grid = interp(XXGrid(:,1)',4); ZZ2Grid = interp(ZZGrid(1,:),4);
%XX2Grid = [320:4:728]; ZZ2Grid = [20:4:164]';

tempdata = squeeze(utemp(:,1+28,:));
dispu_y29 = gridfit(reshape(XXGrid,M*L,1),reshape(ZZGrid,M*L,1),reshape(tempdata,M*L,1),XX2Grid,ZZ2Grid);
figure, %surf(0.42*XX2Grid,0.425*ZZ2Grid,0.42*dispu_y29,'edgecolor','none');
contourf(0.42*XX2Grid,0.425*ZZ2Grid,0.42*dispu_y29,128,'linestyle','none');
axis equal; axis tight; view(2); colormap(cMap); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex'; title('dispu y29 (units:um)')
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
%if x(M,N) < 200,set(gca,'XTick',[]);end
%if y(M,N) < 200,set(gca,'YTick',[]);end

% tempdata = squeeze(vtemp(:,1+28,:));
% tempdata2 = gridfit(reshape(XXGrid,M*L,1),reshape(ZZGrid,M*L,1),reshape(tempdata,M*L,1),XX2Grid,ZZ2Grid);
% figure, surf(XX2Grid,ZZ2Grid,tempdata2,'edgecolor','none');
% axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex';

% tempdata = squeeze(wtemp(:,1+28,:));
% dispw_y29 = gridfit(reshape(XXGrid,M*L,1),reshape(ZZGrid,M*L,1),reshape(tempdata,M*L,1),XX2Grid,ZZ2Grid);
% figure, surf(XX2Grid,ZZ2Grid,dispw_y29,'edgecolor','none');
% axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex'; title('dispw y29')

% ============================================= 
utemp = reshape(USubpb2(1:3:end),M,N,L); vtemp = reshape(USubpb2(2:3:end),M,N,L); wtemp = reshape(USubpb2(3:3:end),M,N,L);

[YYGrid,ZZGrid] = ndgrid(320:8:728,20:8:164); 
%YY2Grid = interp(YYGrid(:,1)',4); ZZ2Grid = interp(ZZGrid(1,:),4);
YY2Grid = [320:1:728]; ZZ2Grid = [20:1:164]';

% tempdata = squeeze(utemp(1+23,:,:));
% tempdata2 = gridfit(reshape(YYGrid,N*L,1),reshape(ZZGrid,N*L,1),reshape(tempdata,N*L,1),YY2Grid,ZZ2Grid);
% figure, surf(YY2Grid,ZZ2Grid,tempdata2,'edgecolor','none');
% axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex';
% % xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% %if x(M,N) < 200,set(gca,'XTick',[]);end
% %if y(M,N) < 200,set(gca,'YTick',[]);end

tempdata = squeeze(vtemp(1+23,:,:));
dispu_x24 = gridfit(reshape(YYGrid,N*L,1),reshape(ZZGrid,N*L,1),reshape(tempdata,N*L,1),YY2Grid,ZZ2Grid);
figure, %surf(0.42*YY2Grid,0.425*ZZ2Grid,0.42*dispu_x24,'edgecolor','none');
contourf(0.42*YY2Grid,0.425*ZZ2Grid,0.42*dispu_x24,128,'linestyle','none');
axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex'; title('dispv x24  (units:um)')

tempdata = squeeze(wtemp(1+23,:,:));
dispw_x24 = gridfit(reshape(YYGrid,N*L,1),reshape(ZZGrid,N*L,1),reshape(tempdata,N*L,1),YY2Grid,ZZ2Grid);
figure, %surf(0.42*YY2Grid,0.425*ZZ2Grid,0.425*(dispw_x24+0.1335),'edgecolor','none');
contourf(0.42*YY2Grid,0.425*ZZ2Grid,0.425*(dispw_x24+0.1335),128,'linestyle','none');
axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex'; title('dispw x24 (units:um)')

% ============================================= 
% utemp = reshape(USubpb2(1:3:end),M,N,L); vtemp = reshape(USubpb2(2:3:end),M,N,L); wtemp = reshape(USubpb2(3:3:end),M,N,L);
% 
% [XXGrid,YYGrid] = ndgrid(320:8:728,320:8:728); 
% %XX2Grid = interp(XXGrid(:,1)',4); YY2Grid = interp(YYGrid(1,:)',4); 
% XX2Grid = [320:1:728]; YY2Grid = [320:1:728]';
% 
% tempdata = squeeze(utemp(:,:,10));
% dispu_z10 = gridfit(reshape(XXGrid,M*N,1),reshape(YYGrid,M*N,1),reshape(tempdata,M*N,1),XX2Grid,YY2Grid);
% figure, surf(XX2Grid,YY2Grid,dispu_z10,'edgecolor','none');
% axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex'; title('dispu z10')
% % xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% %if x(M,N) < 200,set(gca,'XTick',[]);end
% %if y(M,N) < 200,set(gca,'YTick',[]);end
% 
% tempdata = squeeze(vtemp(:,:,10));
% dispv_z10 = gridfit(reshape(XXGrid,M*N,1),reshape(YYGrid,M*N,1),reshape(tempdata,M*N,1),XX2Grid,YY2Grid);
% figure, surf(XX2Grid,YY2Grid,dispv_z10,'edgecolor','none');
% axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex'; title('dispv z10')
% 
% tempdata = squeeze(wtemp(:,:,10));
% dispw_z10 = gridfit(reshape(XXGrid,M*N,1),reshape(YYGrid,M*N,1),reshape(tempdata,M*N,1),XX2Grid,YY2Grid);
% figure, surf(XX2Grid,YY2Grid,dispw_z10,'edgecolor','none');
% axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex'; title('dispw z10')


%% Plot slice strain fields
%%
e12temp = reshape(0.5*(FStraintemp(2:9:end)+FStraintemp(4:9:end)),M-2*Rad,N-2*Rad,L-2*Rad);
e13temp = reshape(0.5*(FStraintemp(3:9:end)+FStraintemp(7:9:end)),M-2*Rad,N-2*Rad,L-2*Rad);
e23temp = reshape(0.5*(FStraintemp(6:9:end)+FStraintemp(8:9:end)),M-2*Rad,N-2*Rad,L-2*Rad);
e33temp = reshape(FStraintemp(9:9:end),M-2*Rad,N-2*Rad,L-2*Rad); 
e22temp = reshape(FStraintemp(5:9:end),M-2*Rad,N-2*Rad,L-2*Rad);  
e11temp = reshape(FStraintemp(1:9:end),M-2*Rad,N-2*Rad,L-2*Rad);  
[XXGrid,YYGrid] = ndgrid(320+Rad*winstepsize(1):winstepsize(1):728-Rad*winstepsize(1),...
    320+Rad*winstepsize(2):winstepsize(2):728-Rad*winstepsize(2)); 
XX2Grid = interp(XXGrid(:,1)',4); YY2Grid = interp(YYGrid(1,:)',4); 


tempdata = squeeze(e12temp(:,:,10-Rad));
tempdata2 = gridfit(reshape(XXGrid,(M-2*Rad)*(N-2*Rad),1),reshape(YYGrid,(M-2*Rad)*(N-2*Rad),1),reshape(tempdata,(M-2*Rad)*(N-2*Rad),1),XX2Grid,YY2Grid);
figure, surf(XX2Grid,YY2Grid,tempdata2,'edgecolor','none');
axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex'; 

tempdata = squeeze(e13temp(:,:,10-Rad));
tempdata2 = gridfit(reshape(XXGrid,(M-2*Rad)*(N-2*Rad),1),reshape(YYGrid,(M-2*Rad)*(N-2*Rad),1),reshape(tempdata,(M-2*Rad)*(N-2*Rad),1),XX2Grid,YY2Grid);
figure, surf(XX2Grid,YY2Grid,tempdata2,'edgecolor','none');
axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex'; caxis([-0.025,0.025])

tempdata = squeeze(e23temp(:,:,10-Rad));
tempdata2 = gridfit(reshape(XXGrid,(M-2*Rad)*(N-2*Rad),1),reshape(YYGrid,(M-2*Rad)*(N-2*Rad),1),reshape(tempdata,(M-2*Rad)*(N-2*Rad),1),XX2Grid,YY2Grid);
figure, surf(XX2Grid,YY2Grid,tempdata2,'edgecolor','none');
axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

tempdata = squeeze(e11temp(:,:,10-Rad));
tempdata2 = gridfit(reshape(XXGrid,(M-2*Rad)*(N-2*Rad),1),reshape(YYGrid,(M-2*Rad)*(N-2*Rad),1),reshape(tempdata,(M-2*Rad)*(N-2*Rad),1),XX2Grid,YY2Grid);
figure, surf(XX2Grid,YY2Grid,tempdata2,'edgecolor','none');
axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

tempdata = squeeze(e22temp(:,:,10-Rad));
tempdata2 = gridfit(reshape(XXGrid,(M-2*Rad)*(N-2*Rad),1),reshape(YYGrid,(M-2*Rad)*(N-2*Rad),1),reshape(tempdata,(M-2*Rad)*(N-2*Rad),1),XX2Grid,YY2Grid);
figure, surf(XX2Grid,YY2Grid,tempdata2,'edgecolor','none');
axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

tempdata = squeeze(e33temp(:,:,10-Rad));
tempdata2 = gridfit(reshape(XXGrid,(M-2*Rad)*(N-2*Rad),1),reshape(YYGrid,(M-2*Rad)*(N-2*Rad),1),reshape(tempdata,(M-2*Rad)*(N-2*Rad),1),XX2Grid,YY2Grid);
figure, surf(XX2Grid,YY2Grid,tempdata2,'edgecolor','none');
axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';  caxis([-0.05,0.000])

% ============================================= 
[XXGrid,ZZGrid] = ndgrid(320+Rad*8:8:728-Rad*8,20+Rad*8:8:164-Rad*8); 
XX2Grid = interp(XXGrid(:,1)',4); ZZ2Grid = interp(ZZGrid(1,:),4);

tempdata = squeeze(e33temp(:,1+28-Rad,:));
tempdata2 = gridfit(reshape(XXGrid,(M-2*Rad)*(L-2*Rad),1),reshape(ZZGrid,(M-2*Rad)*(L-2*Rad),1),reshape(tempdata,(M-2*Rad)*(L-2*Rad),1),XX2Grid,ZZ2Grid);
figure, surf(XX2Grid,ZZ2Grid,tempdata2,'edgecolor','none');
axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

tempdata = squeeze(e12temp(:,1+28-Rad,:));
tempdata2 = gridfit(reshape(XXGrid,(M-2*Rad)*(L-2*Rad),1),reshape(ZZGrid,(M-2*Rad)*(L-2*Rad),1),reshape(tempdata,(M-2*Rad)*(L-2*Rad),1),XX2Grid,ZZ2Grid);
figure, surf(XX2Grid,ZZ2Grid,tempdata2,'edgecolor','none');
axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

tempdata = squeeze(e13temp(:,1+28-Rad,:));
tempdata2 = gridfit(reshape(XXGrid,(M-2*Rad)*(L-2*Rad),1),reshape(ZZGrid,(M-2*Rad)*(L-2*Rad),1),reshape(tempdata,(M-2*Rad)*(L-2*Rad),1),XX2Grid,ZZ2Grid);
figure, surf(XX2Grid,ZZ2Grid,tempdata2,'edgecolor','none');
axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

% ============================================= 
[YYGrid,ZZGrid] = ndgrid(320+Rad*8:8:728-Rad*8,20+Rad*8:8:164-Rad*8); 
YY2Grid = interp(YYGrid(:,1)',4); ZZ2Grid = interp(ZZGrid(1,:),4);

tempdata = squeeze(e33temp(1+23-Rad,:,:));
tempdata2 = gridfit(reshape(YYGrid,(N-2*Rad)*(L-2*Rad),1),reshape(ZZGrid,(N-2*Rad)*(L-2*Rad),1),reshape(tempdata,(N-2*Rad)*(L-2*Rad),1),YY2Grid,ZZ2Grid);
figure, surf(YY2Grid,ZZ2Grid,tempdata2,'edgecolor','none');
axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

tempdata = squeeze(e12temp(1+23-Rad,:,:));
tempdata2 = gridfit(reshape(YYGrid,(N-2*Rad)*(L-2*Rad),1),reshape(ZZGrid,(N-2*Rad)*(L-2*Rad),1),reshape(tempdata,(N-2*Rad)*(L-2*Rad),1),YY2Grid,ZZ2Grid);
figure, surf(YY2Grid,ZZ2Grid,tempdata2,'edgecolor','none');
axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

tempdata = squeeze(e23temp(1+23-Rad,:,:));
tempdata2 = gridfit(reshape(YYGrid,(N-2*Rad)*(L-2*Rad),1),reshape(ZZGrid,(N-2*Rad)*(L-2*Rad),1),reshape(tempdata,(N-2*Rad)*(L-2*Rad),1),YY2Grid,ZZ2Grid);
figure, surf(YY2Grid,ZZ2Grid,tempdata2,'edgecolor','none');
axis equal; axis tight; view(2); colormap(coolwarm(128)); box on; set(gca,'fontsize',18); 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

    