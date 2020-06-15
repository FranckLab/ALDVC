% U0World = U0; % U0World(2:3:end) = -U0(2:3:end);  
% y0World = (size(Img{1},2)+1-xyz0.y); coordinatesFEMWorld = [coordinatesFEM(:,1),size(Img{1},2)+1-coordinatesFEM(:,2),coordinatesFEM(:,3)];
% close all; Plotuv(U0World,x0,y0World); % Plotdisp_show(U,coordinatesFEM,elementsFEM); % Plotuv(UExact,x0,y0);
% % figure(1); view([-140,60]);% colormap(coolwarm(32)); % caxis([-1.0,1.0]); % caxis([-0.45,0.45]);
% % figure(2); view([-140,60]);% colormap(coolwarm(32)); % caxis([-0.085, 0.005]); % caxis([-0.085,0.015 ]); 
% Plotdisp_show(U0World,coordinatesFEMWorld,elementsFEM); % Plot initial values
% figure(3); axis equal;% colormap(coolwarm(32)); colorbar; view(2); axis tight; set(gca,'fontsize',18); box on; % caxis([-1.0,1.0]); % caxis([-0.45,0.45]);
% figure(4); axis equal; %colormap(coolwarm(32)); colorbar; view(2); axis tight; set(gca,'fontsize',18); box on; % caxis([-0.085, 0.005]); % caxis([-0.085,0.015 ]); 

Plotdisp_show3(U0,coordinatesFEM,elementsFEM);


%% ------ Exact deformation ------
% % save(['imposed_disp_series','.mat'],'u','x0');
% % Generate model
% dvcExact = createpde(1);
% 
% % Apply mesh
% DTExact = delaunayTriangulation(x0{1}(:,1), x0{1}(:,2), x0{1}(:,3));
% geometryFromMesh(dvcExact,DTExact.Points',DTExact.ConnectivityList');
% 
% % Plot dvcExact 
% figure, pdeplot3D(dvcExact,'ColorMapData',u{4}(:,1),'FaceAlpha',0.5);
% title('Exact $x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
% set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex';
% set(gcf,'color','w'); colormap jet; colorbar; box on
% 
% 
% figure, pdeplot3D(dvcExact,'ColorMapData',u{4}(:,2),'FaceAlpha',0.5);
% title('Exact $y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
% set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex';
% set(gcf,'color','w'); colormap jet; colorbar; box on
% 
% 
% figure, pdeplot3D(dvcExact,'ColorMapData',u{4}(:,3),'FaceAlpha',0.5);
% title('Exact $z-$displacement $w$','FontWeight','Normal','Interpreter','latex');
% set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex';
% set(gcf,'color','w'); colormap jet; colorbar; box on
