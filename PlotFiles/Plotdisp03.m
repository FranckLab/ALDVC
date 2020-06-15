% ==============================================
% function Plotdisp_show3 for 3D FEM data
% ==============================================
function Plotdisp03(U,coordinatesFEM,elementsFEM,IndiOrAll)

% Generate model
dvcVOI = createpde(1);

% Apply mesh
DT = delaunayTriangulation(coordinatesFEM(:,1), coordinatesFEM(:,2), coordinatesFEM(:,3));
geometryFromMesh(dvcVOI,DT.Points',DT.ConnectivityList');
% ------ FEMesh has structure ------
% FEMesh with properties:
%
%              Nodes: [3x10003 double]
%           Elements: [10x5774 double]
%     MaxElementSize: 9.7980
%     MinElementSize: 4.8990
%      MeshGradation: 1.5000
%     GeometricOrder: 'quadratic'
% ----------------------------------


if strcmp(IndiOrAll,'Individual') == 1
    % plot dvcZOI domain
    % figure, pdegplot(dvcZOI,'FaceLabels','on','FaceAlpha',0.5); set(gca,'fontsize',18);title('ZOI body')
    figure,
    % ------------------ u -----------------------
    %subplot(2,3,1),
    pdeplot3D(dvcVOI,'ColorMapData',U(1:3:end),'FaceAlpha',1.0);
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
    % xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    set(gcf,'color','w'); colormap jet; %coolwarm(32);
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar('Location','eastoutside'); b.TickLabelInterpreter = 'latex';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a.Position = a.Position - [0 0 .03 .03];
    b.Position = b.Position + [.05 0 0 0];
    %zticks('')
    %set(gca,'fontsize',20)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xlim([min(coordinatesFEM(:,1)),max(coordinatesFEM(:,1))]);
    ylim([min(coordinatesFEM(:,2)),max(coordinatesFEM(:,2))]);
    zlim([min(coordinatesFEM(:,3)),max(coordinatesFEM(:,3))]);
    
    %
    % pause;
    
    % ------------------ v -----------------------
    figure, %subplot(2,3,2),
    pdeplot3D(dvcVOI,'ColorMapData',U(2:3:end),'FaceAlpha',1.0);
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
    % xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    set(gcf,'color','w'); colormap jet %coolwarm(32);
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    a.Position = a.Position - [0 0 .03 .03];
    b.Position = b.Position + [.05 0 0 0];
    
    xlim([min(coordinatesFEM(:,1)),max(coordinatesFEM(:,1))]);
    ylim([min(coordinatesFEM(:,2)),max(coordinatesFEM(:,2))]);
    zlim([min(coordinatesFEM(:,3)),max(coordinatesFEM(:,3))]);
    
    
    
    % ------------------ w -----------------------
    figure, %subplot(2,3,3),
    pdeplot3D(dvcVOI,'ColorMapData',U(3:3:end),'FaceAlpha',1.0);
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('$z-$displacement $w$','FontWeight','Normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
    % xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    set(gcf,'color','w'); colormap jet %coolwarm(32);
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    a.Position = a.Position - [0 0 .03 .03];
    b.Position = b.Position + [.05 0 0 0];
    
    xlim([min(coordinatesFEM(:,1)),max(coordinatesFEM(:,1))]);
    ylim([min(coordinatesFEM(:,2)),max(coordinatesFEM(:,2))]);
    zlim([min(coordinatesFEM(:,3)),max(coordinatesFEM(:,3))]);
    
else
    
    % plot dvcZOI domain
    % figure, pdegplot(dvcZOI,'FaceLabels','on','FaceAlpha',0.5); set(gca,'fontsize',18);title('ZOI body')
    figure,
    % ------------------ u -----------------------
    subplot(2,3,1),
    pdeplot3D(dvcVOI,'ColorMapData',U(1:3:end),'FaceAlpha',1.0);
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
    % xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    set(gcf,'color','w'); colormap jet; %coolwarm(32);
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([min(coordinatesFEM(:,1)),max(coordinatesFEM(:,1))]);
    ylim([min(coordinatesFEM(:,2)),max(coordinatesFEM(:,2))]);
    zlim([min(coordinatesFEM(:,3)),max(coordinatesFEM(:,3))]);
    
    %
    % pause;
    
    % ------------------ v -----------------------
    subplot(2,3,2),
    pdeplot3D(dvcVOI,'ColorMapData',U(2:3:end),'FaceAlpha',1.0);
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
    % xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    set(gcf,'color','w'); colormap jet %coolwarm(32);
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([min(coordinatesFEM(:,1)),max(coordinatesFEM(:,1))]);
    ylim([min(coordinatesFEM(:,2)),max(coordinatesFEM(:,2))]);
    zlim([min(coordinatesFEM(:,3)),max(coordinatesFEM(:,3))]);
    
    
    
    % ------------------ w -----------------------
    subplot(2,3,3),
    pdeplot3D(dvcVOI,'ColorMapData',U(3:3:end),'FaceAlpha',1.0);
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('$z-$displacement $w$','FontWeight','Normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
    % xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    set(gcf,'color','w'); colormap jet %coolwarm(32);
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([min(coordinatesFEM(:,1)),max(coordinatesFEM(:,1))]);
    ylim([min(coordinatesFEM(:,2)),max(coordinatesFEM(:,2))]);
    zlim([min(coordinatesFEM(:,3)),max(coordinatesFEM(:,3))]);
    
end




%% ---------------------------------
% figure, pdeplot3D(dvcZOI,'ColorMapData',U(2:3:end),'FaceAlpha',0.5);
%
% title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
% set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
%
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex';
%
% set(gcf,'color','w'); colormap jet; colorbar; box on
%
%
% % ---------------------------------
% figure, pdeplot3D(dvcZOI,'ColorMapData',U(3:3:end),'FaceAlpha',0.5);
%
% title('$z-$displacement $w$','FontWeight','Normal','Interpreter','latex');
% set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
%
% a = gca; a.TickLabelInterpreter = 'latex';
% b = colorbar; b.TickLabelInterpreter = 'latex';
%
% set(gcf,'color','w'); colormap jet; colorbar; box on
%



