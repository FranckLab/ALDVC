%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slice plot 3D displacement using grid data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plotstrain03(FLocal,x,y,z,sizeOfImg,IndiOrAll)

% ------ Load coordinates ------
M = size(x,1); N = size(y,2); L = size(z,3);
unitx = (-x(1,1,1)+x(M,N,L))/(M-1); unity = (-y(1,1,1)+y(M,N,L))/(N-1); unitz = (-z(1,1,1)+z(M,N,L))/(L-1);

% ------ Load strain info ------


% Generate model
dvcVOI = createpde(1);
% Apply mesh
DT = delaunayTriangulation(x(:), y(:), z(:));
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
    
    % ------------------ e11 -----------------------
    figure, %subplot(1,1,1),
    h=pdeplot3D(dvcVOI,'ColorMapData',FLocal(1:9:end),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{11}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
    % xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex'); zlabel('$z$ (pixels)','Interpreter','latex');
    set(gcf,'color','w'); %colormap coolwarm(32);  %caxis([-0.01,0.01]);
    caxis auto;
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([x(1,1,1),x(M,N,L)]); ylim([y(1,1,1),y(M,N,L)]); zlim([z(1,1,1),z(M,N,L)]);
    
    % pause;
    
    % ------------------ e22 -----------------------
    figure, %subplot(2,3,2),
    pdeplot3D(dvcVOI,'ColorMapData',FLocal(5:9:end),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{22}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
    %xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    
    set(gcf,'color','w'); %colormap coolwarm(32);  %caxis([-0.01,0.01]);
    caxis auto;
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([x(1,1,1),x(M,N,L)]); ylim([y(1,1,1),y(M,N,L)]); zlim([z(1,1,1),z(M,N,L)]);
    
    
    % ------------------ e33 -----------------------
    figure, %subplot(2,3,3),
    pdeplot3D(dvcVOI,'ColorMapData',FLocal(9:9:end),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{33}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;   axis equal;  % view([90 -90])
    %xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    
    set(gcf,'color','w'); %colormap coolwarm(32);  %caxis([-0.01,0.01]);
    caxis auto;
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([x(1,1,1),x(M,N,L)]); ylim([y(1,1,1),y(M,N,L)]); zlim([z(1,1,1),z(M,N,L)]);
    
    
    % ------------------ e12 -----------------------
    figure, %subplot(2,3,4),
    pdeplot3D(dvcVOI,'ColorMapData',0.5*(FLocal(2:9:end)+FLocal(4:9:end)),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{12}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;   axis equal;  % view([90 -90])
    %xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    
    set(gcf,'color','w'); %colormap coolwarm(32);  %caxis([-0.01,0.01]);
    caxis auto;
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([x(1,1,1),x(M,N,L)]); ylim([y(1,1,1),y(M,N,L)]); zlim([z(1,1,1),z(M,N,L)]);
    
    
    % ------------------ e13 -----------------------
    figure, %subplot(2,3,5),
    pdeplot3D(dvcVOI,'ColorMapData',0.5*(FLocal(3:9:end)+FLocal(7:9:end)),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{13}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;   axis equal;  % view([90 -90])
    %xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    
    set(gcf,'color','w');
    %colormap coolwarm(32);  %caxis([-0.01,0.01]);
    caxis auto;
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([x(1,1,1),x(M,N,L)]); ylim([y(1,1,1),y(M,N,L)]); zlim([z(1,1,1),z(M,N,L)]);
    
    
    % ------------------ e23 -----------------------
    figure, %subplot(2,3,6),
    pdeplot3D(dvcVOI,'ColorMapData',0.5*(FLocal(6:9:end)+FLocal(8:9:end)),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{23}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;   axis equal;  % view([90 -90])
    %xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    
    set(gcf,'color','w');
    %colormap coolwarm(32);  %caxis([-0.01,0.01]);
    caxis auto;
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([x(1,1,1),x(M,N,L)]); ylim([y(1,1,1),y(M,N,L)]); zlim([z(1,1,1),z(M,N,L)]);
    
    
else
    
    
    figure,
    % ------------------ e11 -----------------------
    subplot(2,3,1),
    h=pdeplot3D(dvcVOI,'ColorMapData',FLocal(1:9:end),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{11}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
    % xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex'); zlabel('$z$ (pixels)','Interpreter','latex');
    set(gcf,'color','w'); %colormap coolwarm(32);  %caxis([-0.01,0.01]);
    caxis auto; % caxis([0.098,0.102]) % 
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([x(1,1,1),x(M,N,L)]); ylim([y(1,1,1),y(M,N,L)]); zlim([z(1,1,1),z(M,N,L)]);
    
    % pause;
    
    
    % ------------------ e22 -----------------------
    subplot(2,3,2),
    pdeplot3D(dvcVOI,'ColorMapData',FLocal(5:9:end),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{22}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
    %xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    
    set(gcf,'color','w'); %colormap coolwarm(32);  %caxis([-0.01,0.01]);
    caxis auto;
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([x(1,1,1),x(M,N,L)]); ylim([y(1,1,1),y(M,N,L)]); zlim([z(1,1,1),z(M,N,L)]);
    
    
    % ------------------ e33 -----------------------
    subplot(2,3,3),
    pdeplot3D(dvcVOI,'ColorMapData',FLocal(9:9:end),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{33}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;   axis equal;  % view([90 -90])
    %xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    
    set(gcf,'color','w'); %colormap coolwarm(32);  %caxis([-0.01,0.01]);
    caxis auto;
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([x(1,1,1),x(M,N,L)]); ylim([y(1,1,1),y(M,N,L)]); zlim([z(1,1,1),z(M,N,L)]);
    
    
    % ------------------ e12 -----------------------
    subplot(2,3,4),
    pdeplot3D(dvcVOI,'ColorMapData',0.5*(FLocal(2:9:end)+FLocal(4:9:end)),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{12}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;   axis equal;  % view([90 -90])
    %xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    
    set(gcf,'color','w');%colormap coolwarm(32);  %caxis([-0.01,0.01]);
    caxis auto;
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([x(1,1,1),x(M,N,L)]); ylim([y(1,1,1),y(M,N,L)]); zlim([z(1,1,1),z(M,N,L)]);
    
    
    % ------------------ e13 -----------------------
    subplot(2,3,5),
    pdeplot3D(dvcVOI,'ColorMapData',0.5*(FLocal(3:9:end)+FLocal(7:9:end)),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{13}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;   axis equal;  % view([90 -90])
    %xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    
    set(gcf,'color','w');
    %colormap coolwarm(32);  %caxis([-0.01,0.01]);
    caxis auto;
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([x(1,1,1),x(M,N,L)]); ylim([y(1,1,1),y(M,N,L)]); zlim([z(1,1,1),z(M,N,L)]);
    
    
    % ------------------ e23 -----------------------
    subplot(2,3,6),
    pdeplot3D(dvcVOI,'ColorMapData',0.5*(FLocal(6:9:end)+FLocal(8:9:end)),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{23}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;   axis equal;  % view([90 -90])
    %xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    
    set(gcf,'color','w');
    colormap coolwarm(32);  %colormap coolwarm(32);  %caxis([-0.01,0.01]);
    caxis auto;
    try
        load('colormap_RdYlBu.mat'); colormap(cMap)
    catch
    end
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([x(1,1,1),x(M,N,L)]); ylim([y(1,1,1),y(M,N,L)]); zlim([z(1,1,1),z(M,N,L)]);
    
    
end


