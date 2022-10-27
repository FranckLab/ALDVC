%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slice plot 3D displacement using grid data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plotstrain3(FStrainTensor,x,y,z,varargin)
%FUNCTION Plotstrain3(F,x,y,z,varargin)
% To plot DVC 3D volumetric strain fields
% -------------------------------------------------
%
%   INPUT: FStrainTensor    Strain tensor was re-written as a vector: 
%                           FStrainTensor = [..., F11_nodei, F21_nodei, F31_nodei, ...
%                                                 F12_nodei, F22_nodei, F32_nodei, ...
%                                                 F13_nodei, F23_nodei, F33_nodei, ...]';
%          x,y,z            x-, y-, & z-coordinates of nodal points 
%          individualOrAll  Parameter to plot each individual displacement
%                           component or plot all components in the same figure
%          colorMap         User selected colormap matrix (i.e., a [Nx3] matrix) or
%                           a colormap name (i.e., 'jet')
%          caxisLimit       To set colormap limits, it can be a [1x2] array or a [3,2] matrix as 
%                           [caxis_limit_1_plot_dispx, caxis_limit_2_plot_dispx; 
%                            caxis_limit_1_plot_dispy, caxis_limit_2_plot_dispy;
%                            caxis_limit_1_plot_dispz, caxis_limit_2_plot_dispz].  
%
%   OUTPUT: Plots of strain fields.
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jin.yang@austin.utexas.edu -or- aldicdvc@gmail.com
% Last time updated: 10/2022.
% ==============================================

[individualOrAll,colorMap,caxisLimit] = parseargs(varargin);

%%
% ------ Load coordinates ------
M = size(x,1); N = size(y,2); L = size(z,3);
 
 
% Generate model
dvcVOI = createpde(1);

% Apply mesh
DT = delaunayTriangulation(x(:), y(:), z(:));
geometryFromMesh(dvcVOI,DT.Points',DT.ConnectivityList');
% ------ FEMesh has structure ------
% FEMesh with properties:
%
%              Nodes: [3xN double]
%           Elements: [4xM double]
%     GeometricOrder: 'linear'
% ----------------------------------

%%
if (strcmp(IndiOrAll,'Individual')==1) || (strcmp(IndiOrAll,'individual')==1)
    
    %% %%%%% Plot e11 %%%%%%%
    figure,  
    h = pdeplot3D(dvcVOI,'ColorMapData',FStrainTensor(1:9:end),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{11}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w');
    % xlabel('$x$ (pixels)','Interpreter','latex'); 
    % ylabel('$y$ (pixels)','Interpreter','latex'); 
    % zlabel('$z$ (pixels)','Interpreter','latex');
    
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(1,:)); end
    
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar('Location','eastoutside'); b.TickLabelInterpreter = 'latex';
    
    %%%%%% Re-locate axis and tick positions %%%%%%%
    a.Position = a.Position - [0 0 .03 .03];
    b.Position = b.Position + [.05 0 0 0];
     
    %%%%%% Make the view look tight %%%%%%%
    xlim([x(1,1,1),x(M,N,L)]); 
    ylim([y(1,1,1),y(M,N,L)]); 
    zlim([z(1,1,1),z(M,N,L)]);
     
    %% %%%%% Plot e22 %%%%%%%
    figure,  
    pdeplot3D(dvcVOI,'ColorMapData',FStrainTensor(5:9:end),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{22}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w');
     
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(1,:)); end
    
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar('Location','eastoutside'); b.TickLabelInterpreter = 'latex';
    
    %%%%%% Re-locate axis and tick positions %%%%%%%
    a.Position = a.Position - [0 0 .03 .03];
    b.Position = b.Position + [.05 0 0 0];
     
    %%%%%% Make the view look tight %%%%%%%
    xlim([x(1,1,1),x(M,N,L)]); 
    ylim([y(1,1,1),y(M,N,L)]); 
    zlim([z(1,1,1),z(M,N,L)]);
    
    
    %% %%%%% Plot e33 %%%%%%%
    figure,  
    pdeplot3D(dvcVOI,'ColorMapData',FStrainTensor(9:9:end),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{33}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w');
    
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(1,:)); end
    
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar('Location','eastoutside'); b.TickLabelInterpreter = 'latex';
    
    %%%%%% Re-locate axis and tick positions %%%%%%%
    a.Position = a.Position - [0 0 .03 .03];
    b.Position = b.Position + [.05 0 0 0];
     
    %%%%%% Make the view look tight %%%%%%%
    xlim([x(1,1,1),x(M,N,L)]); 
    ylim([y(1,1,1),y(M,N,L)]); 
    zlim([z(1,1,1),z(M,N,L)]);
    
    
    %% %%%%% Plot e12 %%%%%%%
    figure, 
    pdeplot3D(dvcVOI,'ColorMapData',0.5*(FStrainTensor(2:9:end)+FStrainTensor(4:9:end)),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{12}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w');
    
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(1,:)); end
    
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar('Location','eastoutside'); b.TickLabelInterpreter = 'latex';
    
    %%%%%% Re-locate axis and tick positions %%%%%%%
    a.Position = a.Position - [0 0 .03 .03];
    b.Position = b.Position + [.05 0 0 0];
     
    %%%%%% Make the view look tight %%%%%%%
    xlim([x(1,1,1),x(M,N,L)]); 
    ylim([y(1,1,1),y(M,N,L)]); 
    zlim([z(1,1,1),z(M,N,L)]);
    
    
    %% %%%%% Plot e13 %%%%%%%
    figure,  
    pdeplot3D(dvcVOI,'ColorMapData',0.5*(FStrainTensor(3:9:end)+FStrainTensor(7:9:end)),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{13}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w');
    
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(1,:)); end
    
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar('Location','eastoutside'); b.TickLabelInterpreter = 'latex';
    
    %%%%%% Re-locate axis and tick positions %%%%%%%
    a.Position = a.Position - [0 0 .03 .03];
    b.Position = b.Position + [.05 0 0 0];
     
    %%%%%% Make the view look tight %%%%%%%
    xlim([x(1,1,1),x(M,N,L)]); 
    ylim([y(1,1,1),y(M,N,L)]); 
    zlim([z(1,1,1),z(M,N,L)]);
    
    
    %% %%%%% Plot e23 %%%%%%%
    figure,  
    pdeplot3D(dvcVOI,'ColorMapData',0.5*(FStrainTensor(6:9:end)+FStrainTensor(8:9:end)),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{23}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w');
    
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(1,:)); end
    
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar('Location','eastoutside'); b.TickLabelInterpreter = 'latex';
    
    %%%%%% Re-locate axis and tick positions %%%%%%%
    a.Position = a.Position - [0 0 .03 .03];
    b.Position = b.Position + [.05 0 0 0];
     
    %%%%%% Make the view look tight %%%%%%%
    xlim([x(1,1,1),x(M,N,L)]); 
    ylim([y(1,1,1),y(M,N,L)]); 
    zlim([z(1,1,1),z(M,N,L)]);
    
    
else %%%%%% Plot all figures on the same page
    
    %% %%%%% Plot e11 %%%%%%%%
    figure, subplot(2,3,1),
    h=pdeplot3D(dvcVOI,'ColorMapData',FStrainTensor(1:9:end),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{11}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w');
    
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(1,:)); end

    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';

    xlim([x(1,1,1),x(M,N,L)]); 
    ylim([y(1,1,1),y(M,N,L)]); 
    zlim([z(1,1,1),z(M,N,L)]);


    %% %%%%% Plot e22 %%%%%%%%
    hold on; subplot(2,3,2),
    pdeplot3D(dvcVOI,'ColorMapData',FStrainTensor(5:9:end),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{22}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w');
    
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(1,:)); end

    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';

    xlim([x(1,1,1),x(M,N,L)]); 
    ylim([y(1,1,1),y(M,N,L)]); 
    zlim([z(1,1,1),z(M,N,L)]);
    
    %% %%%%% Plot e33 %%%%%%%
    hold on; subplot(2,3,3),
    pdeplot3D(dvcVOI,'ColorMapData',FStrainTensor(9:9:end),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{33}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w');
    
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(1,:)); end

    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';

    xlim([x(1,1,1),x(M,N,L)]); 
    ylim([y(1,1,1),y(M,N,L)]); 
    zlim([z(1,1,1),z(M,N,L)]);
 
    %% %%%%% Plot e12 %%%%%%%
    hold on; subplot(2,3,4),
    pdeplot3D(dvcVOI,'ColorMapData',0.5*(FStrainTensor(2:9:end)+FStrainTensor(4:9:end)),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{12}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w');
    
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(1,:)); end

    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';

    xlim([x(1,1,1),x(M,N,L)]); 
    ylim([y(1,1,1),y(M,N,L)]); 
    zlim([z(1,1,1),z(M,N,L)]);
    
    %% %%%%% Plot e13 %%%%%%%
    hold on; subplot(2,3,5),
    pdeplot3D(dvcVOI,'ColorMapData',0.5*(FStrainTensor(3:9:end)+FStrainTensor(7:9:end)),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{13}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w');
    
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(1,:)); end

    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';

    xlim([x(1,1,1),x(M,N,L)]); 
    ylim([y(1,1,1),y(M,N,L)]); 
    zlim([z(1,1,1),z(M,N,L)]);
    
    %% %%%%% Plot e23 %%%%%%%
    hold on; subplot(2,3,6),
    pdeplot3D(dvcVOI,'ColorMapData',0.5*(FStrainTensor(6:9:end)+FStrainTensor(8:9:end)),'FaceAlpha',1.0,'FaceColor','interp');
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('Strain $e_{23}$','fontweight','normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w');
    
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(1,:)); end

    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';

    xlim([x(1,1,1),x(M,N,L)]); 
    ylim([y(1,1,1),y(M,N,L)]); 
    zlim([z(1,1,1),z(M,N,L)]);
    
end


