function Plotdisp3(U,coordinatesXYZ,varargin)
%FUNCTION Plotdisp3(U,coordinatesXYZ,varargin)
% To plot DVC 3D volumetric displacement fields 
% -------------------------------------------------
%
%   INPUT: U                Displacement vector: U = [Ux_node1, Uy_node1, Uz_node1, ...
%                                                     Ux_node2, Uy_node2, Uz_node2, ...   
%                                                     ..., 
%                                                     Ux_nodeN, Uy_nodeN, Uz_nodeN]';
%          coordinatesXYZ   Coordinates of nodal points: 
%                           [x-coord_node1, y-coord_node1, z-coord_node1;  
%                            x-coord_node2, y-coord_node2, z-coord_node2;  
%                            ...;
%                            x-coord_nodeN, y-coord_nodeN, z-coord_nodeN];
%          individualOrAll  Parameter to plot each individual displacement
%                           component or plot all components in the same figure
%          colorMap         User selected colormap matrix (i.e., a [Nx3] matrix) or
%                           a colormap name (i.e., 'jet')
%          caxisLimit       To set colormap limits, it can be a [1x2] array or a [3,2] matrix as 
%                           [caxis_limit_1_plot_dispx, caxis_limit_2_plot_dispx; 
%                            caxis_limit_1_plot_dispy, caxis_limit_2_plot_dispy;
%                            caxis_limit_1_plot_dispz, caxis_limit_2_plot_dispz].  
%
%   OUTPUT: Plots of x-, y-, and z-displacement fields.
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jin.yang@austin.utexas.edu -or- aldicdvc@gmail.com
% Last time updated: 10/2022.
% ==============================================

[individualOrAll,colorMap,caxisLimit] = parseargs(varargin);


%%
% Generate model
dvcVOI = createpde(1);

% Apply mesh
DT = delaunayTriangulation(coordinatesXYZ(:,1), coordinatesXYZ(:,2), coordinatesXYZ(:,3));
geometryFromMesh(dvcVOI,DT.Points',DT.ConnectivityList');
% ------ FEMesh has structure ------
% FEMesh with properties:
%
%              Nodes: [3xN double]
%           Elements: [4xM double]
%     GeometricOrder: 'linear'
% ----------------------------------

%%
if (strcmp(individualOrAll,'Individual')==1) || (strcmp(individualOrAll,'individual')==1)
    
    %%%%%% Plot dvcZOI body %%%%%%
    % figure, pdegplot(dvcVOI,'FaceLabels','on','FaceAlpha',0.5); 
    % set(gca,'fontsize',18);title('ZOI body')

    %% %%%% Plot x-disp %%%%%%%
    figure, pdeplot3D(dvcVOI,'ColorMapData',U(1:3:end),'FaceAlpha',1.0);
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;  axis equal; set(gcf,'color','w');  
    % xlabel('$x$ (pixels)','Interpreter','latex'); 
    % ylabel('$y$ (pixels)','Interpreter','latex');
    
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(1,:)); end
    
    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar('Location','eastoutside'); b.TickLabelInterpreter = 'latex';
    
    %%%%%% Re-locate axis and tick positions %%%%%%%
    a.Position = a.Position - [0 0 .03 .03];
    b.Position = b.Position + [.05 0 0 0];
     
    %%%%%% Make the view look tight %%%%%%%
    xlim([min(coordinatesXYZ(:,1)),max(coordinatesXYZ(:,1))]);
    ylim([min(coordinatesXYZ(:,2)),max(coordinatesXYZ(:,2))]);
    zlim([min(coordinatesXYZ(:,3)),max(coordinatesXYZ(:,3))]);
     
    
    %% %%%% Plot y-disp %%%%%%%
    figure, pdeplot3D(dvcVOI,'ColorMapData',U(2:3:end),'FaceAlpha',1.0);
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;  axis equal; set(gcf,'color','w'); 
    % xlabel('$x$ (pixels)','Interpreter','latex'); 
    % ylabel('$y$ (pixels)','Interpreter','latex');
     
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(2,:)); end

    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    %%%%%% Re-locate axis and tick positions %%%%%%%
    a.Position = a.Position - [0 0 .03 .03];
    b.Position = b.Position + [.05 0 0 0];
    
    %%%%%% Make the view look tight %%%%%%%
    xlim([min(coordinatesXYZ(:,1)),max(coordinatesXYZ(:,1))]);
    ylim([min(coordinatesXYZ(:,2)),max(coordinatesXYZ(:,2))]);
    zlim([min(coordinatesXYZ(:,3)),max(coordinatesXYZ(:,3))]);
    
    
    
    %% %%%% Plot z-disp %%%%%%%
    figure, pdeplot3D(dvcVOI,'ColorMapData',U(3:3:end),'FaceAlpha',1.0);
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('$z-$displacement $w$','FontWeight','Normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on;  axis equal; set(gcf,'color','w'); 
    % xlabel('$x$ (pixels)','Interpreter','latex'); 
    % ylabel('$y$ (pixels)','Interpreter','latex');
    
    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(3,:)); end

    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    %%%%%% Re-locate axis and tick positions %%%%%%%
    a.Position = a.Position - [0 0 .03 .03];
    b.Position = b.Position + [.05 0 0 0];
    
    %%%%%% Make the view look tight %%%%%%%
    xlim([min(coordinatesXYZ(:,1)),max(coordinatesXYZ(:,1))]);
    ylim([min(coordinatesXYZ(:,2)),max(coordinatesXYZ(:,2))]);
    zlim([min(coordinatesXYZ(:,3)),max(coordinatesXYZ(:,3))]);
    

else %%%%%% Plot all figures on the same page
    
    %% %%%% Plot x-disp %%%%%%%
    figure, subplot(2,3,1),
    pdeplot3D(dvcVOI,'ColorMapData',U(1:3:end),'FaceAlpha',1.0);
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w');

    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(1,:)); end

    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([min(coordinatesXYZ(:,1)),max(coordinatesXYZ(:,1))]);
    ylim([min(coordinatesXYZ(:,2)),max(coordinatesXYZ(:,2))]);
    zlim([min(coordinatesXYZ(:,3)),max(coordinatesXYZ(:,3))]);
     
    %% %%%% Plot y-disp %%%%%%%
    hold on, subplot(2,3,2),
    pdeplot3D(dvcVOI,'ColorMapData',U(2:3:end),'FaceAlpha',1.0);
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w'); 

    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(2,:)); end

    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([min(coordinatesXYZ(:,1)),max(coordinatesXYZ(:,1))]);
    ylim([min(coordinatesXYZ(:,2)),max(coordinatesXYZ(:,2))]);
    zlim([min(coordinatesXYZ(:,3)),max(coordinatesXYZ(:,3))]);
    
    %% %%%% Plot z-disp %%%%%%%
    hold on, subplot(2,3,3),
    pdeplot3D(dvcVOI,'ColorMapData',U(2:3:end),'FaceAlpha',1.0);
    hold on, pdegplot(dvcVOI,'FaceLabel','off','FaceAlpha',0.3);
    title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
    set(gca,'fontsize',18); axis on; axis equal; set(gcf,'color','w'); 

    try colormap(colorMap); catch; end
    if ~isempty(caxisLimit), caxis(caxisLimit(3,:)); end

    a = gca; a.TickLabelInterpreter = 'latex';
    b = colorbar; b.TickLabelInterpreter = 'latex';
    
    xlim([min(coordinatesXYZ(:,1)),max(coordinatesXYZ(:,1))]);
    ylim([min(coordinatesXYZ(:,2)),max(coordinatesXYZ(:,2))]);
    zlim([min(coordinatesXYZ(:,3)),max(coordinatesXYZ(:,3))]);
    
end

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [individualOrAll,colorMap,caxisLimit] = parseargs(vargin)

individualOrAll = []; colorMap = []; caxisLimit = [];
 
try 
    individualOrAll = vargin{1};     
catch
    individualOrAll = 'all';
end

try 
    colorMap = vargin{2};
    if isempty(colorMap), colorMap = 'turbo'; end
catch
    colorMap = 'turbo';
end

try 
    caxisLimit = vargin{3};
    if size(caxisLimit,1)==1
        caxisLimit = repmat(caxisLimit,3,1);
    end
catch     
end

end


