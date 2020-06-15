% ==============================================
% function Plotdisp_show3 for 3D FEM data
% ==============================================
function Plotstrain_show3(F,coordinatesFEM,elementsFEM)

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

% plot dvcZOI domain
% ---------------------------------
figure,
subplot(2,3,1), pdeplot3D(dvcVOI,'ColorMapData',F(1:9:end),'FaceAlpha',0.5);

title('Strain $e_{11}$','fontweight','normal','Interpreter','latex');  
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
set(gcf,'color','w'); colormap jet; colorbar; box on
 
% ---------------------------------
subplot(2,3,2), pdeplot3D(dvcVOI,'ColorMapData',F(5:9:end),'FaceAlpha',0.5);

title('Strain $e_{22}$','fontweight','normal','Interpreter','latex');  
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
set(gcf,'color','w'); colormap jet; colorbar; box on

% ---------------------------------
subplot(2,3,3), pdeplot3D(dvcVOI,'ColorMapData',F(9:9:end),'FaceAlpha',0.5);

title('Strain $e_{33}$','fontweight','normal','Interpreter','latex');  
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
set(gcf,'color','w'); colormap jet; colorbar; box on

% ---------------------------------
subplot(2,3,4), pdeplot3D(dvcVOI,'ColorMapData',0.5*(F(2:9:end)+F(4:9:end)),'FaceAlpha',0.5);

title('Strain $e_{12}$','fontweight','normal','Interpreter','latex');  
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
set(gcf,'color','w'); colormap jet; colorbar; box on

% ---------------------------------
subplot(2,3,5), pdeplot3D(dvcVOI,'ColorMapData',0.5*(F(3:9:end)+F(7:9:end)),'FaceAlpha',0.5);

title('Strain $e_{13}$','fontweight','normal','Interpreter','latex');  
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
set(gcf,'color','w'); colormap jet; colorbar; box on

% ---------------------------------
subplot(2,3,6), pdeplot3D(dvcVOI,'ColorMapData',0.5*(F(6:9:end)+F(8:9:end)),'FaceAlpha',0.5);

title('Strain $e_{23}$','fontweight','normal','Interpreter','latex');  
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
set(gcf,'color','w'); colormap jet; colorbar; box on




