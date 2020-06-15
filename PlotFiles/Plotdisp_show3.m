% ==============================================
% function Plotdisp_show3 for 3D FEM data
% ==============================================
function Plotdisp_show3(U,coordinatesFEM,elementsFEM)

% Generate model
dvcZOI = createpde(1);

% Apply mesh
DT = delaunayTriangulation(coordinatesFEM(:,1), coordinatesFEM(:,2), coordinatesFEM(:,3)); 
geometryFromMesh(dvcZOI,DT.Points',DT.ConnectivityList');
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
% figure, pdegplot(dvcZOI,'FaceLabels','on','FaceAlpha',0.5); set(gca,'fontsize',18);title('ZOI body')
figure, 
subplot(2,3,1)
pdeplot3D(dvcZOI,'ColorMapData',U(1:3:end),'FaceAlpha',0.5);

title('$x$-displacement $u$','FontWeight','Normal','Interpreter','latex');
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

set(gcf,'color','w'); colormap jet; colorbar; box on
 

% ---------------------------------
subplot(2,3,2), pdeplot3D(dvcZOI,'ColorMapData',U(2:3:end),'FaceAlpha',0.5);

title('$y$-displacement $v$','FontWeight','Normal','Interpreter','latex');
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

set(gcf,'color','w'); colormap jet; colorbar; box on

 
% ---------------------------------
subplot(2,3,3), pdeplot3D(dvcZOI,'ColorMapData',U(3:3:end),'FaceAlpha',0.5);

title('$z$-displacement $w$','FontWeight','Normal','Interpreter','latex');
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
 
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

set(gcf,'color','w'); colormap jet; colorbar; box on




