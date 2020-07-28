function plotCone3(xyz,uvw)
% Cone plot solved displacements
%
% -------------------------------------------
% Author: Jin Yang, aldicdvc@gmail.com
% Date: 2020.07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=xyz(:,1); y=xyz(:,2); z=xyz(:,3);
u=uvw(1:3:end); v=uvw(2:3:end); w=uvw(3:3:end);
 

% === Cone plot the solution ===
figure,  u_mag = sqrt(u.^2 + v.^2 + w.^2); Umag_max = max(u_mag);
hc = coneplot(x,y,z,u,v,w,0.04,'nointerp'); 
  caxis([0,Umag_max]); fvc = repmat(u_mag(:)',[42 1]);
set(hc, 'facecolor', 'flat', 'facevertexcdata', fvc(:));
hc.EdgeColor = 'none'; hc.AmbientStrength = 0.6; hc.DiffuseStrength = 0.75;  
hc.SpecularStrength = 0.4; hc.SpecularExponent = 3; 
h_color = colorbar; set(h_color, 'ylim', [0, Umag_max]); set(h_color,'fontsize',14);

lighting phong;
title('Displacement (unit: vx)','fontweight','normal')

view([60,30]); set(gca,'fontsize',14);  axis on; 
ylabel('y'); xlabel('x'); zlabel('z');  colormap(jet);
 grid minor; grid on; set(gcf,'color','w');