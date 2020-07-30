function plotStreamline3(xGrid,yGrid,zGrid,u3x_meas_Grid,u3y_meas_Grid,u3z_meas_Grid,xGridsl,yGridsl,zGridsl)
% Plot stream lines of solved displacements
%
% References
% [1] Matlab Exchange File: https://www.mathworks.com/matlabcentral/fileexchange/24049-streamcolor
% -------------------------------------------
% Author: Jin Yang, aldicdvc@gmail.com
% Date: 2020.07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


u3_mag_meas_Grid = sqrt(u3x_meas_Grid.^2 + u3y_meas_Grid.^2 + u3z_meas_Grid.^2);

figure,streamcolor(  xGrid,  yGrid, zGrid, ...
    (u3x_meas_Grid), u3y_meas_Grid,  u3z_meas_Grid, ...
     (xGridsl), (yGridsl), (zGridsl), u3_mag_meas_Grid  );
 

u3_mag_max = max(u3_mag_meas_Grid(:));
cb=colorbar('Ticks',[0,0.25,0.5,0.75,1],'TickLabels', ...
    {num2str(0,'%4.0f'),num2str(0.25*u3_mag_max,'%4.0f'),num2str(0.5*u3_mag_max,'%4.0f'), ...
    num2str(0.75*u3_mag_max,'%4.0f'),num2str(u3_mag_max,'%4.0f')});


view([60,30]); set(gca,'fontsize',14); set(cb,'fontsize',16);
ylabel('x axis'); xlabel('y axis'); zlabel('z axis');  
grid minor; grid on; set(gcf,'color','w'); 
