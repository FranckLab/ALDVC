%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slice plot 3D displacement using grid data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plotuvw(ULocal,uvw,xyz,Phi)

% ------ Load coordinates ------
x = xyz.x; y = xyz.y; z = xyz.z;
M = size(x,1); N = size(y,2); L = size(z,3);
unitx = (-x(1,1,1)+x(M,N,L))/(M-1); unity = (-y(1,1,1)+y(M,N,L))/(N-1); unitz = (-z(1,1,1)+z(M,N,L))/(L-1);

% ------ Load displacement info ------
if isempty(ULocal) ~= 1
    u = reshape(ULocal(1:3:end),M,N,L);
    v = reshape(ULocal(2:3:end),M,N,L);
    w = reshape(ULocal(3:3:end),M,N,L);
else
   u = uvw.u; v = uvw.v; w = uvw.w;
end
 
% ------ Plot figures ------
[imgx,imgy,imgz] = meshgrid(x(1,1,1):unitx:x(M,N,L),y(1,1,1):unity:y(M,N,L),z(1,1,1):unitz:z(M,N,L));
figure; axes3D = axes; daspect([1 1 1]);
axis([x(1,1,1) x(M,N,L) y(1,1,1) y(M,N,L) z(1,1,1) z(M,N,L)]); hold on
% add slice planes
hSlice = slice( imgx,imgy,imgz,u,[ x(1,1,1) ],[ y(M,N,L) ], [ z(1,1,1), 0.5*(z(1,1,1)+z(M,N,L)) ] );
set(hSlice,'EdgeColor','none','FaceColor','interp');
title('x Displacement','fontweight','normal'); set(gca,'fontsize',18); view([50,25]); colorbar;

figure; axes3D = axes; daspect([1 1 1]);
axis([x(1,1,1) x(M,N,L) y(1,1,1) y(M,N,L) z(1,1,1) z(M,N,L)]); hold on
% add slice planes
hSlice = slice( imgx,imgy,imgz,v,[ x(1,1,1) ],[ y(M,N,L) ], [ z(1,1,1) 0.5*(z(1,1,1)+z(M,N,L)) ] );
set(hSlice,'EdgeColor','none','FaceColor','interp');
title('y Displacement','fontweight','normal'); set(gca,'fontsize',18); view([50,25]); colorbar;

figure; axes3D = axes; daspect([1 1 1]);
axis([x(1,1,1) x(M,N,L) y(1,1,1) y(M,N,L) z(1,1,1) z(M,N,L)]); hold on
% add slice planes
hSlice = slice( imgx,imgy,imgz,w,[ x(1,1,1) ],[ y(M,N,L) ], [ z(1,1,1) 0.5*(z(1,1,1)+z(M,N,L)) ] );
set(hSlice,'EdgeColor','none','FaceColor','interp');
title('z Displacement','fontweight','normal'); set(gca,'fontsize',18); view([50,25]); colorbar;


% ------ Plot cross correlation function ------ 
if isempty(Phi) ~= 1
    figure; axes3D = axes; daspect([1 1 1]);
    axis([x(1,1,1) x(M,N,L) y(1,1,1) y(M,N,L) z(1,1,1) z(M,N,L)]); hold on
    % add slice planes
    hSlice = slice( imgx,imgy,imgz,Phi,[ x(1,1,1) ],[ y(M,N,L) ], [ z(1,1,1), 0.5*(z(1,1,1)+z(M,N,L)) ] );
    set(hSlice,'EdgeColor','none','FaceColor','interp');
    title('Normalized cross correlation','fontweight','normal'); set(gca,'fontsize',18); view([50,25]); colorbar;
end


