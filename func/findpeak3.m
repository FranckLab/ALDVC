function [xpeak,ypeak,zpeak,max_f] = findpeak3(f,subpixel)

%FINDPEAK Find extremum of matrix.
%   [XPEAK,YPEAK,ZPEAK,MAX_F] = FINDPEAK(F,SUBPIXEL) finds the extremum of F,
%   MAX_F, and its location (XPEAK, YPEAK, ZPEAK). F is a three dimensional matrix. 
%   MAX_F is the maximum value of F, or an estimate of the extremum if a subpixel
%   extremum is requested.
%
%   SUBPIXEL is a boolean that controls if FINDPEAK attempts to estimate the
%   extremum location to subpixel precision. If SUBPIXEL is false, FINDPEAK
%   returns the coordinates of the maximum absolute value of F and MAX_F is
%   max(abs(F(:))). If SUBPIXEL is true, FINDPEAK fits a 2nd order
%   polynomial to the maximum absolute value of F. 
%   In this case, MAX_F is the absolute value of the polynomial evaluated
%   at its extremum.
%
%   Note: Even if SUBPIXEL is true, there are some cases that result
%   in FINDPEAK returning the coordinates of the maximum absolute value
%   of F:
%   * When the maximum absolute value of F is on the edge of matrix F.
%   * When the coordinates of the estimated polynomial extremum would fall
%     outside the coordinates of the points used to constrain the estimate.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision $  $Date: 2004/10/20 17:54:47 $
%   $Revision $  $Date: 2017/10/22 21:15:00 $


% get peak pixel
[max_f, imax] = max(( f(:) ));
[ypeak, xpeak, zpeak] = ind2sub(size(f),imax);

% disp(['ypeak=',num2str(ypeak),' xpeak=',num2str(xpeak),' zpeak=',num2str(zpeak)]);

    
if ~subpixel || ...
    ypeak==1 || ypeak==size(f,1) || xpeak==1 || xpeak==size(f,2) || zpeak==1 || zpeak==size(f,3) % on edge
    return % return absolute peak
    
else
    
    % fit a 2nd order polynomial to 27 points  
    % using 27 pixels centered on irow,jcol    
    u = f(ypeak-1:ypeak+1, xpeak-1:xpeak+1, zpeak-1:zpeak+1);
    u = u(:);
    x = [-1 -1 -1  0  0  0  1  1  1 -1 -1 -1  0  0  0  1  1  1  -1 -1 -1  0  0  0  1  1  1  ]';
    y = [-1  0  1 -1  0  1 -1  0  1 -1  0  1 -1  0  1 -1  0  1  -1  0  1 -1  0  1 -1  0  1  ]';
    z = [-1 -1 -1 -1 -1 -1 -1 -1 -1  0  0  0  0  0  0  0  0  0   1  1  1  1  1  1  1  1  1  ]';
    
    % u(x,y,z) = A(1) + A(2)*x + A(3)*y + A(4)*z + A(5)*x*y + A(6)*x*z+ ...
    % A(7)*y*z + A(8)*x^2 + A(9)*y^2 + A(10)*z^2
    X = [ones(27,1), x, y, z, x.*y, x.*z, y.*z, x.^2, y.^2, z.^2];
    
    % u = X*A
    % opts.LT = false; opts.UT = false; opts.UHESS = false; opts.SYM = false; opts.RECT = true; opts.TRANSA = false;
    % A = linsolve(X,u,opts);
    % A = X\u;
    A = (X'*X)\(X'*u);
    
    % get maximum, where du/dx = du/dy = du/dz = 0
    x_offset = -(A(2)*A(7)^2 - A(3)*A(6)*A(7) - A(4)*A(5)*A(7) + 2*A(3)*A(5)*A(10) + 2*A(4)*A(6)*A(9) - 4*A(2)*A(9)*A(10))/(2*(A(10)*A(5)^2 - A(5)*A(6)*A(7) + A(9)*A(6)^2 + A(8)*A(7)^2 - 4*A(8)*A(9)*A(10)));
    y_offset = -(A(3)*A(6)^2 - A(2)*A(6)*A(7) - A(4)*A(5)*A(6) + 2*A(2)*A(5)*A(10) + 2*A(4)*A(7)*A(8) - 4*A(3)*A(8)*A(10))/(2*(A(10)*A(5)^2 - A(5)*A(6)*A(7) + A(9)*A(6)^2 + A(8)*A(7)^2 - 4*A(8)*A(9)*A(10)));
    z_offset = -(A(4)*A(5)^2 - A(2)*A(5)*A(7) - A(3)*A(5)*A(6) + 2*A(2)*A(6)*A(9) + 2*A(3)*A(7)*A(8) - 4*A(4)*A(8)*A(9))/(2*(A(10)*A(5)^2 - A(5)*A(6)*A(7) + A(9)*A(6)^2 + A(8)*A(7)^2 - 4*A(8)*A(9)*A(10)));
 
    if abs(x_offset)>1 || abs(y_offset)>1 || abs(z_offset)>1
        % adjusted peak falls outside set of 9 points fit,
        return % return absolute peak
    end
    
     % return only one-thousandth of a pixel precision
    % x_offset = round(1000*x_offset)/1000;
    % y_offset = round(1000*y_offset)/1000;    
    
    xpeak = xpeak + x_offset;
    ypeak = ypeak + y_offset; 
    zpeak = zpeak + z_offset;
    
    % Calculate extremum of fitted function
    max_f = [1 x_offset y_offset z_offset x_offset*y_offset x_offset*z_offset ...
        y_offset*z_offset x_offset^2 y_offset^2 z_offset^2] * A;
    
    
end
    
   