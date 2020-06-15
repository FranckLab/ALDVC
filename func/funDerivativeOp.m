 function A = funDerivativeOp(M,N,h)
% Generate first order gradient operator A in 2D case: {F} = {A}{U}
 
DIM = 2;
A = sparse(DIM^2*M*N, DIM*M*N); sizeA = size(A); 


[XX,YY ] = ndgrid(1:M,1:N ); XX = XX(:); YY = YY(:);  
INDEXI = ones(2*DIM^2*M*N ,1); INDEXJ = INDEXI; INDEXVAL = INDEXI;

% ------ Inside assign values ------
 for tempi = 1:M*N 
    tempx = XX(tempi); tempy = YY(tempi);  
%     for tempy = 1:N
%         for tempx = 1:M
            
            % Determine BorderOrNot
            BorderOrNot = ones(1,DIM*2);  
            if tempy == 1, BorderOrNot(2) = 0; end; if tempy == N, BorderOrNot(1) = 0; end
            if tempx == 1, BorderOrNot(4) = 0; end; if tempx == M, BorderOrNot(3) = 0; end
            
            index = tempx + M*(tempy-1)  ; % find the index of the point (tempx,tempy,tempz);
            indexNeighbors = index*ones(1,DIM*2)+BorderOrNot.*[M,-M,1,-1 ]; %{Back,Front,Right,Left,Top,Bottom}
            % indexBack = index+M; indexFront = index-M; indexLeft = index-1;
            % indexRight = index+1; indexTop = index+M*N; indexBottom = index-M*N;
            
            % Find index of affine deformation gradient tensor {F}
            indexFrow = DIM^2*index*ones(1,DIM^2)+[-(DIM^2-1):1:0]; %{F11,F21,F31,F12,F22,F32,F13,F23,F33};
            
            % % Just some help lines
            % indexF11col1 = 3*indexLeft-2; indexF11col2 = 3*indexRight-2;
            % indexF12col1 = 3*indexFront-2; indexF12col2 = 3*indexBack-2;
            % indexF13col1 = 3*indexBottom-2; indexF13col2 = 3*indexTop-2;
            %
            % indexF21col1 = 3*indexLeft-1; indexF21col2 = 3*indexRight-1;
            % indexF22col1 = 3*indexFront-1; indexF22col2 = 3*indexBack-1;
            % indexF23col1 = 3*indexBottom-1; indexF23col2 = 3*indexTop-1;
            %
            % indexF31col1 = 3*indexLeft; indexF31col2 = 3*indexRight;
            % indexF32col1 = 3*indexFront; indexF32col2 = 3*indexBack;
            % indexF33col1 = 3*indexBottom; indexF33col2 = 3*indexTop;
            %
            % indexFcol1 = [3*indexLeft-2,3*indexLeft-1,3*indexLeft, ...
            % 3*indexFront-2,3*indexFront-1,3*indexFront,...
            % 3*indexBottom-2,3*indexBottom-1,3*indexBottom];
            
            indexFcol1 = [DIM*indexNeighbors(4)*ones(1,DIM)+[-(DIM-1):1:0], DIM*indexNeighbors(2)*ones(1,DIM)+[-(DIM-1):1:0] ];
            indexFcol2 = [DIM*indexNeighbors(3)*ones(1,DIM)+[-(DIM-1):1:0], DIM*indexNeighbors(1)*ones(1,DIM)+[-(DIM-1):1:0] ];
   
            % indexF1 = sub2ind(sizeA,indexFrow,indexFcol1);
            % indexF2 = sub2ind(sizeA,indexFrow,indexFcol2);
            
            INDEXI(2*4*tempi-7:2*4*tempi-0) = [indexFrow,indexFrow]';
            INDEXJ(2*4*tempi-7:2*4*tempi-0) = [indexFcol1,indexFcol2]';
            INDEXVAL(2*4*tempi-7:2*4*tempi-4) = -ones(4,1);
            
            if BorderOrNot(3)*BorderOrNot(4) == 0, INDEXVAL(2*4*tempi-7:2*4*tempi-6) = -2*ones(2,1); INDEXVAL(2*4*tempi-3:2*4*tempi-2) = 2*ones(2,1); end
            if BorderOrNot(1)*BorderOrNot(2) == 0, INDEXVAL(2*4*tempi-5:2*4*tempi-4) = -2*ones(2,1); INDEXVAL(2*4*tempi-1:2*4*tempi ) = 2*ones(2,1); end
             
            
            % A(indexF1) = -1; A(indexF2) = 1;
            % if BorderOrNot(3)*BorderOrNot(4) == 0, A(indexF1(1:DIM)) = -2; A(indexF2(1:DIM)) = 2; end
            %if BorderOrNot(1)*BorderOrNot(2) == 0, A(indexF1(DIM+1:2*DIM)) = -2; A(indexF2(DIM+1:2*DIM)) = 2; end
            
%         end
%     end
 end

A =  (1.0/(2*h)) * sparse(INDEXI,INDEXJ,INDEXVAL, sizeA(1),sizeA(2));


%% Old codes
% A = sparse(4*M*N, 2*M*N);
% A(1,1) = -2; A(1,3)= 2; 
% A(3,1) = -2; A(3,2*(M+1)-1) = 2;
% A(2,2) = -2; A(2,4) = 2;
% A(4,2) = -2; A(4,2*(M+1)) = 2;
% 
% % Inside
% for tempx = 2:(M-1)
%     for tempy = 2:(N-1)
%         index = tempx+M*(tempy-1); % Find the point position
%         indexUp = tempx+M*tempy;
%         indexDown = tempx+M*(tempy-2);
%         indexLeft = tempx-1+M*(tempy-1);
%         indexRight = tempx+1+M*(tempy-1);
%         
%         indexF11row = 4*index-3; 
%         indexF21row = 4*index-2;
%         indexF12row = 4*index-1;
%         indexF22row = 4*index;
%         
%         indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(indexRight)-1;
%         indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(indexUp)-1;
%         indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(indexRight);
%         indexF22col1 = 2*(indexDown); indexF22col2 = 2*(indexUp);
%         
%         A(indexF11row,indexF11col1) = -1; A(indexF11row,indexF11col2) = 1;
%         A(indexF21row,indexF21col1) = -1; A(indexF21row,indexF21col2) = 1;
%         A(indexF12row,indexF12col1) = -1; A(indexF12row,indexF12col2) = 1;
%         A(indexF22row,indexF22col1) = -1; A(indexF22row,indexF22col2) = 1; 
%     end
% end
% 
% % Bottom line
% for tempx = 2:(M-1)
%     for tempy = 1
%         index = tempx+M*(tempy-1); % Find the point position
%         indexUp = tempx+M*tempy;
%         % indexDown = tempx+M*(tempy-2);
%         indexLeft = tempx-1+M*(tempy-1);
%         indexRight = tempx+1+M*(tempy-1);
%         
%         indexF11row = 4*index-3; 
%         indexF21row = 4*index-2;
%         indexF12row = 4*index-1;
%         indexF22row = 4*index;
%         
%         indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(indexRight)-1;
%         indexF12col1 = 2*(index)-1; indexF12col2 = 2*(indexUp)-1;
%         indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(indexRight);
%         indexF22col1 = 2*(index); indexF22col2 = 2*(indexUp);
%         
%         A(indexF11row,indexF11col1) = -1; A(indexF11row,indexF11col2) = 1;
%         A(indexF21row,indexF21col1) = -1; A(indexF21row,indexF21col2) = 1;
%         A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
%         A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2; 
%         
%     end
% end
% 
% % Top
% for tempx = 2:(M-1)
%     for tempy = N
%         index = tempx+M*(tempy-1); % Find the point position
%         % indexUp = tempx+M*tempy;
%         indexDown = tempx+M*(tempy-2);
%         indexLeft = tempx-1+M*(tempy-1);
%         indexRight = tempx+1+M*(tempy-1);
%         
%         indexF11row = 4*index-3; 
%         indexF21row = 4*index-2;
%         indexF12row = 4*index-1;
%         indexF22row = 4*index;
%         
%         indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(indexRight)-1;
%         indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(index)-1;
%         indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(indexRight);
%         indexF22col1 = 2*(indexDown); indexF22col2 = 2*(index);
%         
%         A(indexF11row,indexF11col1) = -1; A(indexF11row,indexF11col2) = 1;
%         A(indexF21row,indexF21col1) = -1; A(indexF21row,indexF21col2) = 1;
%         A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
%         A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2; 
%     end
% end
% 
% % Left
% for tempx = 1
%     for tempy = 2:(N-1)
%         index = tempx+M*(tempy-1); % Find the point position
%         indexUp = tempx+M*tempy;
%         indexDown = tempx+M*(tempy-2);
%         % indexLeft = tempx-1+M*(tempy-1);
%         indexRight = tempx+1+M*(tempy-1);
%         
%         indexF11row = 4*index-3; 
%         indexF21row = 4*index-2;
%         indexF12row = 4*index-1;
%         indexF22row = 4*index;
%         
%         indexF11col1 = 2*(index )-1; indexF11col2 = 2*(indexRight)-1;
%         indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(indexUp)-1;
%         indexF21col1 = 2*(index ); indexF21col2 = 2*(indexRight);
%         indexF22col1 = 2*(indexDown); indexF22col2 = 2*(indexUp);
%         
%         A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
%         A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
%         A(indexF12row,indexF12col1) = -1; A(indexF12row,indexF12col2) = 1;
%         A(indexF22row,indexF22col1) = -1; A(indexF22row,indexF22col2) = 1; 
%     end
% end
% 
% % Right
% for tempx = M
%     for tempy = 2:(N-1)
%         index = tempx+M*(tempy-1); % Find the point position
%         indexUp = tempx+M*tempy;
%         indexDown = tempx+M*(tempy-2);
%         indexLeft = tempx-1+M*(tempy-1);
%         % indexRight = tempx+1+M*(tempy-1);
%         
%         indexF11row = 4*index-3; 
%         indexF21row = 4*index-2;
%         indexF12row = 4*index-1;
%         indexF22row = 4*index;
%         
%         indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(index )-1;
%         indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(indexUp)-1;
%         indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(index );
%         indexF22col1 = 2*(indexDown); indexF22col2 = 2*(indexUp);
%         
%         A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
%         A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
%         A(indexF12row,indexF12col1) = -1; A(indexF12row,indexF12col2) = 1;
%         A(indexF22row,indexF22col1) = -1; A(indexF22row,indexF22col2) = 1; 
%     end
% end
% 
% 
% % LeftBottom
% for tempx = 1
%     for tempy = 1
%         index = tempx+M*(tempy-1); % Find the point position
%         indexUp = tempx+M*tempy;
%         % indexDown = tempx+M*(tempy-2);
%         % indexLeft = tempx-1+M*(tempy-1);
%         indexRight = tempx+1+M*(tempy-1);
%         
%         indexF11row = 4*index-3; 
%         indexF21row = 4*index-2;
%         indexF12row = 4*index-1;
%         indexF22row = 4*index;
%         
%         indexF11col1 = 2*(index )-1; indexF11col2 = 2*(indexRight)-1;
%         indexF12col1 = 2*(index )-1; indexF12col2 = 2*(indexUp)-1;
%         indexF21col1 = 2*(index ); indexF21col2 = 2*(indexRight);
%         indexF22col1 = 2*(index ); indexF22col2 = 2*(indexUp);
%         
%         A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
%         A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
%         A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
%         A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2; 
%     end
% end
% 
% % RightBottom
% for tempx = M
%     for tempy = 1
%         index = tempx+M*(tempy-1); % Find the point position
%         indexUp = tempx+M*tempy;
%         %indexDown = tempx+M*(tempy-2);
%         indexLeft = tempx-1+M*(tempy-1);
%         %indexRight = tempx+1+M*(tempy-1);
%         
%         indexF11row = 4*index-3; 
%         indexF21row = 4*index-2;
%         indexF12row = 4*index-1;
%         indexF22row = 4*index;
%         
%         indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(index )-1;
%         indexF12col1 = 2*(index )-1; indexF12col2 = 2*(indexUp)-1;
%         indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(index );
%         indexF22col1 = 2*(index ); indexF22col2 = 2*(indexUp);
%         
%         A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
%         A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
%         A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
%         A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2; 
%     end
% end
% 
% % LeftTop
% for tempx = 1
%     for tempy = N
%         index = tempx+M*(tempy-1); % Find the point position
%         %indexUp = tempx+M*tempy;
%         indexDown = tempx+M*(tempy-2);
%         %indexLeft = tempx-1+M*(tempy-1);
%         indexRight = tempx+1+M*(tempy-1);
%         
%         indexF11row = 4*index-3; 
%         indexF21row = 4*index-2;
%         indexF12row = 4*index-1;
%         indexF22row = 4*index;
%         
%         indexF11col1 = 2*(index )-1; indexF11col2 = 2*(indexRight)-1;
%         indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(index )-1;
%         indexF21col1 = 2*(index ); indexF21col2 = 2*(indexRight);
%         indexF22col1 = 2*(indexDown); indexF22col2 = 2*(index );
%         
%         A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
%         A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
%         A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
%         A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2; 
%     end
% end
% 
% % RightTop
% for tempx = M
%     for tempy = N
%         index = tempx+M*(tempy-1); % Find the point position
%         %indexUp = tempx+M*tempy;
%         indexDown = tempx+M*(tempy-2);
%         indexLeft = tempx-1+M*(tempy-1);
%         %indexRight = tempx+1+M*(tempy-1);
%         
%         indexF11row = 4*index-3; 
%         indexF21row = 4*index-2;
%         indexF12row = 4*index-1;
%         indexF22row = 4*index;
%         
%         indexF11col1 = 2*(indexLeft)-1; indexF11col2 = 2*(index)-1;
%         indexF12col1 = 2*(indexDown)-1; indexF12col2 = 2*(index)-1;
%         indexF21col1 = 2*(indexLeft); indexF21col2 = 2*(index);
%         indexF22col1 = 2*(indexDown); indexF22col2 = 2*(index);
%         
%         A(indexF11row,indexF11col1) = -2; A(indexF11row,indexF11col2) = 2;
%         A(indexF21row,indexF21col1) = -2; A(indexF21row,indexF21col2) = 2;
%         A(indexF12row,indexF12col1) = -2; A(indexF12row,indexF12col2) = 2;
%         A(indexF22row,indexF22col1) = -2; A(indexF22row,indexF22col2) = 2; 
%     end
% end
% 
% 
% A = A/(2*h);