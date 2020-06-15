function A = funDerivativeOp3(M,N,L,h)
% Generate first order gradient operator A in 3D case: {F} = {A}{U}

% A = coder.ceval('funDerivativeOp3',M,N,L,h);

DIM = 3;
A = sparse(DIM^2*M*N*L, DIM*M*N*L); sizeA = size(A); 
  
% ------ Inside assign values ------
[XX,YY,ZZ] = ndgrid(1:M,1:N,1:L); XX = XX(:); YY = YY(:); ZZ = ZZ(:);
INDEXI = zeros(2*DIM^2*M*N*L,1); INDEXJ = INDEXI; INDEXVAL = INDEXI;
for tempi = 1:M*N*L
    tempx = XX(tempi); tempy = YY(tempi); tempz = ZZ(tempi);
% for tempz = 1:L
%     for tempy = 1:N
%         for tempx = 1:M
            
            % Determine BorderOrNot
            BorderOrNot = ones(1,DIM*2); %{Back,Front,Right,Left,Top,Bottom}
            if tempz == 1, BorderOrNot(6) = 0; end; if tempz == L, BorderOrNot(5) = 0; end
            if tempy == 1, BorderOrNot(2) = 0; end; if tempy == N, BorderOrNot(1) = 0; end
            if tempx == 1, BorderOrNot(4) = 0; end; if tempx == M, BorderOrNot(3) = 0; end
            
            index = tempx + M*(tempy-1) + M*N*(tempz-1); % find the index of the point (tempx,tempy,tempz);
            indexNeighbors = index*ones(1,DIM*2)+BorderOrNot.*[M,-M,1,-1,M*N,-M*N]; %{Back,Front,Right,Left,Top,Bottom}
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
            indexFcol1 = [DIM*indexNeighbors(4)*ones(1,DIM)+[-(DIM-1):1:0], DIM*indexNeighbors(2)*ones(1,DIM)+[-(DIM-1):1:0], DIM*indexNeighbors(6)*ones(1,DIM)+[-(DIM-1):1:0]];
            indexFcol2 = [DIM*indexNeighbors(3)*ones(1,DIM)+[-(DIM-1):1:0], DIM*indexNeighbors(1)*ones(1,DIM)+[-(DIM-1):1:0], DIM*indexNeighbors(5)*ones(1,DIM)+[-(DIM-1):1:0]];
            
            %indexF1 = sub2ind(sizeA,indexFrow,indexFcol1);
            %indexF2 = sub2ind(sizeA,indexFrow,indexFcol2);
             
            INDEXI(2*9*tempi-17:2*9*tempi-0) = [indexFrow,indexFrow]';
            INDEXJ(2*9*tempi-17:2*9*tempi-0) = [indexFcol1,indexFcol2]';
            INDEXVAL(2*9*tempi-17:2*9*tempi-9) = -ones(9,1);
            INDEXVAL(2*9*tempi-8:2*9*tempi) = ones(9,1);
            
            if BorderOrNot(3)*BorderOrNot(4) == 0, INDEXVAL(2*9*tempi-17:2*9*tempi-15) = -2*ones(3,1); INDEXVAL(2*9*tempi-8:2*9*tempi-6) = 2*ones(3,1); end
            if BorderOrNot(1)*BorderOrNot(2) == 0, INDEXVAL(2*9*tempi-14:2*9*tempi-12) = -2*ones(3,1); INDEXVAL(2*9*tempi-5:2*9*tempi-3) = 2*ones(3,1); end
            if BorderOrNot(5)*BorderOrNot(6) == 0, INDEXVAL(2*9*tempi-11:2*9*tempi-9) = -2*ones(3,1); INDEXVAL(2*9*tempi-2:2*9*tempi ) = 2*ones(3,1); end
            
            % A(indexF1) = -1; A(indexF2) = 1;
            %if BorderOrNot(3)*BorderOrNot(4) == 0, A(indexF1(1:3)) = -2; A(indexF2(1:3)) = 2; end
            %if BorderOrNot(1)*BorderOrNot(2) == 0, A(indexF1(4:6)) = -2; A(indexF2(4:6)) = 2; end
            %if BorderOrNot(5)*BorderOrNot(6) == 0, A(indexF1(7:9)) = -2; A(indexF2(7:9)) = 2; end
%         end
%     end
% end
end

A =   sparse(INDEXI,INDEXJ,INDEXVAL, sizeA(1),sizeA(2));
A(1:9:end,:) = A(1:9:end,:)*(1.0/(2*h(1))) ;
A(2:9:end,:) = A(2:9:end,:)*(1.0/(2*h(1))) ;
A(3:9:end,:) = A(3:9:end,:)*(1.0/(2*h(1))) ;
A(4:9:end,:) = A(4:9:end,:)*(1.0/(2*h(2))) ;
A(5:9:end,:) = A(5:9:end,:)*(1.0/(2*h(2))) ;
A(6:9:end,:) = A(6:9:end,:)*(1.0/(2*h(2))) ;
A(7:9:end,:) = A(7:9:end,:)*(1.0/(2*h(3))) ;
A(8:9:end,:) = A(8:9:end,:)*(1.0/(2*h(3))) ;
A(9:9:end,:) = A(9:9:end,:)*(1.0/(2*h(3))) ;


%% ------ Bottom plane ------
% tempz = 1;
% for tempy = 2:(N-1)
%     for tempx = 2:(M-1)
%         
%         index = tempx + M*(tempy-1) + M*N*(tempz-1); % find the index of the point (tempx,tempy,tempz);
%         indexNeighbors = index*ones(1,6)+[M,-M,1,-1,M*N,0]; % Change Bottom to be index itself
%         indexFrow = 9*index*ones(1,9)+[-8:1:0]; % same with inside
%         indexFcol1 = [3*indexNeighbors(4)*ones(1,3)+[-2:1:0], 3*indexNeighbors(2)*ones(1,3)+[-2:1:0], 3*indexNeighbors(6)*ones(1,3)+[-2:1:0]];
%         indexFcol2 = [3*indexNeighbors(3)*ones(1,3)+[-2:1:0], 3*indexNeighbors(1)*ones(1,3)+[-2:1:0], 3*indexNeighbors(5)*ones(1,3)+[-2:1:0]];
%             
%         indexF1 = sub2ind(sizeA,indexFrow,indexFcol1);
%         indexF2 = sub2ind(sizeA,indexFrow,indexFcol2);
%         A(indexF1) = -1; A(indexF2) = 1; A(indexF1(7:9)) = -2; A(indexF2(7:9)) = 2;
%         
%     end
% end
% 
% % ------ Top plane ------
% % tempz = L;
% for tempy = 2:(N-1)
%     for tempx = 2:(M-1)
%         
%         index = tempx + M*(tempy-1) + M*N*(tempz-1); % find the index of the point (tempx,tempy,tempz);
%         indexNeighbors = index*ones(1,6)+[M,-M,1,-1,0,-M*N]; % Change Bottom to be index itself
%         indexFrow = 9*index*ones(1,9)+[-8:1:0]; % same with inside
%         indexFcol1 = [3*indexNeighbors(4)*ones(1,3)+[-2:1:0], 3*indexNeighbors(2)*ones(1,3)+[-2:1:0], 3*indexNeighbors(6)*ones(1,3)+[-2:1:0]];
%         indexFcol2 = [3*indexNeighbors(3)*ones(1,3)+[-2:1:0], 3*indexNeighbors(1)*ones(1,3)+[-2:1:0], 3*indexNeighbors(5)*ones(1,3)+[-2:1:0]];
%             
%         indexF1 = sub2ind(sizeA,indexFrow,indexFcol1);
%         indexF2 = sub2ind(sizeA,indexFrow,indexFcol2);
%         A(indexF1) = -1; A(indexF2) = 1; A(indexF1(7:9)) = -2; A(indexF2(7:9)) = 2;
%         
%     end
% end

