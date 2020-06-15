function [UNew,DUDX,DUDY,DUDZ] = PlaneFit3(U,winsize,Rad,index)

[M,N,L]=size(U); DUDY=U; DUDX=U; DUDZ=U; UNew=U; DIM=3;
 
%hbar = parfor_progressbar(length(tempii),'Please wait for PlaneFit2!');
hbar=waitbar(0,'wait for PlaneFit3!');

for kk = (Rad(3)+1):L-Rad(3)
    for jj = (Rad(2)+1):N-Rad(2)
      for ii = (Rad(1)+1):M-Rad(1)
        
        LSMatrix = ones((2*Rad(1)+1)*(2*Rad(2)+1)*(2*Rad(3)+1), 1+DIM);
        for tempk = -Rad(3):Rad(3)
        for tempj = -Rad(2):Rad(2)
            for tempi = -Rad(1):Rad(1)
                ind = (2*Rad(1)+1)*(2*Rad(2)+1)*(tempk+Rad(3)) + (2*Rad(1)+1)*(tempj+Rad(2)) + tempi+Rad(1) + 1;
                LSMatrix(ind,:) = [1 tempi*winsize(1) tempj*winsize(2) tempk*winsize(3)];
            end
        end
        end

        LSb = zeros((2*Rad(1)+1)*(2*Rad(2)+1)*(2*Rad(3)+1), 1);
        for tempk = -Rad(3):Rad(3)
        for tempj = -Rad(2):Rad(2)
            for tempi = -Rad(1):Rad(1)
                ind = (2*Rad(1)+1)*(2*Rad(2)+1)*(tempk+Rad(3)) + (2*Rad(1)+1)*(tempj+Rad(2)) + tempi+Rad(1) + 1;
                LSb(ind,:) = U(ii+tempi, jj+tempj, kk+tempk);
            end
        end
        end
         
        tempVector = (LSMatrix'*LSMatrix)\(LSMatrix'*LSb);
        UNew(ii,jj,kk) = tempVector(1);
        DUDX(ii,jj,kk) = tempVector(2);
        DUDY(ii,jj,kk) = tempVector(3);
        DUDZ(ii,jj,kk) = tempVector(4);
         
        % end
        % hbar.iterate(1);
        tempij = (kk-Rad(3)-1)*(M-2*Rad(1))*(N-2*Rad(2)) + (jj-Rad(2)-1)*(M-2*Rad(1)) + ii-Rad(1) ;
        if index == 1
            waitbar(  tempij/((M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)))/3  );
        elseif index == 2
            waitbar(  tempij/((M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)))/3 + 1/3  );
        elseif index == 3
            waitbar(  tempij/((M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)))/3 + 2/3 );
        else
            waitbar(  tempij/((M-2*Rad(1))*(N-2*Rad(2))*(L-2*Rad(3)))  );
        end
        
      end
    end
end
close(hbar);


