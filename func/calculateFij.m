% Computes the deformation gradient tensor (Fij) and Jacobian of the
% deformation gradient tensor
%
%
%
% INPUTS : 
% -------------------------------------------------------------------------
%
%   u_         = displacement field vector with rigid drift removed. 
%                       Format: cell array, each containing a 3D matrix for each
%                       time point (components in x,y,z)
%                           unew{time}{1} = displacement in x-direction
%                           unew{time}{2} = displacement in y-direction
%                           unew{time}{3} = displacement in z-direction
%                           unew{time}{4} = magnitude
%   dm         = meshgrid spacing
%   m2vx       = meter to pixel conversion of original images in [x y z]
%   type       = spatial differentiation kernel used for gradientN
%                options: 'fb', 'prewitt', 'sobel', 'scharr',
%                       'stencil', or 'optimal#' (# are odd numbers from 5
%                       to 19)  
%                Suggested: 'optimal 9' 
%
%
% OUTPUTS
% -------------------------------------------------------------------------
%   Fij     = deformation gradient tensor calculated on the mesh grid
%               with spacing dm.
%               Format: Fij{time}{3x3 deformation gradient tensor matrix}
%   J       = determinant of the deformation gradient tensor defined on the
%               mesh grid with spacing dm.
%               Format: J{time}
%
%
% NOTES
% -------------------------------------------------------------------------
% none
% 
%
%%
function [Fij, J] = calculateFij(varargin)

%Establish variables and inputs
[u,spacing,m2vx,type] = parseInputs(varargin{:});
maxTime = length(u);    
Fij = cell(maxTime,1);  
J = cell(maxTime,1);    

for i = 1:maxTime
    Fij{i} = funCalculateFij(u,spacing,m2vx,type);
    J{i} = funCalculateJ(Fij{i});
end

end

function Fij = funCalculateFij(u,spacing,m2vx,type)

% Calculate Displacment Gradient
Fij = cell(3,3);
for i = 1:3,
    [Fij{i,1}, Fij{i,2}, Fij{i,3}] = gradientN(u{i}{1}*m2vx(i),type);
end

%Calculate Deformation Gradient
for i = 1:3
    for j = 1:3
        try
        Fij{i,j} = Fij{i,j}/(spacing(j)*m2vx(j));
        catch
            Fij{i,j} = Fij{i,j}/(spacing*m2vx(j));
        end
        
    end
end

% for k = 1:9, Fij{k} = Fij{k}/spacing; end
for i = 1:3, Fij{i,i} = Fij{i,i} + 1; end

end

function J = funCalculateJ(Fij)
%Calculate Jacobian of Deformation gradient
J =     Fij{1,1}.*Fij{2,2}.*Fij{3,3};
J = J + Fij{1,2}.*Fij{2,3}.*Fij{3,1};
J = J + Fij{1,3}.*Fij{2,1}.*Fij{3,2};
J = J - Fij{1,3}.*Fij{2,2}.*Fij{3,1};
J = J - Fij{1,2}.*Fij{2,1}.*Fij{3,3};
J = J - Fij{1,1}.*Fij{2,3}.*Fij{3,2};

negJ = (J<0);
if sum(negJ(:))>0
       J(J<0) = NaN;
       J = inpaint_nans3(J, 1);
       J(J<0) = NaN;
       J = fillmissing(J, 'nearest');
       
end


end

function [u,spacing,m2vx,type] = parseInputs(varargin)
u = varargin{1};    
spacing = varargin{2};
m2vx = varargin{3};

if length(varargin) < 4, type = 'optimal9';
else type = varargin{4};
end

for i = 1:length(u)
    u{i} = cellfun(@double, u{i}, 'UniformOutput', 0);
end

end