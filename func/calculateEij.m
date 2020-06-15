% Computes the lagrangian and infinitesimal strain tensors from the deformation 
% gradient tensor. 
%
% INPUTS : 
% -------------------------------------------------------------------------
%   Fij     = deformation gradient tensor calculated on the mesh grid
%               with spacing dm.
%               Format: Fij{time}{3x3 deformation gradient tensor matrix}
%                     
%
% OUTPUTS
% -------------------------------------------------------------------------
%   Eij     = Lagrangian strain tensor calculated on the mesh grid with
%               spacing dm.
%                   Format: Eij{time}{3x3 Lagrangian strain tensor
%                   matrix}
%   eij     = Infinitesimal strain tensor calculated on the mesh grid with
%               spacing dm. 
%                   Format: eij{time}{3x3 infinitesimal strain tensor
%                   matrix}
%
%
% NOTES
% -------------------------------------------------------------------------
% none
%%

function [Eij, eij] = calculateEij(Fij)

%Set up variables
maxTime = length(Fij);
Eij = cell(maxTime,1);
eij = cell(maxTime,1);

for i = 1:maxTime
  Eij{i} = funCalculateLagrangianEij(Fij{i});
  eij{i} = funCalculateEij(Fij{i});
end


end

function eij = funCalculateEij(Fij)
%Calculate Infinitesimal strain
% A Bower, "solid mechanics", pg 22
eij = cell(3,3);

eij{1,1} = Fij{1,1} - 1;
eij{2,2} = Fij{2,2} - 1;
eij{3,3} = Fij{3,3} - 1;

eij{1,2} = 0.5*(Fij{1,2} + Fij{2,1});
eij{1,3} = 0.5*(Fij{1,3} + Fij{3,1});
eij{2,3} = 0.5*(Fij{2,3} + Fij{3,2});

eij{2,1} = eij{1,2}; eij{3,1} = eij{1,3}; eij{3,2} = eij{2,3};

end

function Eij = funCalculateLagrangianEij(Fij)
%Calcualte Lagrangian Strain
% A. Bower,"Solid mechanics", pg 20

Eij{1,1} = 0.5*(Fij{1,1}.*Fij{1,1} + Fij{2,1}.*Fij{2,1} + Fij{3,1}.*Fij{3,1} - 1);
Eij{2,2} = 0.5*(Fij{1,2}.*Fij{1,2} + Fij{2,2}.*Fij{2,2} + Fij{3,2}.*Fij{3,2} - 1);
Eij{3,3} = 0.5*(Fij{1,3}.*Fij{1,3} + Fij{2,3}.*Fij{2,3} + Fij{3,3}.*Fij{3,3} - 1);
Eij{1,2} = 0.5*(Fij{1,1}.*Fij{1,2} + Fij{2,1}.*Fij{2,2} + Fij{3,1}.*Fij{3,2});
Eij{1,3} = 0.5*(Fij{1,1}.*Fij{1,3} + Fij{2,1}.*Fij{2,3} + Fij{3,1}.*Fij{3,3});
Eij{2,3} = 0.5*(Fij{1,2}.*Fij{1,3} + Fij{2,2}.*Fij{2,3} + Fij{3,2}.*Fij{3,3});

Eij{3,2} = Eij{2,3}; Eij{2,1} = Eij{1,2}; Eij{3,1} = Eij{1,3};

end



