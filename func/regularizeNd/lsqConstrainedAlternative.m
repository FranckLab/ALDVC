function x = lsqConstrainedAlternative(C,d,AInequality,bInequality)
% Solve the least squares inequality constraint problem by reformulation to minimal distance problem
%
%   x = lsqConstrainedAlternative(C,d,AInequality,bInequality)
%
%% Inputs
% C - Multiplier matrix, specified as a matrix of doubles. C represents the multiplier of the solution x in the expression C*x - d. C is M-by-N, where M is the number of equations, and N is the number of elements of x.
% d - Constant vector, specified as a vector of doubles. d represents the additive constant term in the expression C*x - d. d is M-by-1, where M is the number of equations.
% AInequality - Linear inequality constraint matrix, specified as a matrix of doubles. AInequality represents the linear coefficients in the constraints AInequality*x ? bInequality. A has size Mineq-by-N, where Mineq is the number of constraints and N is the number of elements of x. To save memory, pass A as a sparse matrix.
% bInequality - Linear inequality constraint vector, specified as a vector of doubles. bInequality represents the constant vector in the constraints AInequality*x ? bInequality. b has length Mineq, where A is Mineq-by-N.
% 
%% Outputs
% x - Solution, returned as a vector that minimizes the norm of C*x-d subject to the constraints.
%
%% Description
% Solves:
% minimize 1/2*||C*x-d||_2^2 
%     x
%
% subject to AInequality*x <= bInequality
%
% lsqConstrainedAlternative reformulates a linear optimization problem with inequality constraint to a minimal distance problem and uses lsqnonneg to solve the
% problem. Minimal distance problem looks like this:
%    minimize ||x||^2
%       x
%
%     subject to   Abar*x <= bbar 
% 
% See the following link for discussion: https://www.mathworks.com/matlabcentral/answers/402953-reformulate-a-constrained-linear-least-square-problem?s_tid=prof_contriblnk
% For more information see: Lawson, C. L., and R. J. Hanson. "Solving Least Squares Problems, Classics in Applied Mathematics, SIAM, 1995."
%
%% Example
% C = [0.9501    0.7620    0.6153    0.4057
%     0.2311    0.4564    0.7919    0.9354
%     0.6068    0.0185    0.9218    0.9169
%     0.4859    0.8214    0.7382    0.4102
%     0.8912    0.4447    0.1762    0.8936];
% d = [0.0578
%     0.3528
%     0.8131
%     0.0098
%     0.1388];
% A = [0.2027    0.2721    0.7467    0.4659
%     0.1987    0.1988    0.4450    0.4186
%     0.6037    0.0152    0.9318    0.8462];
% b = [0.5251
%     0.2026
%     0.6721];
%
% x = lsqConstrainedAlternative(C, d, A, b)
% % Compare with lsqlin(C, d, A, b)
%

% Author(s): Jason Nicholson
% $Revision: 1.0 $  $Date: 2019/01/15 19:35:00 $

% Transform into minimal distance
[Q,R] = qr(C,0);
dbar = Q'*d;
Abar = AInequality/R;
bbar = bInequality - Abar*dbar;

% Get min-distance into lsqnonneg form 
n = size(Abar,2);
E = [Abar';
    bbar'];
f = [zeros(n,1); -1];
[~,~,residual] = lsqnonneg(E,f);
xbar = -residual(1:n)/residual(end);

% Map back
x = R\(xbar+dbar);


end

