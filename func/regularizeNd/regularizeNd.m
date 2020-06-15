function yGrid = regularizeNd(x, y, xGrid, smoothness, interpMethod, solver, maxIterations, solverTolerance)
% regularizeNd  Fits a nD lookup table with smoothness to scattered data.
%
%   yGrid = regularizeNd(x, y, xGrid)
%   yGrid = regularizeNd(x, y, xGrid, smoothness)
%   yGrid = regularizeNd(x, y, xGrid, smoothness, interpMethod)
%   yGrid = regularizeNd(x, y, xGrid, smoothness, interpMethod, solver)
%   yGrid = regularizeNd(x, y, xGrid, smoothness, interpMethod, solver, maxIterations)
%   yGrid = regularizeNd(x, y, xGrid, smoothness, interpMethod, solver, maxIterations, solverTolerance)
%
%% Inputs
%      x - column vector or matrix of column vectors, containing scattered
%          data. Each row contains one point. Each column corresponds to a
%          dimension.
%
%      y - vector containing the corresponds values to x. y has the same
%          number of rows as x.
%
%  xGrid - cell array containing vectors defining the nodes in the grid in
%          each dimension. xGrid{1} corresponds with x(:,1) for instance.
%          Unequal spacing in the grid vectors is allowed. The grid vectors
%          must completely span x. For instance the values of x(:,1) must
%          be within the bounds of xGrid{1}. If xGrid does not span x,
%          an error is thrown. 
%
%  smoothness - scalar or vector. - The numerical "measure" of what we want
%          to achieve along an axis/dimension, regardless of the
%          resolution, the aspect ratio between axes, or the scale of the
%          overall problem. The ratio of smoothness to fidelity of the
%          output surface, i.e. ratio of smoothness to "goodness of
%          fit." This must be a positive real number. If it is a vector, it
%          must have same number of elements as columns in x.
%
%          A smoothness of 1 gives equal weight to fidelity (goodness of fit)
%          and smoothness of the output surface.  This results in noticeable
%          smoothing.  If your input data has little or no noise, use
%          0.01 to give smoothness 1% as much weight as goodness of fit.
%          0.1 applies a little bit of smoothing to the output surface.
%
%          If this parameter is a vector, then it defines the relative
%          smoothing to be associated with each axis/dimension. This allows
%          the user to apply a different amount of smoothing in the each
%          axis/dimension.
%
%          DEFAULT: 0.01
%
%   interpMethod - character, denotes the interpolation scheme used
%          to interpolate the data.
%
%          Even though there is a computational complexity difference between
%          linear, nearest, and cubic interpolation methods, the
%          interpolation method is not the dominant factor in the
%          calculation time in regularizeNd. The dominant factor in
%          calculation time is the size of the grid and the solver used. So
%          in general, do not choose your interpolation method based on
%          computational complexity. Choose your interpolation method because
%          of accuracy and shape that you are looking to obtain.
%
%          'linear' - Uses linear interpolation within the grid. linear
%                     interpolation requires that extrema occur at the grid
%                     points. linear should be smoother than nearest for
%                     the same grid. As the number of dimension grows,
%                     the number of grid points used to interpolate at a
%                     query point grows with 2^nDimensions. i.e. 2d needs 4
%                     points, 3d needs 8 points, 4d needs 16 points per
%                     query point. In general, linear can use smaller
%                     smoothness values than cubic and still be well
%                     conditioned.
%
%          'nearest' - Nearest neighbor interpolation. Nearest should
%                      be the least complex but least smooth.
%
%          'cubic' - Uses Lagrange cubic interpolation. Cubic interpolation
%                    allows extrema to occur at other locations besides the
%                    grid points. Cubic should provide the most flexible
%                    relationship for a given xGrid. As the number of
%                    dimension grows, the number of grid points used to
%                    interpolate at a query point grows with 4^nDimensions.
%                    i.e. 2d needs 16 points, 3d needs 64 points, 4d needs
%                    256 points per query point. cubic has good properties
%                    of accuracy and smoothness but is the most complex
%                    interpMethod to calculate.
%
%          DEFAULT: 'linear'
%
%
%   solver - string that denotes the solver used for the
%            resulting linear system. The default is most often the best
%            choice.
%
%          What solver should you use? The short answer is use 'normal' as
%          a first guess. '\' may be best numerically for most smoothness
%          parameters and high extents of extrapolation. If you receive
%          rank deficiency warnings with 'normal', try the '\' solver.
%          Otherwise, use the 'normal' solver because it is usually faster
%          than the '\' solver.
%
%          The larger the numbers of grid points, the larger the solve time.
%          Since the equations generated tends to be well conditioned, the
%          'normal' solver is  a good choice. Beware using 'normal' when a
%          small smoothing parameter is used, since this will make the
%          equations less well conditioned. The 'normal' solver for large
%          grids is 3x faster than the '\'.
%
%          Use the 'pcg', 'symmlq', or 'lsqr' solver when the 'normal' and
%          '\' fail. Out of memory errors with 'normal' or '\' are reason to
%          try the iterative solvers. These errors are rare however they
%          happen. Start with the 'pcg' solver. Then 'symmlq'. Finally try
%          'lsqr' solver. The 'lsqr' solver is usually slow compared to the
%          'pcg' and 'symmlq' solver.
%
%          '\' - uses matlab's backslash operator to solve the sparse
%                system.
%
%          'lsqr' - Uses the MATLAB lsqr solver. This solver is not
%                   recommended. Try 'pcg' or 'symmlq' first and use
%                   'lsqr' as a last resort. Experiments have shown that
%                   'pcg' and 'symmlq' solvers are faster and just as
%                   accurate as 'lsqr' for the matrices generated by
%                   regularizeNd. The same preconditioner as
%                   the 'pcg' solver is used.
%
%          'normal' - Constructs the normal equation and solves.
%                     x = (A'A)\(A'*y). From testing, this seems to be a well
%                     conditioned and faster way to solve this type of
%                     equation system than backslash x = A\y. Testing shows
%                     that the normal equation is 3x faster than the '\'
%                     solver for this type of problem. A'*A preserves the
%                     sparsity and is symmetric positive definite. Often
%                     A'*A will have less nonzero elements than A. i.e.
%                     nnz(A'*A) < nnz(A).
%                 
%          'pcg' - Calls the MATLAB pcg iterative solver that solves the
%                  normal equation, (A'A)*x = A'*y, for x. Use this solver
%                  first when 'normal' and '\' fail. The 'pcg' solver tries
%                  to generate the Incomplete Cholesky Factorization
%                  (ichol) as a preconditioner. If Incomplete Cholesky
%                  Factorization fails, then diagonal compensation is
%                  added. There may be a case where the preconditioner just
%                  cannot be calculated and thus no preconditioner is used.
%
%          'symmlq' - Calls the MATLAB symlq iterative solver that solves
%                     the normal equation, (A'A)*x = A'*y, for x. Use this
%                     solver if 'pcg' has issues. 'symmlq' uses the same
%                     preconditioner as 'pcg'.
%
%          DEFAULT: 'normal'
%
%
%   maxIterations - Only used if the solver is set to the iterative
%                   solvers, 'lsqr', 'pcg', or 'symmlq'. Reducing this will
%                   speed up the solver at the cost of accuracy. Increasing
%                   it will increase accuracy at the cost of time. The
%                   default value is the smaller of 100,000 and the number
%                   of nodes in the grid.
%
%          DEFAULT: min(1e5,  nTotalGridPoints)
%
%
%   solverTolerance - Only used if the solver is set to the iterative
%                     solvers, 'lsqr', 'pcg', or 'symmlq'. The
%                     solverTolerance is used with 'lsqr', 'pcg', or
%                     'symmlq'. Smaller increases accuracy and reduces
%                     speed. Larger decreases accuracy and increases speed.
%
%          DEFAULT: 1e-11*abs(max(y) - min(y))
%
%
%% Output
%  yGrid   - array containing the fitted surface or hypersurface
%            corresponding to the grid points xGrid. yGrid is in the ndgrid
%            format. In 2d, ndgrid format is the transpose of meshgrid
%            format.
%
%% Description
% regularizeNd answers the question what is the best possible lookup table
% that the scattered data input x and output y in the least squares sense
% with smoothing? regularizeNd is meant to calculate a smooth lookup table
% given n-D scattered data. regularizeNd supports extrapolation from a
% scattered data set as well.
%
% The calculated lookup table yGrid is meant to be used with
% griddedInterpolant class with the conservative memory form. Call
% griddedInterpolant like F = griddedInterpolant(xGrid, yGrid).
% 
% Desirable properties of regularizeNd
%     - Calculates a relationship between the input x and the output y
%       without definition of the functional form of x to y.
%     - Often the fit is superior to polynomial type fitting without 
%       the wiggles.
%     - Extrapolation is possible from a scattered data set. 
%     - After creating the lookup table yGrid and using it with
%       griddedInterpolant, as the query point moves away from the
%       scattered data, the relationship between the input x and output y
%       becomes more linear because of the smoothness equations and no
%       nearby fidelity equations. The linear relationship is a good
%       choice when the relationship between x and y is unknown in
%       extrapolation.
%     - regularizeNd can handle 1D, 2D, nD input data to 1D output data.
%        RegularizeData3D and gridfit can only handle 2D input and 1D out
%       (total 3D). 
%     - regularizeNd can handle setting the smoothness to 0 in any
%        axis/dimension. This means no smoothing is applied in a particular
%        axis/dimension and the data is just a least squares fit of a lookup
%        table in that axis/dimension.
%
%  For an introduction on how regularization works, start here:
%  https://mathformeremortals.wordpress.com/2013/01/29/introduction-to-regularizing-with-2d-data-part-1-of-3/
%
%% Acknowledgement
% Special thanks to Peter Goldstein, author of RegularizeData3D, for his
% coaching and help through writing regularizeNd.
%
%% Example
%
% % setup some input points, output points, and noise
% x = 0.5:0.1:4.5;
% y = 0.5:0.1:5.5;
% [xx,yy] = ndgrid(x,y);
% z = tanh(xx-3).*sin(2*pi/6*yy);
% noise = (rand(size(xx))-0.5).*xx.*yy/30;
% zNoise = z + noise;
% 
% % setup the grid for lookup table
% xGrid = linspace(0,6,210);
% yGrid = linspace(0,6.6,195);
% gridPoints = {xGrid, yGrid};
% 
% % setup some difference in scale between the different dimensions/axes
% xScale = 100;
% x = xScale*x;
% xx=xScale*xx;
% xGrid = xScale*xGrid;
% gridPoints{1} = xGrid; 
%
% % smoothness parameter. i.e. fit is weighted 1000 times greater than
% % smoothness.
% smoothness = 0.001;
% 
% % regularize
% zGrid = regularizeNd([xx(:), yy(:)], zNoise(:), gridPoints, smoothness);
% % Note this s the same as 
% % zGrid = regularizeNd([xx(:), yy(:)], zNoise(:), gridPoints, smoothness, 'linear', 'normal');
%
% % create girrdedInterpolant function
% F = griddedInterpolant(gridPoints, zGrid, 'linear');
% 
% % plot and compare
% surf(x,y,z', 'FaceColor', 'g')
% hold all;
% surf(x,y,zNoise','FaceColor', 'm')
% surf(xGrid, yGrid, zGrid', 'FaceColor', 'r')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% legend({'Exact', 'Noisy', 'regularizeNd'},'location', 'best');
%

% Author(s): Jason Nicholson
% $Revision: 1.6 $  $Date: 2017/12/10 18:04:00 $

%% Input Checking and Default Values
narginchk(3, 8);
nargoutchk(0,1);

% helper function used mostly for when variables are renamed
getname = @(x) inputname(1);

% Set default smoothness or check smoothness
if nargin() < 4 || isempty(smoothness)
    smoothness = 0.01;
else
    assert(all(smoothness>=0), '%s must be positive in all components.', getname(smoothness));
end

% check that the grid is of the right type
assert(iscell(xGrid), '%s must be a cell array.', getname(xGrid));

% calculate the number of dimension
nDimensions = size(x,2);

% check for the matching dimensionality
assert(nDimensions == numel(xGrid), 'Dimensionality mismatch. The number of columns in %s does not match the number of cells in %s.', getname(x), getname(xGrid));

% Check if smoothness is a scalar. If it is, convert it to a vector
if isscalar(smoothness)
    smoothness = ones(nDimensions,1).*smoothness;
end

% Set default interp method or check method
if nargin() < 5 || isempty(interpMethod)
    interpMethod = 'linear';
else
    interpMethodsPossible = {'cubic', 'linear', 'nearest'};
    assert(any(strcmpi(interpMethod, interpMethodsPossible)), '%s is not a possible interpolation method. Check your spelling.', interpMethod);
    interpMethod = lower(interpMethod);
end

% Set default solver or check the solver
if nargin() < 6 || (nargin()==6 && isempty(solver))
    solver = 'normal';
else
    directSolvers = {'\', 'normal'};
    iterativeSolvers = {'lsqr', 'pcg', 'symmlq'};
    solversPossible = horzcat(directSolvers, iterativeSolvers);
    assert(any(strcmpi(solver, solversPossible)), '%s is not an acceptable %s. Check spelling and try again.', solver, getname(solver));
end

% Check the grid vectors is a cell array
assert(iscell(xGrid), '%s must be a cell array where the ith cell contains nodes for ith dimension', getname(xGrid));

% arrange xGrid as a row cell array. This helps with cellfun and arrayfun
% later because the shape is always the same. From here on, the shape is
% known.
xGrid = reshape(xGrid,1,[]);

% arrange the grid vector as column vectors. This is helpful with arrayfun
% and cellfun calls because the shape is always the same. From here on the
% grid vectors shape is a known.
xGrid = cellfun(@(u) reshape(u,[],1), xGrid, 'UniformOutput', false);

% calculate the number of points in each dimension of the grid
nGrid = cellfun(@(u) length(u), xGrid);
nTotalGridPoints = prod(nGrid);

% check maxIterations if the solver is iterative
if nargin() == 6 && any(strcmpi(iterativeSolvers, solver))
    maxIterations = min(1e5, nTotalGridPoints);
elseif nargin() == 7 && any(strcmpi(iterativeSolvers, solver))
    message = sprintf('%s must be a positive scalar integer.', getname(maxIterations));
    assert(isscalar(maxIterations), message);
    assert(fix(maxIterations)==maxIterations, message);
    assert(maxIterations>0, message);
else
    % Do Nothing. maxIterations is unused.
end

% check solverTolerance if the solver is iterative
if any(nargin() == [6 7]) && any(strcmpi(iterativeSolvers, solver))
    solverTolerance = abs(max(y)-min(y))*1e-11;
elseif nargin() ==8 && any(strcmpi(iterativeSolvers, solver))
    message = sprintf('%s must be a positive scalar.', getname(solverTolerance));
    assert(isscalar(solverTolerance), message);
    assert(solverTolerance>0, message);
else
    % Do Nothing. maxIterations is unused.
end


% Check y rows matches the number in x
nScatteredPoints = size(x,1);
assert( nScatteredPoints == size(y, 1), '%s must have same number of rows as %s',getname(x), getname(y));



% Check input points are within min and max of grid.
xGridMin = cellfun(@(u) min(u), xGrid);
xGridMax = cellfun(@(u) max(u), xGrid);
assert(all(all(bsxfun(@ge, x, xGridMin))) & all(all(bsxfun(@le, x, xGridMax))), 'All %s points must be within the range of the grid vectors', getname(x));

% calculate the difference between grid points for each dimension
dx = cellfun(@(uGrid) diff(uGrid), xGrid, 'UniformOutput', false);

% Check for monotonic increasing grid points in each dimension
assert(all(cellfun(@(du) ~any(du<=0), dx)), 'All grid points in %s must be monotonically increasing.', getname(xGrid));

% Check that there are enough points to form an output surface. Linear and
% nearest interpolation types require 3 points in each output grid
% dimension because of the numerical 2nd derivative needs three points.
% Cubic interpolation requires 4 points per dimension.
switch interpMethod
    case 'linear'
        minGridVectorLength = 3;
    case 'nearest'
        minGridVectorLength = 3;
    case 'cubic'
        minGridVectorLength = 4;
    otherwise
        error('Code should never reach this otherwise there is a bug.')
end
assert(all(nGrid >= minGridVectorLength), 'Not enough grid points in each dimension. %s interpolation method and numerical 2nd derivatives requires %d points.', interpMethod, minGridVectorLength);
        

%% Calculate Fidelity Equations

switch interpMethod
    case 'nearest' % nearest neighbor interpolation in a cell
        
        % Preallocate before loop
        xWeightIndex = cell(1, nDimensions);
        
        for iDimension = 1:nDimensions
            % Find cell index
            % determine the cell the x-points lie in the xGrid
            % loop over the dimensions/columns, calculating cell index
            [~,xIndex] = histc(x(:,iDimension), xGrid{iDimension});
            
            % For points that lie ON the max value of xGrid{iDimension} (i.e. the
            % last value), histc returns an index that is equal to the length of
            % xGrid{iDimension}. xGrid{iDimension} describes nGrid(iDimension)-1
            % cells. Therefore, we need to find when a cell has an index equal to
            % the length of nGrid(iDimension) and reduce the index by 1.
            xIndex(xIndex == nGrid(iDimension))=nGrid(iDimension)-1;
            
            % Calculate the cell fraction. This corresponds to a value between 0 and 1.
            % 0 corresponds to the beginning of the cell. 1 corresponds to the end of
            % the cell. The min and max functions help ensure the output is always
            % between 0 and 1.
            cellFraction = min(1,max(0,(x(:,iDimension) - xGrid{iDimension}(xIndex))./dx{iDimension}(xIndex)));
            
            % calculate the index of nearest point
            xWeightIndex{iDimension} = round(cellFraction)+xIndex;
        end
         
        % clean up a little
        clear(getname(cellFraction), getname(xIndex));
        
        % calculate linear index
        xWeightIndex = subscript2index(nGrid, xWeightIndex{:});
        
        % the weight for nearest interpolation is just 1
        weight  = 1;
        
        % Form the sparse Afidelity matrix for fidelity equations
        Afidelity = sparse((1:nScatteredPoints)', xWeightIndex, weight, nScatteredPoints, nTotalGridPoints);
        
    case 'linear'  % linear interpolation
        
        % This will be needed below
        % Each cell has 2^nDimension nodes. The local dimension index label is 1 or 2 for each dimension. For instance, cells in 2d
        % have 4 nodes with the following indexes:
        % node label  =  1  2  3  4
        % index label = [1, 1, 2, 2;
        %                1, 2, 1, 2]
        % Said in words, node 1 is one, one. node 2 is one, two. node
        % three is two, one. node 4 is two, two.
        localCellIndex = (arrayfun(@(digit) str2double(digit), dec2bin(0:2^nDimensions-1))+1)';
        
        % preallocate
        weight = ones(nScatteredPoints, 2^nDimensions);
        xWeightIndex = cell(1, nDimensions);
        
        % loop over dimensions calculating subscript index in each dimension for
        % scattered points.
        for iDimension = 1:nDimensions
            % Find cell index
            % determine the cell the x-points lie in the xGrid
            % loop over the dimensions/columns, calculating cell index
            [~,xIndex] = histc(x(:,iDimension), xGrid{iDimension});
            
            % For points that lie ON the max value of xGrid{iDimension} (i.e. the
            % last value), histc returns an index that is equal to the length of
            % xGrid{iDimension}. xGrid{iDimension} describes nGrid(iDimension)-1
            % cells. Therefore, we need to find when a cell has an index equal to
            % the length of nGrid(iDimension) and reduce the index by 1.
            xIndex(xIndex == nGrid(iDimension))=nGrid(iDimension)-1;
            
            % Calculate the cell fraction. This corresponds to a value between 0 and 1.
            % 0 corresponds to the beginning of the cell. 1 corresponds to the end of
            % the cell. The min and max functions help ensure the output is always
            % between 0 and 1.
            cellFraction = min(1,max(0,(x(:,iDimension) - xGrid{iDimension}(xIndex))./dx{iDimension}(xIndex)));
            
            % In linear interpolation, there is two weights per dimension
            %                                weight 1      weight 2
            weightsCurrentDimension = [1-cellFraction, cellFraction];
            
            % Calculate weights
            % After the for loop finishes, the rows of weight sum to 1 as a check.
            % multiply the weights from each dimension
            weight = weight.*weightsCurrentDimension(:, localCellIndex(iDimension,:));
            
            % compute the index corresponding to the weight
            xWeightIndex{iDimension} = bsxfun(@plus, xIndex, localCellIndex(iDimension,:)-1);
        end
        
        % clean up a little
        clear(getname(cellFraction), getname(xIndex), getname(weightsCurrentDimension), getname(localCellIndex));
        
        % calculate linear index
        xWeightIndex = subscript2index(nGrid, xWeightIndex{:});
        
        % Form the sparse Afidelity matrix for fidelity equations
        Afidelity = sparse(repmat((1:nScatteredPoints)',1,2^nDimensions), xWeightIndex, weight, nScatteredPoints, nTotalGridPoints);
        
    case 'cubic'
        
        % This will be needed below.
        % Each cubic interpolation has 4^nDimension nodes. The local 
        % dimension index label is 1, 2, 3, or 4 for each dimension. For 
        % instance, cubic interpolation in 2d has 16 nodes with the 
        % following indexes:
        %    node label  =  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
        % localCellIndex = [1 1 1 1 2 2 2 2 3 3  3  3  4  4  4  4;
        %                   1 2 3 4 1 2 3 4 1 2  3  4  1  2  3  4]
        localCellIndex = (arrayfun(@(digit) str2double(digit), dec2base(0:4^nDimensions-1,4))+1)';
        
        % Preallocate before loop
        weight = ones(nScatteredPoints, 4^nDimensions);
        xWeightIndex = cell(1, nDimensions);
        
        for iDimension = 1:nDimensions
            % Find cell index. Determine the cell the x-points lie in the
            % current xGrid dimension.
            [~,xIndex] = histc(x(:,iDimension), xGrid{iDimension});
            
            % Calculate low index used in cubic interpolation. 4 points are
            % needed  for cubic interpolation. The low index corresponds to
            % the smallest grid point used in the interpolation. The min
            % and max ensures that the boundaries of the grid are
            % respected. For example, give a point x = 1.6 and a xGrid =
            % [0,1,2,3,4,5]. The points used for cubic interpolation would
            % be [0,1,2,3]. If x = 0.5, the points used would be [0,1,2,3];
            % this respects the bounds of the grid. If x = 4.9, the points
            % used would be [2,3,4,5]; again this respects the bounds of
            % the grid. 
            xIndex = min(max(xIndex-1,1), nGrid(iDimension)-3);
            
            % Setup to calculate the 1d weights in the current dimension.
            % The 1d weights are based on cubic Lagrange polynomial
            % interpolation. The alphas and betas below help keep the
            % calculation readable and also save on a few floating point
            % operations at the cost of memory. There are 4 cubic Lagrange
            % polynomials that correspond to the weights. They have the
            % following form
            %
            % p1(x) = (x-x2)/(x1-x2)*(x-x3)/(x1-x3)*(x-x4)/(x1-x4) 
            % p2(x) = (x-x1)/(x2-x1)*(x-x3)/(x2-x3)*(x-x4)/(x2-x4) 
            % p3(x) = (x-x1)/(x3-x1)*(x-x2)/(x3-x2)*(x-x4)/(x3-x4) 
            % p4(x) = (x-x1)/(x4-x1)*(x-x2)/(x4-x2)*(x-x3)/(x4-x3)
            % 
            % The alphas and betas are defined as follows
            % alpha1 = x - x1
            % alpha2 = x - x2
            % alpha3 = x - x3
            % alpha4 = x - x4
            %
            % beta12 = x1 - x2
            % beta13 = x1 - x3
            % beta14 = x1 - x4
            % beta23 = x2 - x3
            % beta24 = x2 - x4
            % beta34 = x3 - x4
            alpha1 = x(:,iDimension) - xGrid{iDimension}(xIndex);
            alpha2 = x(:,iDimension) - xGrid{iDimension}(xIndex+1);
            alpha3 = x(:,iDimension) - xGrid{iDimension}(xIndex+2);
            alpha4 = x(:,iDimension) - xGrid{iDimension}(xIndex+3);
            beta12 = xGrid{iDimension}(xIndex) - xGrid{iDimension}(xIndex+1);
            beta13 = xGrid{iDimension}(xIndex) - xGrid{iDimension}(xIndex+2);
            beta14 = xGrid{iDimension}(xIndex) - xGrid{iDimension}(xIndex+3);
            beta23 = xGrid{iDimension}(xIndex+1) - xGrid{iDimension}(xIndex+2);
            beta24 = xGrid{iDimension}(xIndex+1) - xGrid{iDimension}(xIndex+3);
            beta34 = xGrid{iDimension}(xIndex+2) - xGrid{iDimension}(xIndex+3);
            
            weightsCurrentDimension = [ alpha2./beta12.*alpha3./beta13.*alpha4./beta14, ...
                                       -alpha1./beta12.*alpha3./beta23.*alpha4./beta24, ...
                                        alpha1./beta13.*alpha2./beta23.*alpha4./beta34, ...
                                       -alpha1./beta14.*alpha2./beta24.*alpha3./beta34];

            % Accumulate the weight contribution for each dimension by
            % multiplication. After the for loop finishes, the rows of
            % weight sum to 1 as a check
            weight = weight.*weightsCurrentDimension(:, localCellIndex(iDimension,:));
            
            % compute the index corresponding to the weight
            xWeightIndex{iDimension} = bsxfun(@plus, xIndex, localCellIndex(iDimension,:)-1);
        end
        
        % clean up a little
        clear(getname(alpha1), getname(alpha2), getname(alpha3), getname(alpha4), ...
              getname(beta12), getname(beta13), getname(beta14), getname(beta23), ...
              getname(beta24), getname(beta34), getname(weightsCurrentDimension), ...
              getname(xIndex), getname(localCellIndex));
        
        % convert linear index
        xWeightIndex = subscript2index(nGrid, xWeightIndex{:});

         % Form the sparse Afidelity matrix for fidelity equations
        Afidelity = sparse(repmat((1:nScatteredPoints)',1,4^nDimensions), xWeightIndex, weight, nScatteredPoints, nTotalGridPoints);
        
    otherwise
        error('Code should never reach this point. If it does, there is a bug.');
end

% clean up
clear(getname(dx), getname(weight), getname(x), getname(xWeightIndex));

%% Smoothness Equations

%%% calculate the number of smoothness equations in each dimension

% nEquations is a square matrix where the ith row contains
% number of smoothing equations in each dimension. For instance, if the
% nGrid is [ 3 6 7 8] and ith row is 2, nEquationPerDimension contains
% [3 4 7 8]. Therefore, the nSmoothnessEquations is 3*4*7*8=672 for 2nd dimension (2nd row).
nEquationsPerDimension = repmat(nGrid, nDimensions,1);
nEquationsPerDimension = nEquationsPerDimension - 2*eye(nDimensions);
nSmoothnessEquations = prod(nEquationsPerDimension,2);

% Calculate the total number of Smooth equations
nTotalSmoothnessEquations = sum(nSmoothnessEquations);

%%% Calculate regularization matrices

% Preallocate the regularization equations
Lreg = cell(nDimensions, 1);

% compute the index multiplier for each dimension. This is used for
% calculating the linear index.
multiplier = cumprod(nGrid);

% loop over each dimension. calculate numerical 2nd derivatives weights.
for iDimension=1:nDimensions
    if smoothness(iDimension) == 0
        nTotalSmoothnessEquations = nTotalSmoothnessEquations - nSmoothnessEquations(iDimension);
        Lreg{iDimension} = [];
        
        % In the special case you try to fit a lookup table with no
        % smoothing, index1, index2, and index3 do not exist. The clear
        % statement later would throw an error if index1, index2, and
        % index3 did not exist.
        index1=[];
        index2=[];
        index3=[];
    else
        % initialize the index for the first grid vector
        if iDimension==1
            index1 = (1:nGrid(1)-2)';
            index2 = (2:nGrid(1)-1)';
            index3 = (3:nGrid(1))';
        else
            index1 = (1:nGrid(1))';
            index2 = index1;
            index3 = index1;
        end
        
        % loop over dimensions accumulating the contribution to the linear
        % index vector in each dimension. Note this section of code works very
        % similar to combining ndgrid and sub2ind. Basically, inspiration came
        % from looking at ndgrid and sub2ind.
        for iCell = 2:nDimensions
            if iCell == iDimension
                index1 = reshape(bsxfun(@(indx, currentIndex) indx + (currentIndex-1)*multiplier(iCell-1), index1, (1:nGrid(iCell)-2)), [], 1);
                index2 = reshape(bsxfun(@(indx, currentIndex) indx + (currentIndex-1)*multiplier(iCell-1), index2, (2:nGrid(iCell)-1)), [], 1);
                index3 = reshape(bsxfun(@(indx, currentIndex) indx + (currentIndex-1)*multiplier(iCell-1), index3, (3:nGrid(iCell))), [], 1);
            else
                currentDimensionIndex = 1:nGrid(iCell);
                index1 = reshape(bsxfun(@(indx, currentIndex) indx + (currentIndex-1)*multiplier(iCell-1), index1, currentDimensionIndex), [], 1);
                index2 = reshape(bsxfun(@(indx, currentIndex) indx + (currentIndex-1)*multiplier(iCell-1), index2, currentDimensionIndex), [], 1);
                index3 = reshape(bsxfun(@(indx, currentIndex) indx + (currentIndex-1)*multiplier(iCell-1), index3, currentDimensionIndex), [], 1);
            end
        end
        
        
        % Scales as if there is the same number of residuals along the
        % current dimension as there are fidelity equations total; use the
        % square root because the residuals will be squared to minimize
        % squared error.
        smoothnessScale = sqrt(nScatteredPoints/nSmoothnessEquations(iDimension));
        
        % Axis Scaling. This is equivalent to normalizing the current axis
        % to 0 to 1. i.e. If you scale one axis, the same smoothness factor
        % can be used to get similar shaped topology.
        axisScale = (xGridMax(iDimension) - xGridMin(iDimension)).^2;

        
        % Create the Lreg for each dimension and store it a cell array.
        Lreg{iDimension} = sparse(repmat((1:nSmoothnessEquations(iDimension))',1,3), ...
            [index1, index2, index3], ...
            smoothness(iDimension)*smoothnessScale*axisScale*secondDerivativeWeights(xGrid{iDimension},nGrid(iDimension), iDimension, nGrid), ...
            nSmoothnessEquations(iDimension), ...
            nTotalGridPoints);
    end
end

% clean up and free up memory
clear(getname(index1), getname(index2), getname(index3), getname(xGrid));

%% Assemble and Solve the Overall Equation System

% concatenate the fidelity equations and smoothing equations together
A = vertcat(Afidelity, Lreg{:});

% clean up
clear(getname(Afidelity), getname(Lreg)); 

% solve the full system
switch solver
    case '\'
            yGrid = A\[y;sparse(nTotalSmoothnessEquations,1)];
    case 'normal'
            yGrid = (A'*A)\(A'*[y;sparse(nTotalSmoothnessEquations,1)]);
    case {'lsqr', 'pcg', 'symmlq'}    
        switch solver
            case {'pcg', 'symmlq'}
                % setup needed normal equation matrices
                AA = A'*A;
                d = A'*[y;sparse(nTotalSmoothnessEquations,1)];
                
                % clean up
                clear(getname(A), getname(y));
                
                % calculate preconditioner if possible
                [M, preconditioner] = calculatePreconditioner(AA);
                
                % Call pcg or symmlq differently depending on the preconditioner
                switch preconditioner
                    case 'none'
                        [yGrid, solverExitFlag] = feval(solver, AA, d, solverTolerance, maxIterations);
                    case 'ichol'
                        [yGrid, solverExitFlag] = feval(solver, AA, d, solverTolerance, maxIterations, M, M');
                    otherwise
                        error('Code should never reach this. Something is wrong with the preconditioner switch statement. Fix it.');
                end % end pcg, symmlq preconditioner switch statement
                
            case 'lsqr'
                % calculate preconditioner if possible
                [M, preconditioner] = calculatePreconditioner(A'*A);
                
                % Call lsqr differently depending on the preconditioner
                switch preconditioner
                    case 'none'
                        [yGrid, solverExitFlag] = lsqr(A,[y;sparse(nTotalSmoothnessEquations,1)], solverTolerance, maxIterations);
                    case 'ichol'
                        [yGrid, solverExitFlag] = lsqr(A,[y;sparse(nTotalSmoothnessEquations,1)], solverTolerance, maxIterations, M');
                    otherwise
                        error('Code should never reach this. Something is wrong with the preconditioner switch statement. Fix it.');
                end % end lsqr preconditioner switch statement
                
            otherwise
                error('Code should never reach this. Something is wrong with iterative solver switch statement.');
        end % end iterative solver switch block
        
        % Check the iterative solver flag
        switch solverExitFlag
            case 0
                % Do nothing. This is good.
            case 1
                warning('%s iterated %d times but did not converge.', solver, maxIterations);
            case 2
                warning('The %s preconditioner was ill-conditioned.', solver);
            case 3
                warning('%s stagnated. (Two consecutive iterates were the same.)', solver);
            case 4
                warning('During %s solving, one of the scalar quantities calculated during pcg became too small or too large to continue computing.', solver);
            otherwise
                error('Code should never reach this. Something is wrong with iterative flag switch block.');
        end % iterative solver flag switch block
        
    otherwise
        error('Code should never reach this line. If it does, there is a bug.');
end  % switch solver

% convert to a full column vector
yGrid = full(yGrid);

% reshape if needed
if nDimensions > 1
    yGrid = reshape(yGrid, nGrid);
end

end %


%%
function weights = secondDerivativeWeights(x, nX, dim, arraySize)
% calculates the weights for a 2nd order numerical 2nd derivative
%
% Inputs
% x - grid vector
% nX - The length of x.
% dim - The dimension for which the numerical 2nd derivative is calculated
% arraySize - The size of the grid.
% Outputs
% weights  - weights of the numerical second derivative in a column vector
% form

% Calculate the numerical second derivative weights.
% The weights come from differentiating the parabolic Lagrange polynomial twice.
%
% parabolic Lagrange polynomial through 3 points:
% y = [(x-x2)*(x-x3)/((x1-x2)*(x1-x3)), (x-x1)*(x-x3)/((x2-x1)*(x2-x3)), (x-x1)*(x-x2)/((x3-x1)*(x3-x2))]*[y1;y2;y3];
%
% differentiating twice:
% y'' = 2./[(x1-x2)*(x1-x3), (x2-x1)*(x2-x3), (x3-x1)*(x3-x2)]*[y1;y2;y3];
%
x1 = x(1:nX-2);
x2 = x(2:nX-1);
x3 = x(3:nX);
weights = 2./[(x1-x3).*(x1-x2), (x2-x1).*(x2-x3), (x3-x1).*(x3-x2)];

% expand the weights across other dimensions and convert to  column vectors
weights = [reshape(ndGrid1D(weights(:,1), dim, arraySize),[], 1), ...
           reshape(ndGrid1D(weights(:,2), dim, arraySize),[], 1), ...
           reshape(ndGrid1D(weights(:,3), dim, arraySize),[], 1)];
end
 
 %%
 function xx = ndGrid1D(x, dim,  arraySize)
 % copies x along all dimensions except the dimension dim
 %
 % Inputs
 % x - column vector
 % dim - The dimension that x is not copied
 % arraySize - The size of the output array. arraySize(dim) is not used.
 %
 % Outputs
 % xx - array with size arraySize except for the dimension dim. The length
 % of dimension dim is numel(x).
 %
 % Description
 % This is very similar to ndgrid except that ndgrid returns all arrays for
 % each input vector. This algorithm returns only one array. The nth output
 % array of ndgrid is same as this algorithm when dim = n. For instance, if
 % ndgrid is given three input vectors, the output size will be arraySize.
 % Calling ndGrid1D(x,3, arraySize) will return the same values as the 3rd
 % output of ndgrid.
 %
 
 % reshape x into a vector with the proper dimensions. All dimensions are 1
 % expect the dimension dim.
 s = ones(1,length(arraySize));
 s(dim) = numel(x);
 if length(arraySize) == 1
     xx = x;
 else
     x = reshape(x,s);
     % expand x along all the dimensions except dim
     arraySize(dim) = 1;
     xx = repmat(x, arraySize);
 end % end if
 end % end function
 
 %%
function ndx = subscript2index(siz,varargin)
% Computes the linear index from the subscripts for an n dimensional array
%
% Inputs
% siz - The size of the array.
% varargin - has the same length as length(siz). Contains the subscript in
% each dimension.
% 
% Description
% This algorithm is very similar sub2ind. However, it will work for 1-D and
% all of the extra functionality for other data types is removed.

k = cumprod(siz);

%Compute linear indices
ndx = varargin{1};
for i = 2:length(varargin)
    ndx = ndx + (varargin{i}-1)*k(i-1);
end
end

%%
function [M, preconditioner] = calculatePreconditioner(AA)
% Calculate the incomplete Cholesky decomposition where M*M'~AA. Use
% diagonal compensation if the incomplete Cholesky decomposition does not
% exist for AA so that M*M'~AA + alpha*diag(diag(AA)). This algorithm should 
% be robust at always producing a preconditioner. 

% set preconditioner to none starting out
preconditioner = 'none';
M = [];

% Try to calculate the ichol preconditioner
try
    M = ichol(AA);
    preconditioner = 'ichol';
catch
    % initial calculations for diagonal compensation
    diagonalCompensation0 = full(max(sum(abs(AA),2)./diag(AA)));
    
    % check that the diagonal compensation is not nan or inf
    if ~isfinite(diagonalCompensation0)
        % Not possible to calculate diagonal compensation. Don't compute
        % preconditioner.
    else
        % find the bounds of a good diagonal compensation separated by a
        % factor of 10
        diagonalCompensationNew = diagonalCompensation0;
        diagonalCompensationFailure = [];
        diagonalCompensationSuccess = [];
        MAX_PRECONDITIONER_BOUNDS_RECALCULATIONS = 20;
        for iDiagonalCompensation=1:MAX_PRECONDITIONER_BOUNDS_RECALCULATIONS
            % Try to calculate ichol with diagonal compensation. If ichol
            % is successful, divide the diagonal compensation by 10 and try
            % again. If there is a failure, multiply 10 and try again.
            try
                M = ichol(AA, struct('diagcomp', diagonalCompensationNew));
                diagonalCompensationSuccess = diagonalCompensationNew;
                diagonalCompensationNew = diagonalCompensationNew/10;
            catch
                diagonalCompensationFailure = diagonalCompensationNew;
                diagonalCompensationNew = diagonalCompensationNew*10;
            end
            
            % Check whether we have diagonalCompensationFailure and
            % diagonalCompensationSuccess. Quit the loop if we have both.
            if ~isempty(diagonalCompensationFailure) && ~isempty(diagonalCompensationSuccess)
                break;
            end
        end % end for loop for diagonal compensation separated by a factor 10
        
        % Make sure we have a diagonalCompensationFailure and
        % diagonalCompensationSuccess. Only proceed if we have both.
        if ~isempty(diagonalCompensationFailure) && ~isempty(diagonalCompensationSuccess)
            % Use a binary search to better find a better diagonal compensation
            MAX_PRECONDITIONER_BINARY_RECALCULATIONS = 3;
            for iDiagonalCompensation=1:MAX_PRECONDITIONER_BINARY_RECALCULATIONS
                diagonalCompensationNew = (diagonalCompensationFailure + diagonalCompensationSuccess)/2;
                try
                    M = ichol(AA, struct('diagcomp', diagonalCompensationNew));
                    diagonalCompensationSuccess = diagonalCompensationNew;
                catch
                    diagonalCompensationFailure = diagonalCompensationNew;
                end
                
                % break the loop diagonalCompensationFailure and
                % diagonalCompensationSuccess are close together.
                if diagonalCompensationFailure + 1000*eps(diagonalCompensationFailure) > diagonalCompensationSuccess
                    break;
                end
            end % end binary search for loop
            
            % Make sure we have a preconditioner
            if ~isempty(M)
                preconditioner = 'ichol';
            else
                % Do nothing. No preconditioner.
            end % end check for preconditioner
        end % end check for diagonal compensation bounds
    end % end initial compensation guess
end % end ichol try-catch block
end % end calculatePreconditioner

