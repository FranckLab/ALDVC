function [Afidelity, Lreg] = regularizeNdMatrices(x, xGrid, smoothness, interpMethod)
% Returns the needed matrices used in the regulareNd problem.
%
%   [Afidelity, Lreg] = regularizeNdMatrices(x, xGrid)
%   [Afidelity, Lreg] = regularizeNdMatrices(x, xGrid, smoothness)
%   [Afidelity, Lreg] = regularizeNdMatrices(x, xGrid, smoothness, interpMethod)
%
%% Inputs
%      x - column vector or matrix of column vectors, containing scattered
%          data. Each row contains one point. Each column corresponds to a
%          dimension.
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
%          of the accuracty and shape that you are looking to obtain.
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
%
%% Output
%  Afidelity   - matirx. contains the fidelity equations. size(A,1) ==
%                size(x,1) == size(y,1). The number of rows in A
%                corresponds to the number of points in x,y. The number of
%                columns corresponds to the number points in the grid.
%
% Lreg - Cell array. L{i} corresponds to the scaled 2nd derivative
%        regularization of the ith dimension.
%% Description
% TODO Write a description. This function is useful if you want to solve a
% constrained version of regularizeNd problem.
%
%  For an introduction on how regularization works, start here:
%  https://mathformeremortals.wordpress.com/2013/01/29/introduction-to-regularizing-with-2d-data-part-1-of-3/
%
%% Acknowledgement
% Special thanks to Peter Goldstein, author of RegularizeData3D, for his
% coaching and help through writing regularizeNd.
%
%% Example
% TODO Write an example
%

% Author(s): Jason Nicholson
% $Revision: 1.0 $  $Date: 2018/01/12 16:48:00 $

%% Input Checking and Default Values
narginchk(2, 4);
nargoutchk(0,2);

% helper function used mostly for when variables are renamed
getname = @(x) inputname(1);

% Set default smoothness or check smoothness
if nargin() < 3 || isempty(smoothness)
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
if nargin() < 4 || isempty(interpMethod)
    interpMethod = 'linear';
else
    interpMethodsPossible = {'cubic', 'linear', 'nearest'};
    assert(any(strcmpi(interpMethod, interpMethodsPossible)), '%s is not a possible interpolation method. Check your spelling.', interpMethod);
    interpMethod = lower(interpMethod);
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

% Check y rows matches the number in x
nScatteredPoints = size(x,1);

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
        
        % Form the sparse A matrix for fidelity equations
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
        
        % Form the sparse A matrix for fidelity equations
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
            % operations at the cost of memory. There are 4 cubic lagrange
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

         % Form the sparse A matrix for fidelity equations
        Afidelity = sparse(repmat((1:nScatteredPoints)',1,4^nDimensions), xWeightIndex, weight, nScatteredPoints, nTotalGridPoints);
        
    otherwise
        error('Code should never reach this point. If it does, there is a bug.');
end

% clean up
clear(getname(dx), getname(weight), getname(x), getname(xWeightIndex));

%% Smoothness Equations
if nargout() > 1
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
            % smothing, index1, index2, and index3 do not exist. The clear
            % statement later would through an eror if index1, index2, and
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
 % array of ndgrid is same as this algoritm when dim = n. For instance, if
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


