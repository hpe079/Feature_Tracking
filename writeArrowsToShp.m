function writeArrowsToShp(filename, X, Y, U, V, S)
% writeArrowsToShp - write arrows to a shape file
%
% SYNTAX:
%
%   writeArrowsToShp(filename, X, Y, U, V)
% 
% INPUTS:
%
%   filename:    Name of shapefile to write to
%
%   X, Y, U, V:  Same as the parameters with corresponding names for 
%               the quiver function 
%
%   S:          A scaling factor.  Use 'auto' to automatically scale,
%               otherwise a number to multiply the speeds by to give the
%               arrow length in metres.  E.g. if this is 2 then an arrow
%               for speed 1.5 is drawn 3 m long, one for 2.4 is drawn 4.8
%               etc.  NOTE that this is NOT the same as the factor used in
%               quiver but should allow a consistent value to be used
%               across different outputs (e.g. different timeframes)
%

%----------------
% STEP 0: CHECK INPUTS
%----------------
if nargin < 5
    error('Must supply five inputs');
elseif nargin < 6
    S = 'auto';
end

i_checkInputs(filename, X, Y, U, V, S);

%----------------
% STEP 1: Work out scaling factor if needed
%----------------
if ischar(S)
    scaleFactor = i_calcAutoScale(X, Y, U, V)
else
    scaleFactor = S;
end

% And apply
U = U*scaleFactor;
V = V*scaleFactor;

%----------------
% STEP 2: Create the shape structure
%----------------

% Arrow coordinates and bounding box
valsX = cell(size(U));
valsY = valsX;
bbox = valsX;
geom = valsX;
id = valsX;

for i = 1:numel(valsX)
    [valsX{i}, valsY{i}] = i_makeArrow(X(i), Y(i), U(i), V(i));
    bbox{i} = [min(valsX{i}), min(valsY{i}); max(valsX{i}), max(valsY{i})];
    geom{i} = 'Line';
    id{i} = i;
end

sh = struct('Geometry', geom, ...
    'BoundingBox', bbox, ...
    'X', valsX, ...
    'Y', valsY, ...
    'ID', id);

%----------------
% STEP 3: Write the file
%----------------
shapewrite(sh, filename);

%--------------------------------------------------------------------------
function [arrowX, arrowY] = i_makeArrow(x, y, u, v)
% Make arrow...

% Easiest thing to do is make an arrow of the right length pointing north
% and then rotate it as needed...  So need to find the right angle and the
% length
arrLength = sqrt(u.^2 + v.^2);
arrAngle = atan2d(u, v);

% Arrow itself...
headLength = 0.3;
headWidth = 0.4;

arrowE = [0, 0, -headWidth*arrLength/2, 0, headWidth*arrLength/2, NaN];
arrowN = [0, arrLength, arrLength*(1-headLength), arrLength, arrLength*(1-headLength), NaN];

% And rotate...
arrowX = arrowE * cosd(arrAngle) + arrowN * sind(arrAngle);
arrowY = -arrowE * sind(arrAngle) + arrowN * cosd(arrAngle);

% And move...
arrowX = arrowX + x;
arrowY = arrowY + y;

%--------------------------------------------------------------------------
function factor = i_calcAutoScale(X, Y, U, V)
% Work out auto scaling for arrows.  Can't do exactly the same as quiver as
% not entirely clear what it does and low level code is pcoded

% Some measure of spacing of X and Y points?

% Start with the total area
dX = (max(X(:))-min(X(:)));
dY = (max(Y(:))-min(Y(:)));
totalArea = dY * dX;

% And so each point has this
areaPerPoint = totalArea / numel(X);

% And then split according to X/Y ratio
ratio = dY/dX;

% Basically solving:  xPerPoint * (xPerPoint * ratio) = areaPerPoint
xPerPoint = sqrt(areaPerPoint / ratio);
yPerPoint = areaPerPoint / xPerPoint;

% Find maximum ratio of U/X for both positive and negative directions?
maxURatio = max(abs(U(:)) / xPerPoint);
maxVRatio = max(abs(V(:)) / yPerPoint);

% And so...
factor = 1/max(maxURatio, maxVRatio);

% But just to be careful
if ~isfinite(factor)
    factor = 1;
end

%--------------------------------------------------------------------------
function i_checkInputs(filename, X, Y, U, V, S)
% Check inputs to the function

% Filename
if ~ischar(filename) || isempty(filename) || size(filename, 1) ~= 1
    error('First input (filename) must be a 1 x n string');
end

% X, Y, U and V must all be double and the same size?
if ~isa(X, 'double') 
    error('Second input (X) must be numeric');
elseif ~isa(Y, 'double') 
    error('Third input (Y) must be numeric');
elseif ~isa(U, 'double') 
    error('Fourth input (U) must be numeric');
elseif ~isa(V, 'double') 
    error('Fifth input (V) must be numeric');
elseif ~isequal(size(X), size(Y)) || ~isequal(size(X), size(U)) || ~isequal(size(X), size(V))
    error('Inputs two to five (X, Y, U and V) must all be the same size');
end


% And S must be scalar and non-negative
if ischar(S) && size(S, 1) == 1 && strcmpi(S, 'auto')
    % This is ok...
elseif ~isnumeric(S) || numel(S) ~= 1 || ~isfinite(S) || S <= 0
    error('Sixth input (S) must be ''auto'' or a positive value');
end
