function [A, iA, B, iB] = featureTrackingReadImages(fileA, fileB)
% featureTrackingReadImages - read two geotiff images for feature tracking
% 
% SYNTAX:
%
%   [A, Ia, B, Ib] = featureTrackingReadImages(fileA, fileB)
% 
% INPUTS:
%
%   fileA:  Name of first file to input
%
%   fileB:  Name of second file to input
%
% OUTPUTS:
%
%   A:  The image data for the first image (as a TWO dimensional single 
%       precision numeric array)
%
%   Ia: geotiff information corresponding to the first image
%
%   B:  The image data for the second image (as a TWO dimensional single 
%       precision numeric array)
%
%   Ib: geotiff information corresponding to the second image
%
% NOTES: 
%
%   In addition to loading the images, the following checks are also
%   carried out:
%
%       i) That the geographical coordinates are in METRES rather than
%          latitude and longitude
%
%       ii) The two images have almost exactly the same resolution and take
%           up the same space
%
% EXAMPLE:
%
%   % Read in files
%   fA = 'G:\Isung_2015to2016\Isung290715_DEM_1.35m_clipRep.tif';
%   fB = 'G:\Isung_2015to2016\Isung120716_DEM_1.35m_clipRep.tif';
%   [A, Ia, B, Ib] = featureTrackingReadImages(fA, fB);
%   
%   % And display
%   figure; 
%   ax1 = subplot (1,2,1);
%   image (A);
%   ax2 = subplot (1,2,2);
%   image (B); 
%   linkaxes ([ax1 ax2]) 
%

%----------------
% STEP 0:  Check inputs
%----------------
if nargin < 2
    error('featureTrackingReadImages:notEnoughInputs', 'Must supply two inputs');
end
i_checkInputs(fileA, fileB);

%----------------
% STEP 1:  Load Files
%----------------

% Pass on to load and do a few quick checks
[A, iA] = i_loadAndCheckFile(fileA);
[B, iB] = i_loadAndCheckFile(fileB);

% And some cross checks?
i_doCrossChecks(iA, iB);

%--------------------------------------------------------------------------
function i_doCrossChecks(iA, iB)
% Cross checks between the two files

% Check A:  EXACTLY the same size??
if ~isequal(iA.RasterSize, iB.RasterSize)
    error('featureTrackingReadImages:differentSizes', 'Images are different sizes (%s and %s)', mat2str(iA.RasterSize), mat2str(iB.RasterSize));
end

% Check B:  Same resolution or just check the overall extent of the image?
dX = iA.XWorldLimits - iB.XWorldLimits;
dY = iA.YWorldLimits - iB.YWorldLimits;

% And some strings
if dX(1) <= 0
    leftStr = sprintf('Left edge of second image is %sm to the right of left edge of first image', num2str(-dX(1)));
else
    leftStr = sprintf('Left edge of second image is %sm to the left of left edge of first image', num2str(dX(1)));
end

if dX(2) <= 0
    rightStr = sprintf('Right edge of second image is %sm to the right of right edge of first image', num2str(-dX(2)));
else
    rightStr = sprintf('Right edge of second image is %sm to the left of right edge of first image', num2str(dX(2)));
end

if dY(1) <= 0
    bottomStr = sprintf('Bottom edge of second image is %sm above bottom edge of first image', num2str(-dY(1)));
else
    bottomStr = sprintf('Bottom edge of second image is %sm below bottom edge of first image', num2str(dY(1)));
end

if dY(2) <= 0
    topStr = sprintf('Top edge of second image is %sm above top edge of first image', num2str(-dY(2)));
else
    topStr = sprintf('Top edge of second image is %sm below top edge of first image', num2str(dY(2)));
end

% Are they errors or displays?!
maxTol = 1;
if dX(1) > maxTol
    error(leftStr);
else
    disp(leftStr);
end
if dX(2) > maxTol
    error(rightStr);
else
    disp(rightStr);
end
if dY(1) > maxTol
    error(bottomStr);
else
    disp(bottomStr);
end
if dY(2) > maxTol
    error(topStr);
else
    disp(topStr);
end

%--------------------------------------------------------------------------
function [A, iA] = i_loadAndCheckFile(fileA)
% Load and check an individual file

%-----------------
% STEP 1:  Load file
%-----------------

% Just pass on, but rewrap the error
try
    [A, iA] = geotiffread(fileA);
catch lErr
    error('featureTrackingReadImages:errorLoadingFile', 'Unexpected error loading "%s":  %s', fileA, lErr.message);
end

%-----------------
% STEP 2:  Make consistent
%-----------------

% FOR THE MOMENT JUST SUPPORTING TWO DIMENSIONAL ARRAYS THAT ARE ALREADY
% SINGLE
if isa(A, 'uint16') && size(A, 3) == 1
    A = single(A)/65535;
elseif ~isa(A, 'single') || size(A, 3) ~= 1
    error('featureTrackingReadImages:unsupportedImageType', 'Image not yet supported - please see Dave to make necessary changes!! (File: "%s")', fileA);
end

% IN FUTURE WILL PROBABLY NEED TO ADD UINT8 OR UINT16 FOR NON-DEM IMAGES AS
% WELL AS [R, G, B] IMAGES

% AND check for any missing values?
imInfo = imfinfo(fileA);
if isfield(imInfo, 'GDAL_NODATA') && ischar(imInfo.GDAL_NODATA)
    missingVal = str2double(imInfo.GDAL_NODATA);
    A(A == missingVal) = NaN;
end

%-----------------
% STEP 3:  Check geographical data
%-----------------

% Check A:  Make sure it is in metres rather than in lat/lon
if ~isprop(iA, 'CellExtentInWorldX') || ~isprop(iA, 'CellExtentInWorldY')
    error('featureTrackingReadImages:badDimensionType', 'Files must use UTM coordanates not Lat/Lon (File: "%s")', fileA);
end

% Check B:  Check they're reasonable??
cX = iA.CellExtentInWorldX;
cY = iA.CellExtentInWorldY;
minVal = 0.1;
maxVal = 100;
if cX < minVal || cX > maxVal
    error('featureTrackingReadImages:xOutsideRange', 'Pixel size in X direction (%sm) is outside expected range (%s - %s m) (File: "%s")', num2str(cX), num2str(minVal), num2str(maxVal), fileA);
elseif cY < minVal || cY > maxVal
    error('featureTrackingReadImages:yOutsideRange', 'Pixel size in Y direction (%sm) is smaller than minimum expected (%s - %s m) (File: "%s")', num2str(cY), num2str(minVal), num2str(maxVal), fileA);
end

% Check C:  Check they're both really close to each other?  Make sure
% within a cm
maxDiff = 0.01;
if abs(cX - cY) > maxDiff
    error('featureTrackingReadImages:xyDifferent', 'Pixel size in X direction (%sm) and Y direction (%sm) are more than %sm different (File: "%s")', num2str(cX), num2str(cY), num2str(maxDiff), fileA);
end


%--------------------------------------------------------------------------
function i_checkInputs(fileA, fileB)
% Check the two inputs are ok...

% Must both be strings (to be filenames...)
if ~ischar(fileA) || ~ischar(fileB)
    error('featureTrackingReadImages:badInputs', 'Inputs must both be strings');
elseif exist(fileA, 'file') ~= 2
    error('featureTrackingReadImages:badInput1', 'Could not find first file ("%s")', fileA);
elseif exist(fileB, 'file') ~= 2
    error('featureTrackingReadImages:badInput2', 'Could not find second file ("%s")', fileB);
end
