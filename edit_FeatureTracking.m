%% User parameters:  To be changed between each run

% ------------------ 
% Processing resolution (i.e. how many pixels get tracked)
processRes = 100; 

% "Super sampling" for ImGRAFT (unsure what this does yet so leaving as
% the default)
super = 2;

% Width and height of processing template (window) in ImGRAFT, in PIXELS
whtemplate=120; 

% Specify the maximum expected displacement (m)
maxdisplacement=80; %metres

% Signal-noise ratio threshold for deciding what to keep AFTER the ImGRAFT
% processing 
sigNoiseThreshold = 2; %2.5; 

% Peak correlation threshold for deciding what to keep AFTER the ImGRAFT
% processing 
correlationThreshold = 0.4; %0.5; 

% Name to out in plot titles
processTitle = 'Store Glacier'; 

% --------- 
% Input Image files

%NB/ Make sure the images are the same resolution and dimension 

% First file to be read in (earliest date) 
fileA = '/Volumes/Seagate/Store_Glacier/Data_WORKNG/UAV_DEMS/2014.06.13_DEM_0.5m_Clip.tif';
dateA = '13-June-2014';

% Second file to be read in (latest date) 
fileB = '/Volumes/Seagate/Store_Glacier/Data_WORKNG/UAV_DEMS/2014.06.16_DEM_0.5m_Clip.tif';
dateB = '16-June-2014';

% --------- 
% Output parameters

% Arrows in the shape files output
scaleFactor = 20;   % This is how many METRES one unit of velocity will be drawn as
downSelect = 2;  % This is how many arrows to keep from the original grid.  1 = all of them, 2 = every second, 3 = every third, etc.

% -------------
% Output files 
 
% Directory to output into 
outputDir = '/Volumes/Seagate/Store_Glacier/Data_WORKNG/Resampled_Strain/Final'; 

% Filename stem to use for all files 
%   e.g. if 'Run1' then outputs might be: 
%       'Run1_velo.tif'
%       'Run1_strain.tif' - keep a log elsewhere of what parameters were
%       used for each run.
outputFileStem = '2014-06-13_2014-06-16_100_2_120_80_2_0.4_0.5m';  

%% Read images

% Pass over to a useful helper function to read in, make consistent, and do
% some rudimentary checking that the files as loaded appear consistent with
% our expectations
[A, Ia, B, Ib] = featureTrackingReadImages(fileA, fileB);

% Number of days between the images being compared 
daysBetweenImages = datenum (dateB, 'dd-mm-yyyy') - datenum (dateA, 'dd-mm-yyyy');

%% Sharpen as per CSJ code
A = (A - min(A(:))) / (max(A(:)) - min(A(:)));
A = imsharpen (A, 'Radius', 2, 'Amount', 20); 
A(A < 0) = 0;
A(A > 1) = 1;

B = (B - min(B(:))) / (max(B(:)) - min(B(:)));
B = imsharpen (B, 'Radius', 2, 'Amount', 20); 
B(B < 0) = 0;
B(B > 1) = 1;


%% Coordinates corresponding to the image

% We checked in the read images code that they both have the same grid, so
% simply get the values from one grid or the other
resolution = Ia.CellExtentInWorldX;
x = linspace(Ia.XWorldLimits(1), Ia.XWorldLimits(2), Ia.RasterSize(2));
y = linspace(Ia.YWorldLimits(2), Ia.YWorldLimits(1), Ia.RasterSize(1));

%% Display images
% featureTrackingCompareImages(x, y, A, B);

%% Create Grid To Track
% Produce a regular grid of points to track 
%   1:25:size(B,2) means track every 25th point in the x direction 
%   1:25:size(B,1) means track every 25th point in the y direction 
[pu, pv]=meshgrid(1:processRes:size(B,2), 1:processRes:size(B,1)); 
pu = pu(:); 
pv = pv(:);

%% Carry out feature tracking 

% Specify or calculate search window size 
whsearch = round(whtemplate + 2*maxdisplacement/resolution);

% Scale Images to be between 0 and 1 
minImVal = min(min(A(:)), min(B(:)));
maxImVal = max(max(A(:)), max(B(:)));

A = (A - minImVal) / (maxImVal - minImVal);
B = (B - minImVal) / (maxImVal - minImVal);

[du, dv, C, Cnoise]=templatematch(A, B, pu(:), pv(:),... % This is IMGRAFT code - do not touch!
    'TemplateWidth', whtemplate, ... 
    'TemplateHeight', whtemplate, ...
    'SearchWidth', whsearch, ...
    'SearchHeight', whsearch, ...
    'SuperSample', super, ...
    'Initialdv', 0, ...
    'Initialdv', 0, ...
    'ShowProgress', true, ...
    'Method', 'myncc'); 

%% Apply thresholds

% Work out what to keep based on correlation threshold and signal to noise
% thresholds
signal2noise = C ./ Cnoise; 
keep = (signal2noise > sigNoiseThreshold) & (C > correlationThreshold); 

% And mark other values as NaN
du(~keep) = NaN;
dv(~keep) = NaN;

% Combine the du and dv values
dxy = [du(:), dv(:)]; 

% Need to double check this, but a POSITIVE dxy(:, 2) actually means
% movement DOWNWARDS (i.e. down through the rows of the image) but the
% orientation of our image is UPWARDS
dxy(:, 2) = -dxy(:, 2);

% Calculate horizontal displacement in m/day 
V = (dxy * resolution) / daysBetweenImages; 

% Calculate total displacement 
Vn = sqrt(sum(V.^2, 2)); 

%% Produce signal to noise and correlation figures
featureTrackingDiagnosticFigure(x, y, A, pu, pv, C, correlationThreshold, signal2noise, sigNoiseThreshold);

%% Produce velocity figure
titleStr = [processTitle, 'velocity', dateA, 'to', dateB];
cLim = [0 20];
featureTrackingVelocityFigure(x, y, A, pu, pv, V, Vn, titleStr, cLim);

%% GUI for finding where points moved to
% featureTrackingDisplacementFigure(x, y, A, B, pu, pv, dxy, resolution);

%% Produce Tiff file for the velocity

% Produce UTM22 coordinates of the points that have been feature tracked 
idxUse = ~isnan(Vn);
xy = pix2map(Ia, pv(idxUse), pu(idxUse));

% Interpolate scattered point data to orginal grid 
[xGrid, yGrid] = meshgrid(x, y);
grid = griddata((xy(:,1)), (xy(:,2)), Vn(keep), xGrid, yGrid, 'natural');

% Write velocity data to tif file 
veloFile = fullfile(outputDir, [outputFileStem, '_velo.tif']); 

info = geotiffinfo(fileA);
geotiffwrite(veloFile, single(grid), Ia, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);

%% Write shape file for arrows
 arrowFile = fullfile(outputDir, [outputFileStem, '_arrows.shp']);

% Do the down selection...

% Start by dropping any NaNs
 idxUse = ~isnan(Vn);
 uvA=[pu(:) pv(:)]; 
 imPixUsed = uvA(idxUse, :);
 vUsed = V(idxUse, :);

% Then keep only those according to the downSelect parameter
 idxForArrows = all(mod(imPixUsed, downSelect*processRes)==1, 2);

% And write to shape file...
 writeArrowsToShp(arrowFile, xy(idxForArrows, 1),xy(idxForArrows, 2),vUsed(idxForArrows, 1), vUsed(idxForArrows,2), scaleFactor);

 
%% Calculate strain rates
diamondSize = 100;
maxArrowVal = 0.01 * 365.25;
strainFile = fullfile(outputDir, [outputFileStem, '_strain.shp']);
strainScaleFactor = 250;
strainDownSelect = 2;

featureTrackingStrainRates(x, y, A, B, pu, pv, dxy, resolution, daysBetweenImages, diamondSize, strainDownSelect, strainScaleFactor, maxArrowVal, strainFile, false)

