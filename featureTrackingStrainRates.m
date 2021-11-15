function featureTrackingStrainRates(x, y, A, B, pu, pv, dxy, resolution, daysBetweenImages, diamondSize, downSelect, scaleFactor, maxArrowVal, strainFile, drawPictures)
% Calculate strain rates and corresponding figures...

includeTestFigures = drawPictures;

% dxy is PIXEL movement, but need metre displacement hereafter
dxy = dxy*resolution;

%--------------
% STEP -1: Cross check for final values
%--------------

if includeTestFigures

    % Using the method of "Surface velocity and mass balance of Ice Strems D
    % and E, West Antarctica", Bindschadler et al. (1996) to cross check the
    % answers
    
    % Do this before the down select below??
    [strainRatesBind, strainDirsBins] = i_calcStrainRatesBindschadler(x, y, pu, pv, dxy, daysBetweenImages, diamondSize);
    
    % Apply down select
    idxPu = mod(pu, downSelect) == 1;
    idxPv = mod(pv, downSelect) == 1;
    idx = idxPu & idxPv;
    
    titleStr = 'COMPARISON Rates Using Bindschadler et. al. (Plus MATLAB Eigenvalues)';
    i_plotStrainRates(x, y, pu(idx), pv(idx), A, strainRatesBind(idx, :), strainDirsBins(idx, :), scaleFactor, maxArrowVal, titleStr);
end

%--------------
% STEP 0: Down select data
%--------------

% Apply down select
idxPu = mod(pu, downSelect) == 1;
idxPv = mod(pv, downSelect) == 1;
idx = idxPu & idxPv;
pu = pu(idx);
pv = pv(idx);
dxy = dxy(idx, :);


%--------------
% STEP 1:  Work out diamonds
%--------------

% For each tracked point, lay a diamond on top, orientation the same as the
% "down glacier" value according to the dxy at that point...
[allDiamondX, allDiamondY] = i_layDiamonds(x, y, pu, pv, dxy, diamondSize);

if includeTestFigures
    i_checkDiamondPlot(allDiamondX, allDiamondY, x, y, A, pu, pv, dxy);
end

%--------------
% STEP 2:  Work out where they move to!
%--------------

% Start with a gridded interpolant so we can reuse if needed....
idxUse = all(~isnan(dxy), 2);
xUse = x(pu(idxUse));
yUse = y(pv(idxUse));
giDx = scatteredInterpolant(xUse(:), yUse(:), dxy(idxUse, 1), 'linear', 'none');
giDy = scatteredInterpolant(xUse(:), yUse(:), dxy(idxUse, 2), 'linear', 'none');

% Work out the displacements for the diamond
diamondDx = giDx(allDiamondX, allDiamondY);
diamondDy = giDy(allDiamondX, allDiamondY);

% And so the new diamonds...
newDiamondX = allDiamondX + diamondDx;
newDiamondY = allDiamondY + diamondDy;

% And again, a diagnostic plot if needed....
if includeTestFigures
    i_prePostDiamondPlot(allDiamondX, allDiamondY, newDiamondX, newDiamondY, x, y, A, B);
end

% Try the same with circles?!
if includeTestFigures
    [circlesPreX, circlesPreY] = i_layCircles(x, y, pu, pv, diamondSize);
    circleDx = giDx(circlesPreX, circlesPreY);
    circleDy = giDy(circlesPreX, circlesPreY);
    circlePostX = circlesPreX + circleDx;
    circlePostY = circlesPreY + circleDy;
    i_prePostCirclePlot(circlesPreX, circlesPreY, circlePostX, circlePostY, x, y, A, B);
end

%--------------
% STEP 3:  Possibly filter out weird diamonds here!!!
%--------------

% E.g. where the centre point isn't in the middle any more, the orientation
% is weird, etc...

% TODO - should probably filter the displacements first then see what is
% needed here...

%--------------
% STEP 4:  Work out strain rates
%--------------
[strainRates, strainDirs] = i_calcStrainRatesNye(allDiamondX, allDiamondY, newDiamondX, newDiamondY, daysBetweenImages, dxy);

%--------------
% STEP 5:  And some kind of plot...
%--------------
titleStr = 'Strain Rates Using Nye';
[plotDataX, plotDataY, strainRates, strainDirs] = i_plotStrainRates(x, y, pu, pv, A, strainRates, strainDirs, scaleFactor, maxArrowVal, titleStr, includeTestFigures);

%--------------
% STEP 6: Write the strain file
%--------------

% And partial strain files?!
[p, f, e] = fileparts(strainFile);

% Want to split into positive and negative strain arrows in each file...
%   This is cheating to bypass a proper rewrite of earlier bits that we
%   know already work...  All arrows have 24 points.  First 12 are for eps3
%   (i.e. first column of strain rates).  Last 12 are for eps1, i.e. second
%   column of strain rates, so...
idxCol1 = 1:12;
idxCol2 = 13:24;

% Now we need to split accordingly into those that are both extensive, only
% one extensive, and neither extensive...
idxExt = strainRates >= 0;
idxAllExt = idxExt(:, 1) & idxExt(:, 2);
idxCol1Ext = idxExt(:, 1) & ~idxExt(:, 2);
idxCol2Ext = idxExt(:, 2) & ~idxExt(:, 1);
idxNoneExt = ~idxExt(:, 1) & ~idxExt(:, 2);

% Preallocate
xExt = cell(size(plotDataX));
yExt = xExt;

xComp = xExt;
yComp = xExt;

% And split accordingly...

% For the "All extensive" ones, just copy everything to the extensive and
% leave compressive as is...
xExt(idxAllExt) = plotDataX(idxAllExt);
yExt(idxAllExt) = plotDataY(idxAllExt);

% For only first column extensive...  Copy the first column arrows to
% extensive and the second column to compressive...

% COULD use a cellfun here but REALLY SLOW!
for i = 1:numel(xExt)
    if idxCol1Ext(i)
        xExt{i} = plotDataX{i}(idxCol1);
        yExt{i} = plotDataY{i}(idxCol1);
        xComp{i} = plotDataX{i}(idxCol2);
        yComp{i} = plotDataY{i}(idxCol2);
    elseif idxCol2Ext(i)
        xExt{i} = plotDataX{i}(idxCol2);
        yExt{i} = plotDataY{i}(idxCol2);
        xComp{i} = plotDataX{i}(idxCol1);
        yComp{i} = plotDataY{i}(idxCol1);
    end
end

% And all compressive
xComp(idxNoneExt) = plotDataX(idxNoneExt);
yComp(idxNoneExt) = plotDataY(idxNoneExt);

idxEmptyExt = cellfun('isempty', xExt);
i_writeStrainFile(fullfile(p, [f, '_extensive', e]), xExt(~idxEmptyExt), yExt(~idxEmptyExt), strainRates(~idxEmptyExt, :), strainDirs(~idxEmptyExt, :));

idxEmptyComp = cellfun('isempty', xComp);
i_writeStrainFile(fullfile(p, [f, '_compressive', e]), xComp(~idxEmptyComp), yComp(~idxEmptyComp), strainRates(~idxEmptyComp, :), strainDirs(~idxEmptyComp, :));

i_writeStrainFile(strainFile, plotDataX, plotDataY, strainRates, strainDirs);

%--------------------------------------------------------------------------
function [fcnC1, fcnC2] = i_makeDownSelFcns(idxCol1, idxCol2)
fcnC1 = @(x)x(idxCol1);
fcnC2 = @(x)x(idxCol2);

%--------------------------------------------------------------------------
function [strainRates, strainDirs] = i_calcStrainRatesBindschadler(x, y, pu, pv, dxy, daysBetweenImages, diamondSize)
% Use the calculations in appendix A of "Surface velocity and mass balance 
% of Ice Strems D and E, West Antarctica", Bindschadler et. al. 1996 to 
% cross check the answers

% Start by reducing to an x-y grid of points?
[unPu, ~, idxInUnPu] = unique(pu);
[unPv, ~, idxInUnPv] = unique(pv);

gridX = NaN(numel(unPv), numel(unPu));
gridY = gridX;
gridDx = gridX;
gridDy = gridY;

idxInsert = sub2ind(size(gridX), idxInUnPv, idxInUnPu);

gridX(idxInsert) = x(pu);
gridY(idxInsert) = y(pv);
gridDx(idxInsert) = dxy(:, 1);
gridDy(idxInsert) = dxy(:, 2);

gridXPost = gridX + gridDx;
gridYPost = gridY + gridDy;

resolution = mean([abs(unique(diff(gridX, 1, 2))); abs(unique(diff(gridY, 1, 1)))]);

% Calculate the a, b, c, d....
%   NOTE that the b/d calculations take into acount that the image
%   rows increase in the opposite direction to the y values...
epsDot_0100 = i_calcEpsDotBind(0, 1, 0, 0, gridX, gridY, gridXPost, gridYPost, daysBetweenImages, resolution, diamondSize);
epsDot_0m100 = i_calcEpsDotBind(0, -1, 0, 0, gridX, gridY, gridXPost, gridYPost, daysBetweenImages, resolution, diamondSize);
epsDot_0110 = i_calcEpsDotBind(0, 1, 1, 0, gridX, gridY, gridXPost, gridYPost, daysBetweenImages, resolution, diamondSize);
epsDot_0m1m10 = i_calcEpsDotBind(0, -1, -1, 0, gridX, gridY, gridXPost, gridYPost, daysBetweenImages, resolution, diamondSize);
epsDot_m1000 = i_calcEpsDotBind(-1, 0, 0, 0, gridX, gridY, gridXPost, gridYPost, daysBetweenImages, resolution, diamondSize);
epsDot_1000 = i_calcEpsDotBind(1, 0, 0, 0, gridX, gridY, gridXPost, gridYPost, daysBetweenImages, resolution, diamondSize);
epsDot_m1001 = i_calcEpsDotBind(-1, 0, 0, 1, gridX, gridY, gridXPost, gridYPost, daysBetweenImages, resolution, diamondSize);
epsDot_0m110 = i_calcEpsDotBind(0, -1, 1, 0, gridX, gridY, gridXPost, gridYPost, daysBetweenImages, resolution, diamondSize);

% And so a, b, c, d
a = 0.5*(epsDot_0100 + epsDot_0m100);
b = 0.5*(epsDot_0110 + epsDot_0m1m10);
c = 0.5*(epsDot_m1000 + epsDot_1000);
d = 0.5*(epsDot_m1001 + epsDot_0m110);

% So the x and y values...
epsDotX = 0.25 * (b + d - a) + 0.75 * c;
epsDotXY = 0.5 * b - 0.5 * d;
epsDotY = 0.75 * a + 0.25 * (b + d - c);

% Eigenvalues and vectors for the symetric matrix...  Could just bung in a
% generic calculation here, but because we're trying to be independent of
% the Nye calculation (i.e. we're just using this for testing) we'll go
% with the MATLAB function...
nGrid = numel(gridX);
strainRates = NaN(nGrid, 2);
strainDirs = strainRates;
for i = 1:nGrid
    arrToUse = [epsDotX(i), epsDotXY(i); epsDotXY(i), epsDotY(i)];
    if any(~isfinite(arrToUse(:)))
        continue
    end
    [eigVecs, eigVals] = eig(arrToUse);
    strainRates(i, :) = [eigVals(1, 1), eigVals(2, 2)];
    thisDirs = atan2d(eigVecs(1, :), eigVecs(2, :));
    strainDirs(i, :) = thisDirs(:).';
end

% And lastly, UNGRID back out to the original x(pu), y(pv) size?
strainDirs = strainDirs(idxInsert, :);
strainRates = strainRates(idxInsert, :);

% 
% 
% handles.hf = figure('units', 'normalized', 'position', [0.1 0.1 0.8 0.8], 'color', 'w'); 
% 
% handles.ax = axes('box', 'on', 'parent', handles.hf, 'ydir', 'normal', 'nextplot', 'add', 'dataaspectratio', [1 1 1]);
% 
% o = 4;
% for r = 1:downSelect:size(gridX, 1)
%     for c = 1:downSelect:size(gridX, 2)
%         idxR = [r, r-o, NaN, r-o, r, NaN, r, r+o, NaN, r+o, r, NaN, r, r, NaN, r, r, NaN, r, r+o, NaN, r, r-o];
%         idxC = [c-o, c, NaN c, c+o, NaN, c+o, c, NaN, c, c-o, NaN, c-0, c, NaN, c, c+o, NaN, c, c, NaN, c, c];
%         
%         idxDrop = idxR < 1 | idxR > size(gridX, 1) | idxC < 1 | idxC > size(gridY, 2);
%         idxR(idxDrop) = NaN;
%         idxC(idxDrop) = NaN;
%         
%         idxPlot = NaN(size(idxC));
%         idxPlot(~isnan(idxC)) = sub2ind(size(gridX), idxR(~isnan(idxR)), idxC(~isnan(idxC)));
%         
%         xPlot = NaN(size(idxC));
%         yPlot = NaN(size(idxC));
%         xPlot(~isnan(idxPlot)) = gridX(idxPlot(~isnan(idxPlot)));
%         yPlot(~isnan(idxPlot)) = gridY(idxPlot(~isnan(idxPlot)));
%         
%         plot(handles.ax, xPlot, yPlot, '.r-');
%         
%         xPlot = NaN(size(idxC));
%         yPlot = NaN(size(idxC));
%         xPlot(~isnan(idxPlot)) = gridXPost(idxPlot(~isnan(idxPlot)));
%         yPlot(~isnan(idxPlot)) = gridYPost(idxPlot(~isnan(idxPlot)));
%         
%         plot(handles.ax, xPlot, yPlot, '.b-');
%         
%     end
% end
% 
% 1;

%--------------------------------------------------------------------------
function epsDot = i_calcEpsDotBind(k, l, m, n, gridX, gridY, gridXPost, gridYPost, daysBetweenImages, resolution, diamondSize)
% Calculate the "l" value from the Bindschadler paper.  Note that it
% is called "l" there, but calling it epsDot here for consistencey with Nye

% Find the Lf and Li according to spcification...
%   Start by preallocating with NaN...
xDistA = NaN(size(gridX));
yDistA = xDistA;
xDistB = xDistA;
yDistB = xDistA;

% Work out the indices of the grid to put into the calculation
[nRows, nCols] = size(gridX); 
[colsATo, colsAFrom] = i_getBindShift(k, nCols, resolution, diamondSize);
[rowsATo, rowsAFrom] = i_getBindShift(l, nRows, resolution, diamondSize);
[colsBTo, colsBFrom] = i_getBindShift(m, nCols, resolution, diamondSize);
[rowsBTo, rowsBFrom] = i_getBindShift(n, nRows, resolution, diamondSize);

% So form the values to actually use for calculation...
xDistA(rowsATo, colsATo) = gridX(rowsAFrom, colsAFrom);
yDistA(rowsATo, colsATo) = gridY(rowsAFrom, colsAFrom);

xDistB(rowsBTo, colsBTo) = gridX(rowsBFrom, colsBFrom);
yDistB(rowsBTo, colsBTo) = gridY(rowsBFrom, colsBFrom);

% And finally, the distances!!
Li = i_calcDist(xDistA, yDistA, xDistB, yDistB);

% Same again for the final...
xDistA(rowsATo, colsATo) = gridXPost(rowsAFrom, colsAFrom);
yDistA(rowsATo, colsATo) = gridYPost(rowsAFrom, colsAFrom);

xDistB(rowsBTo, colsBTo) = gridXPost(rowsBFrom, colsBFrom);
yDistB(rowsBTo, colsBTo) = gridYPost(rowsBFrom, colsBFrom);

Lf = i_calcDist(xDistA, yDistA, xDistB, yDistB);

% Convert to years...
yearsBetweenImages = daysBetweenImages / 365.25;

% And the final answer....
epsDot = (1./yearsBetweenImages) .* log(Lf./Li);

%--------------------------------------------------------------------------
function [rcTo, rcFrom] = i_getBindShift(k, nRc, resolution, diamondSize)
% For a given value of shift, find where to shift from/to

% An inflation factor so we're looking a BIT further away?!
inflation = ceil((diamondSize/2) / resolution);

% Just do the shift..
rcFrom = 1:nRc;
rcTo = rcFrom + -k * inflation;

% And drop anything that is out of range...
idxDrop = rcTo < 1 | rcTo > nRc;
rcFrom(idxDrop) = [];
rcTo(idxDrop) = [];

%--------------------------------------------------------------------------
function i_writeStrainFile(strainFile, plotDataX, plotDataY, strainRates, strainDirs)
% Write strain to shape files...

% Write arrows to shape file....
geom = repmat({'Line'}, size(plotDataX));
bbox = cellfun(@(x, y)[min(x(:)), min(y(:)); max(x(:)), max(y(:))], plotDataX, plotDataY, 'UniformOutput', false);
id = num2cell(1:numel(geom));

strainRates = num2cell(strainRates);
strainDirs = num2cell(mod(strainDirs, 180));

sh = struct('Geometry', geom(:), ...
    'BoundingBox', bbox(:), ...
    'X', plotDataX(:), ...
    'Y', plotDataY(:), ...
    'ID', id(:), ...
    'PSR3', strainRates(:, 1), ...
    'PSD3', strainDirs(:, 1), ...
    'PSR1', strainRates(:, 2), ...
    'PSD1', strainDirs(:, 2));
dbfspec = makedbfspec(sh);
shapewrite(sh, strainFile, 'DbfSpec', dbfspec);



%--------------------------------------------------------------------------
function [plotDataX, plotDataY, strainRates, strainDirs] = i_plotStrainRates(x, y, pu, pv, A, strainRates, strainDirs, scaleFactor, maxArrowVal, titleStr, includeTestFigures)
% Plot the strain rates...

% Start with the standard compare figure
% x1 = x(1);
% y1 = y(1);
% x = x - x(1);
% y = y - y(1);

% And these are the centre points to plot...
centreX = x(pu);
centreY = y(pv);

% And go through plotting arrows for each one
% maxLength = sqrt(abs(x(end) - x(1)) * abs(y(end) - y(1)) / numel(centreY));
% scaleFactor = maxLength / max(abs(strainRates(:)));

% This is a unit length arrow pointing straight up...
headLength = 0;
headWidth = 0;
arrLength = 0.25;
upE = [0, 0, -headWidth*arrLength/2, 0, headWidth*arrLength/2, NaN];
upN = [0, arrLength, arrLength*(1-headLength), arrLength, arrLength*(1-headLength), NaN];

% And straight down...
downE = [0, 0, -headWidth*arrLength/2, 0, headWidth*arrLength/2, NaN];
downN = [arrLength, 0, arrLength*headLength, 0, arrLength*headLength, NaN];

% Loop through all the arrows...
plotDataX = cell(size(strainRates));
plotDataY = plotDataX;

for i = 1:size(strainRates, 1)
    if any(abs(strainRates(i, :)) > maxArrowVal)
        continue
    end
    for j = 1:size(strainRates, 2)
        if isnan(strainRates(i, j))
            continue
        end
   
        % IF we have a contraction, change the arrow to point towards the
        % middle otherwise outwards
        if strainRates(i, j) < 0
            thisE = downE;
            thisN = downN;
        else
            thisE = upE;
            thisN = upN;
        end
        
        % Expand by the appropriate factor?
        thisE = thisE * scaleFactor * abs(strainRates(i, j));
        thisN = thisN * scaleFactor * abs(strainRates(i, j));
        
        % Rotate by the appropriate angle
        arrAngle = strainDirs(i, j);
        thisX = thisE * cosd(arrAngle) + thisN * sind(arrAngle);
        thisY = -thisE * sind(arrAngle) + thisN * cosd(arrAngle);
    
        thisX = [thisX(:); -thisX(:)];
        thisY = [thisY(:); -thisY(:)];
        
        % Shift so the centre is the right place...
        thisX = thisX + centreX(i);
        thisY = thisY + centreY(i);
        
        % And plot....
        plotDataX{i, j} = thisX;
        plotDataY{i, j} = thisY;
    end
    
end

% Concatenate each row of cells so we have one cell per shape
catFcn = @(X)cellfun(@(a, b)[a(:).', b(:).'], X(:, 1), X(:, 2), 'uniformoutput', false);
plotDataX = catFcn(plotDataX);
plotDataY = catFcn(plotDataY);

% And do the actuall plotting...
idxEmpty = cellfun('isempty', plotDataX) | cellfun('isempty', plotDataY);
plotDataX(idxEmpty) = [];
plotDataY(idxEmpty) = [];
strainRates(idxEmpty, :) = [];
strainDirs(idxEmpty, :) = [];


if includeTestFigures
    handles = featureTrackingSingleImage(x, y, A);
    set(handles.ax, 'nextplot', 'add');
    
    nStep = 10000;
    for i = 1:nStep:numel(plotDataX)
        thisIdx = i:min(i+nStep-1, numel(plotDataX));
        thisX = [plotDataX{thisIdx}];
        thisY = [plotDataY{thisIdx}];
        plot(handles.ax, thisX(:), thisY(:), '-b');
    end
    
    title(handles.ax, titleStr);
end

% Correct the outputs to the original scale though...
% plotDataX = cellfun(@(d)d + x1, plotDataX, 'UniformOutput', false);
% plotDataY = cellfun(@(d)d + y1, plotDataY, 'UniformOutput', false);

%--------------------------------------------------------------------------
function [strainRates, strainDirs] = i_calcStrainRatesNye(origX, origY, newX, newY, daysBetweenImages, dxy)
% Calculate strain rates according to Nye ("A method of determining the
% strain-rate tensor at the surface of a glacier", 1959)

%--------------
% STEP 1:  Work out all the a, b, c, d values
%--------------

% Values as per the diagrom in Nyp paper, page 410 (second page of paper)

% For the diamond before...  (Could work out with less calculation from the
% diamond size, but this is ok for now...)
[preA, preB, preC, preD] = i_calcNyeABCD(origX, origY);
[postA, postB, postC, postD] = i_calcNyeABCD(newX, newY);


%---------------
% STEP 2:  Calculate rates
%---------------

epsDot0 = i_calculateEpsDot(preA, postA, daysBetweenImages);
epsDot45 = i_calculateEpsDot(preB, postB, daysBetweenImages);
epsDot90 = i_calculateEpsDot(preC, postC, daysBetweenImages);
epsDot135 = i_calculateEpsDot(preD, postD, daysBetweenImages);

% TODO:  Add consistency check as per Nye that theoretically 
%   epsDot0 + epsDot90 = epsDot45 + epsDot135 
% so plots of the two should look pretty similar

%---------------
% STEP 3:  Calculate strain tensor
%---------------

% Again by Nye calculations...
epsDotX = -0.25 * epsDot0 + 0.25 * epsDot45 + 0.75 * epsDot90 + 0.25 * epsDot135;
epsDotZX = 0.5 * epsDot45 - 0.5 * epsDot135;
epsDotZ = 0.75 * epsDot0 + 0.25 * epsDot45 - 0.25 * epsDot90 + 0.25 * epsDot135;

%---------------
% STEP 4:  Eigenvalues/vectors
%---------------

% More Nye...
epsDot1 = 0.5*(epsDotX + epsDotZ) - sqrt(0.25*(epsDotX - epsDotZ).^2 + epsDotZX.^2);
epsDot3 = 0.5*(epsDotX + epsDotZ) + sqrt(0.25*(epsDotX - epsDotZ).^2 + epsDotZX.^2);

% Angle
phi = 0.5 * atand(2*epsDotZX ./ (epsDotX - epsDotZ));

% Note that this is the anticlockwise direction from Z, which is 90 degrees
% anticlockwise of the downhill direction...

% So translate 
downAngle = atan2d(dxy(:, 1), dxy(:, 2));
phiInGeog = mod(downAngle - 90 - phi, 360);

% And this angle corresponds to either the principle or secondary axes
idxDirIsPrinc = epsDotZ > epsDotX;

% Start with both the same
strainDirs = [phiInGeog, phiInGeog];

% For those where it corresponds to principle, shift around the minor axis
strainDirs(idxDirIsPrinc, 2) = mod(strainDirs(idxDirIsPrinc, 2) + 90, 360);

% And where it corresponds to the minor, shift the major
strainDirs(~idxDirIsPrinc, 1) = mod(strainDirs(~idxDirIsPrinc, 1) - 90, 360);

% Put together the strain rates
strainRates = [epsDot3, epsDot1];


%--------------------------------------------------------------------------
function epsDot = i_calculateEpsDot(pre, post, daysBetweenImages)
% Calcualte strain rates...

% Convert to years...
yearsBetweenImages = daysBetweenImages / 365.25;

% Start with epsdot calculations (remember we have TWO post and pre values
% here for each point....)
epsDot = (1./yearsBetweenImages) .* log(post ./ pre);

% And then take the average...
epsDot = 0.5*(epsDot(:, 1) + epsDot(:, 2));

%--------------------------------------------------------------------------
function [a, b, c, d] = i_calcNyeABCD(x, y)
% Calculate the a, b, c, d values from Nye paper.  Note that our points go
% clockwise from the most down glacier point, with the centre point last...

% Our first point to second corresponds to b2
b2 = i_calcDist(x(:, 1), y(:, 1), x(:, 2), y(:, 2));

% Second to third is d1
d1 = i_calcDist(x(:, 2), y(:, 2), x(:, 3), y(:, 3));

% Third to fourth is b1
b1 = i_calcDist(x(:, 3), y(:, 3), x(:, 4), y(:, 4));

% Fourth to third is d2
d2 = i_calcDist(x(:, 4), y(:, 4), x(:, 1), y(:, 1));

% First to centre (fifth) is c2
c2 = i_calcDist(x(:, 1), y(:, 1), x(:, 5), y(:, 5));

% Second to centre (fifth) is a2
a2 = i_calcDist(x(:, 2), y(:, 2), x(:, 5), y(:, 5));

% Third to centre (fifth) is c1
c1 = i_calcDist(x(:, 3), y(:, 3), x(:, 5), y(:, 5));

% Fourth to centre (fifth) is a1
a1 = i_calcDist(x(:, 4), y(:, 4), x(:, 5), y(:, 5));

% And combine
a = [a1, a2];
b = [b1, b2];
c = [c1, c2];
d = [d1, d2];

%--------------------------------------------------------------------------
function d = i_calcDist(x1, y1, x2, y2)
% Calculate distance between corresponding points in two vectors of [x, y]...
d = sqrt((x1-x2).^2 + (y1-y2).^2);

%--------------------------------------------------------------------------
function i_prePostDiamondPlot(origX, origY, newX, newY, x, y, A, B)
% Diagnostic plot to show where diamonds have moved to...

% Start with the standard compare figure
% origX = origX - x(1);
% origY = origY - y(1);
% newX = newX - x(1);
% newY = newY - y(1);
% x = x - x(1);
% y = y - y(1);

if false
    
    handles = featureTrackingCompareImages(x, y, A, B);
    
    % On the left, add the original diamonds...
    set([handles.ax1, handles.ax2], 'nextplot', 'add');
    i_addDiamonds(handles.ax1, origX, origY);
    i_addDiamonds(handles.ax2, newX, newY);
    
else
    handles = featureTrackingSingleImage(x, y, A);
    set(handles.ax, 'nextplot', 'add');
    i_addDiamonds(handles.ax, origX, origY, 'r', 'r');
    i_addDiamonds(handles.ax,  newX, newY, 'b', 'b');
end

%--------------------------------------------------------------------------
function i_prePostCirclePlot(origX, origY, newX, newY, x, y, A, B)
% Start with the standard compare figure
origX = origX - x(1);
origY = origY - y(1);
newX = newX - x(1);
newY = newY - y(1);
x = x - x(1);
y = y - y(1);

origX(:, end+1) = NaN;
origY(:, end+1) = NaN;
newX(:, end+1) = NaN;
newY(:, end+1) = NaN;

origX = origX.';
origY = origY.';
newX = newX.';
newY = newY.';

handles = featureTrackingCompareImages(x, y, A, B);

% On the left, add the original diamonds...
set([handles.ax1, handles.ax2], 'nextplot', 'add');

plot(handles.ax1, origX(:), origY(:), 'b');
plot(handles.ax2, newX(:), newY(:), 'b');

% Add the first column as a red dot...
plot(handles.ax1, origX(1, :), origY(1, :), '.r');
plot(handles.ax2, newX(1, :), newY(1, :), '.r');

% Similar plot but with the dot at the top in the same place on both so
% that we can check for shape changes only?
handles = featureTrackingSingleImage(x, y, A);

set(handles.ax, 'nextplot', 'add');
plot(handles.ax, origX(:), origY(:), 'r');
plot(handles.ax, origX(1, :), origY(1, :), '.r');

% Shift the x and y...
offsetX = repmat(origX(1, :) - newX(1, :), size(newX, 1), 1);
offsetY = repmat(origY(1, :) - newY(1, :), size(newY, 1), 1);
plot(handles.ax, newX(:) + offsetX(:), newY(:) + offsetY(:), 'b');
plot(handles.ax, newX(1, :)+ offsetX(1, :), newY(1, :)+offsetY(1, :), '.b');
    
%--------------------------------------------------------------------------
function i_addDiamonds(ax, x, y, dotColor, lineColor)
% Plot the diamonds...

if nargin < 4
    dotColor = 'r';
end
if nargin < 5
    lineColor = 'b';
end

% Just join up each necessary segment in turn?  Could just split with NaNs
% but risks making something a bit too big, so for now...
segmentsToJoin = {...
    [1, 2, 3, 4, 1]; ...
    [1, 5, 3]; ...
    [2, 5, 4]};

for i = 1:numel(segmentsToJoin)
    thisX = x(:, segmentsToJoin{i});
    thisY = y(:, segmentsToJoin{i});
    thisX(:, end+1) = NaN; %#ok<AGROW>
    thisY(:, end+1) = NaN; %#ok<AGROW>
    thisX = thisX.';
    thisY = thisY.';
    
    plot(ax, thisX(:), thisY(:), [lineColor, '.-']);
end

% Add the first column as a red dot...
plot(ax, x(:, 1), y(:, 1), [dotColor, '.']);

%--------------------------------------------------------------------------
function i_checkDiamondPlot(allDiamondX, allDiamondY, x, y, A, pu, pv, dxy)
% Quick plot to make sure the orientation of the diamonds do match those of
% the directions...  Can stop running this when the code is properly
% tested...

% Start with the standard compare figure
allDiamondX = allDiamondX - x(1);
allDiamondY = allDiamondY - y(1);
x = x - x(1);
y = y - y(1);

handles = featureTrackingCompareImages(x, y, A, A);

% On the left, draw the diamonds...
set([handles.ax1, handles.ax2], 'nextplot', 'add');
i_addDiamonds(handles.ax1, allDiamondX, allDiamondY);

% And the arrows on the right...
xData = x(pu);
yData = y(pv);

% Scale the arrows to be the same length as we only care about direction
% here
arrLength = sqrt(dxy(:, 1).^2 + dxy(:, 2).^2);
dxy = dxy ./ [arrLength, arrLength];
quiver(handles.ax2, xData(:), yData(:), dxy(:,1), dxy(:,2), 'b')


%--------------------------------------------------------------------------
function [allCirclesX, allCirclesY] = i_layCircles(x, y, pu, pv, diamondSize)
% Add circles around each tracked point....

% Centres
origX = x(pu);
origY = y(pv);

% A circle for correct radius...  Using 30 points...
degUse = linspace(0, 360, 30);
circX = (diamondSize / 2) * sind(degUse);
circY = (diamondSize / 2) * cosd(degUse);

allCirclesX = repmat(origX(:), 1, numel(circX)) + repmat(circX(:).', numel(origX), 1);
allCirclesY = repmat(origY(:), 1, numel(circY)) + repmat(circY(:).', numel(origY), 1);


%--------------------------------------------------------------------------
function [allDiamondX, allDiamondY] = i_layDiamonds(x, y, pu, pv, dxy, diamondSize)
% Add diamonds around each tracked point...

% Find their centers
origX = x(pu);
origY = y(pv);

% Points of a single diamond, centred on [0, 0], clockwise from top
% followed by the centre
diamondX = [0, diamondSize/2, 0, -diamondSize/2, 0];
diamondY = [diamondSize/2, 0, -diamondSize/2, 0, 0];

% Replicate out...
nPoints = numel(origX);
allUnitDiamondX = repmat(diamondX, nPoints, 1);
allUnitDiamondY = repmat(diamondY, nPoints, 1);

% Find the angle we need to rotate by...
arrAngle = atan2d(dxy(:, 1), dxy(:, 2));
arrAngle = repmat(arrAngle, 1, 5);

% And rotate...
allDiamondX = allUnitDiamondX .* cosd(arrAngle) + allUnitDiamondY .* sind(arrAngle);
allDiamondY = -allUnitDiamondX .* sind(arrAngle) + allUnitDiamondY .* cosd(arrAngle);

% Then translate...
allDiamondX = allDiamondX + repmat(origX(:), 1, 5);
allDiamondY = allDiamondY + repmat(origY(:), 1, 5);






