function handles = featureTrackingCompareImages(x, y, A, B)
% Create two linked zoomable images where the colour limits change when
% zoomed in.  No colorbar.  Image drawn as [R, G, B] so independent
% of colormap which can then be used for other things...

aRange = [min(A(:)), max(A(:))];
bRange = [min(B(:)), max(B(:))];

% And draw the images
handles.hf = figure('units', 'normalized', 'position', [0.1 0.1 0.8 0.8], 'color', 'w'); 

[handles.ax1, handles.im1] = i_createSingleImage(x, y, A, aRange, 1, handles.hf);
[handles.ax2, handles.im2] = i_createSingleImage(x, y, B, bRange, 2, handles.hf);

% Link x and y limits
linkaxes ([handles.ax1 handles.ax2]) 

% And set up post pan/zoom
z = zoom(handles.hf);
set(z, 'ActionPostCallback', @(h, edata)i_postZoom(x, y, A, B, handles));

p = pan(handles.hf);
set(p, 'ActionPostCallback', @(h, edata)i_postZoom(x, y, A, B, handles));

%--------------------------------------------------------------------------
function imData = i_calcImData(A, aRange)
% Create RGB data from the original image and range...
imData = (A - aRange(1)) / (aRange(2) - aRange(1));
imData = repmat(imData, [1 1 3]);

%--------------------------------------------------------------------------
function [ax, im] = i_createSingleImage(x, y, A, aRange, idxSubplot, hf)
% Create a single axes/image

% Start with the axes
ax = subplot (1, 2, idxSubplot, 'box', 'on', 'parent', hf, 'ydir', 'normal', 'nextplot', 'add', 'dataaspectratio', [1 1 1]);

% Add the image
imDataA = i_calcImData(A, aRange);
im = image(x, y, imDataA, 'parent', ax);
axis(ax, 'tight')

%--------------------------------------------------------------------------
function i_postZoom(x, y, A, B, handles)

% Find the new range of values across the image...
xLim = get(handles.ax1, 'xlim');
yLim = get(handles.ax1, 'ylim');
xInRange = x >= xLim(1) & x <= xLim(2);
yInRange = y >= yLim(1) & y <= yLim(2);

if ~any(xInRange) || ~any(yInRange)
    aRange = [0 1];
    bRange = [0 1];
else
    aInRange = A(yInRange, xInRange);
    bInRange = B(yInRange, xInRange);
    aRange = [min(aInRange(:)), max(aInRange(:))];
    bRange = [min(bInRange(:)), max(bInRange(:))];
    
    if aRange(1) == aRange(2)
        aRange = aRange + [-0.5, 0.5];
    end
    if bRange(1) == bRange(2)
        bRange = bRange + [-0.5, 0.5];
    end        
    
end

imDataA = i_calcImData(A, aRange);
imDataB = i_calcImData(B, bRange);

set(handles.im1, 'cdata', imDataA);
set(handles.im2, 'cdata', imDataB);
