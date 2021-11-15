function handles = featureTrackingSingleImage(x, y, A)
% Create a single zoomable image where the colour limits change when zoomed
% in.  No colourbar...?  Image drawn as [R, G, B] so independent of
% colormap which can then be used for other things...

aRange = [min(A(:)), max(A(:))];

% And draw the images
handles.hf = figure('units', 'normalized', 'position', [0.1 0.1 0.8 0.8], 'color', 'w'); 

handles.ax = axes('box', 'on', 'parent', handles.hf, 'ydir', 'normal', 'nextplot', 'add', 'dataaspectratio', [1 1 1]);

% Work out the actual image...
imData = i_calcImData(A, aRange);
handles.im = image(x, y, imData, 'parent', handles.ax);
axis(handles.ax, 'tight')
% colorbar('peer', handles.ax, 'location', 'southoutside');

z = zoom(handles.hf);
set(z, 'ActionPostCallback', @(h, edata)i_postZoom(x, y, A, handles));

p = pan(handles.hf);
set(p, 'ActionPostCallback', @(h, edata)i_postZoom(x, y, A, handles));

%--------------------------------------------------------------------------
function imData = i_calcImData(A, aRange)
% Create RGB data from the original image and range...
imData = (A - aRange(1)) / (aRange(2) - aRange(1));
imData = repmat(imData, [1 1 3]);

%--------------------------------------------------------------------------
function i_postZoom(x, y, A, handles)

% Find the new range of values across the image...
xLim = get(handles.ax, 'xlim');
yLim = get(handles.ax, 'ylim');
xInRange = x >= xLim(1) & x <= xLim(2);
yInRange = y >= yLim(1) & y <= yLim(2);

if ~any(xInRange) || ~any(yInRange)
    aRange = [0 1];
else
    aInRange = A(yInRange, xInRange);
    aRange = [min(aInRange(:)), max(aInRange(:))];
    if aRange(1) == aRange(2)
        aRange = aRange + [-0.5, 0.5];
    end
end

imData = i_calcImData(A, aRange);
set(handles.im, 'cdata', imData);
