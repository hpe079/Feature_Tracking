function featureTrackingDisplacementFigure(x, y, A, B, pu, pv, dxy, resolution)
% Figure showing displacements between images

% Start with the two images...
x = x - x(1);
y = y - y(1);
handles = featureTrackingCompareImages(x, y, A, B);
set(handles.hf, 'colormap', jet(64));
title(handles.ax1, 'Before');
title(handles.ax2, 'After');

% Just draw all points as dots on both sides (same colour)
set([handles.ax1, handles.ax2], 'nextplot', 'add');

% Original positions
x1 = x(pu);
y1 = y(pv);
lineA = plot(handles.ax1, x1, y1, '.b', 'markersize', 10);

% New positions
x2 = x1(:) + dxy(:, 1) * resolution;
y2 = y1(:) + dxy(:, 2) * resolution;
lineB = plot(handles.ax2, x2, y2, '.b', 'markersize', 10);

% Add points for selection and callbacks...
selA = line('xdata', [], 'ydata', [], 'parent', handles.ax1, 'linestyle', 'none', 'marker', '.', 'color', 'r', 'markersize', 10);
selB = line('xdata', [], 'ydata', [], 'parent', handles.ax2, 'linestyle', 'none', 'marker', '.', 'color', 'r', 'markersize', 10);

% And some text...
hText = uicontrol('units', 'normalized', 'position', [0 0 1 0.075], ...
    'style', 'text', 'backgroundcolor', 'w', 'string', 'Click On Point To Find Where It Moved To...!', 'fontsize', 12);

set(lineA, 'buttondownfcn', @(h, edata)i_updateSelected(handles.ax1, 1, x1, y1, x2, y2, selA, selB, hText));
set(lineB, 'buttondownfcn', @(h, edata)i_updateSelected(handles.ax2, 2, x1, y1, x2, y2, selA, selB, hText));

%--------------------------------------------------------------------------
function i_updateSelected(hAx, selAx, x1, y1, x2, y2, selA, selB, hText)
% Update the selected point when clicked on...

% Where was clicked?
selPoint = get(hAx, 'currentpoint');
selX = selPoint(1, 1);
selY = selPoint(1, 2);

% Where was it closest to?!
if selAx == 1
    % Clicked on first axes, so look for closest point in x1 and y1
    [~, idxSel] = min((selX - x1).^2 + (selY - y1).^2);
else
    % Clicked on second axes so look for closest point in x2 and y2
    [~, idxSel] = min((selX - x2).^2 + (selY - y2).^2);
end

% And set the data of the selected points appropriately
set(selA, 'xdata', x1(idxSel), 'ydata', y1(idxSel));
set(selB, 'xdata', x2(idxSel), 'ydata', y2(idxSel));

if isnan(x2(idxSel)) || isnan(y2(idxSel))
    textStr = 'Point Failed To Track Within Specified Thresholds';
else
    dX = x2(idxSel) - x1(idxSel);
    dY = y2(idxSel) - y1(idxSel);
    
    if dX < 0 
        strHorz = sprintf('%.2fm Left', -dX);
    else
        strHorz = sprintf('%.2fm Right', dX);
    end
    if dY < 0
        strVert = sprintf('%.2fm Down', -dY);
    else
        strVert = sprintf('%.2fm Up', dY);
    end
    
    textStr = sprintf('Movement:  %s; %s', strHorz, strVert);
end
set(hText, 'string', textStr);
