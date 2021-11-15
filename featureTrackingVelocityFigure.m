function featureTrackingVelocityFigure(x, y, A, pu, pv, V, Vn, titleStr, cLim)
% Velocity figure for feature tracking...

% Create a standard single figure...
x = x - x(1);
y = y - y(1);
handles = featureTrackingSingleImage(x, y, A);
set(handles.ax, 'clim', cLim);
set(handles.hf, 'colormap', jet(64));

% Produce scatter plot where colours represent speed 
set(handles.ax, 'nextplot', 'add');
idxShow = ~isnan(Vn);
scatter(handles.ax, x(pu(idxShow)), y(pv(idxShow)), 500, Vn(idxShow), '.') 

% Produce arrows which show direction 
xData = x(pu);
yData = y(pv);
quiver(handles.ax, xData(:), yData(:), V(:,1), V(:,2),'k')

% Specify the size of scalebarned from using the DEM 
colorbar(handles.ax, 'westoutside'); 
title(handles.ax, titleStr);
