function featureTrackingDiagnosticFigure(x, y, A, pu, pv, C, correlationThreshold, sig2noise, sig2noiseThreshold)
% featureTrackingDiagnosticFigure - Create signal to noise figure

% Start with the standard compare figure
x = x - x(1);
y = y - y(1);
handles = featureTrackingCompareImages(x, y, A, A);
set(handles.hf, 'colormap', jet(64));

%---------------
% Correlation
%---------------

% Add scatter plot for the Correlation values
set(handles.ax1, 'nextplot', 'add', 'clim', [0 1]);
idxShow = ~isnan(C);
scatter(handles.ax1, x(pu(idxShow)), y(pv(idxShow)), 500, C(idxShow), '.');

% And pretty up a bit
title(handles.ax1, sprintf('Correlation (Threshold Used):  %s', num2str(correlationThreshold)));
xlabel(handles.ax1, 'Metres from Left of image');
ylabel(handles.ax1, 'Metres from Bottom of image');
colorbar(handles.ax1, 'southoutside');

%---------------
% Signal2Noise
%---------------

% Add scatter plot for the Signal2noise values
set(handles.ax2, 'nextplot', 'add', 'clim', [0 2.5]);
idxShow = ~isnan(sig2noise);
scatter(handles.ax2, x(pu(idxShow)), y(pv(idxShow)), 500, sig2noise(idxShow), '.');

% And pretty up a bit
title(handles.ax2, sprintf('Signal To Noise (Threshold Used):  %s', num2str(sig2noiseThreshold)));
xlabel(handles.ax2, 'Metres from Left of image');
ylabel(handles.ax2, 'Metres from Bottom of image');
colorbar(handles.ax2, 'southoutside');
 