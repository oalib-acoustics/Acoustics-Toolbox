 function caxisrev( clim )

% sets the color axis limits and reverses the ticks
%
% usage: caxisrev( clim )
%
% MBP January 2010

caxis( clim )
colorbar( 'YDir', 'Reverse' )
drawnow

