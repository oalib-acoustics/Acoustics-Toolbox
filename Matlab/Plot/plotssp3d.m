function plotssp3d( sspfil )
% plotssp3d.m
% Plots the sound speed profile
% Accesses sspfil.ssp

global units jkpsflag
jkpsflag = 0;

iz = 5;   % select the depth-index to plot

[ ~, ~, ext ] = fileparts( sspfil );

if ( ~strcmp( ext, '.ssp' ) )
    sspfil = [ sspfil '.ssp' ]; % append extension
end

%%

% Read the SSPFIL
[ Nx, Ny, Nz, Segx, Segy, Segz, cmat ] = readssp3d( sspfil );

% set labels in m or km
if ( strcmp( units, 'km' ) )
    xlab = 'Range, x (km)';
    ylab = 'Range, y (km)';
else
    xlab = 'Range, x (m)';
    ylab = 'Range, y (m)';
    Segx = 1000 * Segx;
    Segy = 1000 * Segy;
end

%%

% mask out land (where c = 1500) in one case
cmat( cmat == 1500 ) = 1450; %NaN;

%imagesc( rProf, SSP.z, cmat );   % imagesc produces a better PostScript file, using PostScript fonts
%pcolor( rProf, SSP.z, cmat );  ...
%   shading interp; colormap( jet );

figure
imagesc( Segx, Segy, squeeze( cmat( iz, :, : ) ) )
pcolor(  Segx, Segy, squeeze( cmat( iz, :, : ) ) )
shading interp

colormap( jet )
c = colorbar;
c.Label.String = 'Sound speed (m/s)';

%colorbar( 'YDir', 'Reverse' )
%set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature

xlabel( xlab )
ylabel( ylab );
title( [ 'z = ' num2str( Segz( iz ) ) ' m' ] )

%title( deblank( pltitle ) )

% set up axis lengths for publication
if ( jkpsflag )
    set( gcf, 'Units', 'centimeters' )
    set( gca, 'ActivePositionProperty', 'Position', 'Units', 'centimeters' )
    
    set( gca, 'Position', [ 2 2 14.0  7.0 ] )
    %set( gcf, 'PaperPosition', [ 3 3 19.0 10.0 ] )
end