function plotssp3D( FileRoot )
% plotssp.m
% Plots the sound speed profile
% Accesses Fileroot.env and FileRoot.env

global units jkpsflag

sspfil = [ FileRoot '.ssp' ];
envfil = [ FileRoot '.env' ]; % append extension

%[ pltitle, ~, SSP, ~, fid ] = read_env( envfil, 'BELLHOP' ); % read in the environmental file


%%

% Read the SSPFIL
[ Nx, Ny, Nz, Segx, Segy, Segz, cmat ] = readssp3d( 'Taiwan3D.ssp' );

% set labels in m or km
if ( strcmp( units, 'km' ) )
    xlab  = 'Range, x (km)';
    ylab  = 'Range, y (km)';
else
    xlab  = 'Range, x (m)';
    ylab  = 'Range, y (m)';
    Segx = 1000 * Segx;
    Segy = 1000 * Segy;
end

%%

% mask out land (where c= 1500) in one case
cmat( cmat == 1500 ) = 1450; %NaN;

%imagesc( rProf, SSP.z, cmat );   % imagesc produces a better PostScript file, using PostScript fonts
%pcolor( rProf, SSP.z, cmat );  ...
%   shading interp; colormap( jet );

iz = 5;
Segz( iz )

figure
imagesc( Segx, Segy, squeeze( cmat( iz, :, : ) ) )
pcolor(  Segx, Segy, squeeze( cmat( iz, :, : ) ) )
shading interp

colormap( jet )
colorbar

%colorbar( 'YDir', 'Reverse' )
%set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature

xlabel( xlab )
ylabel( ylab );

%title( deblank( pltitle ) )

% set up axis lengths for publication
if ( jkpsflag )
    set( gcf, 'Units', 'centimeters' )
    set( gca, 'ActivePositionProperty', 'Position', 'Units', 'centimeters' )
    
    set( gca, 'Position', [ 2 2 14.0  7.0 ] )
    %set( gcf, 'PaperPosition', [ 3 3 19.0 10.0 ] )
end


