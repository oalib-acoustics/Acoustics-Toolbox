function plotshd3d( filename )

% this doesn't work ...
% contourslice needs x, y, z data to be plaid
% saved as a possible starting point for a corrected version

% plot a single TL surface in dB (polar coordinates)
% usage: plotshd3D( filename )

% open the file and read data

[ PlotTitle, ~, ~, ~, Pos, pressure ] = read_shd( filename );

pressure = squeeze( pressure( :, 1, :, : ) );   % take first source depth

theta = Pos.theta;
rkm   = Pos.r.r / 1000.0;

% shift coordinate system to arbitrary position
xs = 333;	% -12.0
ys = 315;	%  13.8

% tlmin = 75;
% tlmax = 125;

% make plot polar

[ th, r, z ] = meshgrid( theta, rkm, Pos.r.z );
th           = ( 2 * pi / 360. ) * th;   % convert to radians
[ x, y, z ]  = pol2cart( th, r, z );

x = x + xs * ones( size( x ) );
y = y + ys * ones( size( x ) );

tlt = -20.0 * log10( abs( pressure ) );
%tlt = tlt';
tlt = permute( tlt, [ 3 1 2 ] );
% *** plot ***

tej = flipud( colormap( jet( 256 ) ) );

whos

contourslice( x, y, z, tlt, [], [], [ 1000 : 1000 : 5000 ] )
surfc( x, y, tlt ); shading interp
colormap( tej );
colorbar
% caxisrev( [ tlmin, tlmax ] )
%view( 2 )
xlabel( 'Range (km)' )
ylabel( 'Range (km)' )
zlabel( 'Depth (m)' )
title( deblank( PlotTitle ) )
%axis( 'equal' )
