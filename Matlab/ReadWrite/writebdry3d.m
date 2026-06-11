function writebdry3d( bdryfil, interp_type, Bdry )

% Write a 3d bathymetry file for use by BELLHOP3D
%
% useage:
% writebdry3d( bdryfil, interp_type, Bdry )
%    Bdry.X has the x coordinates (km)
%    Bdry.Y has the y coordinates (km)
%    Brdr.z has the z coordinates (m) arranged in a ny x nx matrix
%
% mbp July 2011

% Some of my code uses Bathy.Lonkm ad .Latkm but I think this is confusing
% since they are Eastings and Northings at that point, not lat/longs.

%Bathy.X = Bathy.Lonkm;
%Bathy.Y = Bathy.Latkm;

switch ( interp_type( 1 : 1 ) )
    case ( 'R' )
        %    disp( 'Piecewise-linear approximation to boundary' )
    case ( 'C' )
        %    disp( 'Curvilinear approximation to boundary' )
    otherwise
        fclose all;
        disp( interp_type )
        error( 'Fatal error: Unknown option for boundary/interpolation type' )
end

Bdry.depth( isnan( Bdry.depth ) ) = 0.0;   % remove NaNs

nx = length( Bdry.X );
ny = length( Bdry.Y );

fid = fopen( bdryfil, 'wt' );
fprintf( fid, '''%s'' \n', interp_type );

fprintf( fid, '%i \n', nx );
%fprintf( fid, '%f %f /', Bathy.X( 1 ), Bathy.X( end ) );

fprintf( fid, '%f ', Bdry.X( 1 : end ) );
fprintf( fid, '\n');

fprintf( fid, '%i \n', ny );
%fprintf( fid, '%f %f /', Bathy.Y( 1 ), Bathy.Y( end ) );
fprintf( fid, '%f ', Bdry.Y( 1 : end ) );
fprintf( fid, '\n');

for iy = 1 : ny
    fprintf( fid, '%9.3f ', Bdry.depth( iy, : ) );
    fprintf( fid, '\n');
end

% if 'long' format append a matrix with province types

if ( length( interp_type ) > 1 )
    if ( interp_type( 2 : 2 ) == 'L' )
        for iy = 1 : ny
            fprintf( fid, '%3i ', Bdry.province( iy, : ) );
            fprintf( fid, '\n');
        end

        NProvinces = max( max( Bdry.province ) );
        fprintf( fid, '%i \n', NProvinces );
        for iProv = 1 : NProvinces
            fprintf( fid, '%f ', Bdry.geotype( iProv, : ) );
            fprintf( fid, '\n' );
        end
    end
end

fclose( fid );
