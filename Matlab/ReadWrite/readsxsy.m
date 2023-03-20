function [ sx, sy, Nsx, Nsy ] = readsxsy( fid )

% Read source x and y coordinates
fprintf( '\n_______________________ \n' )

%%
% x coordinates

[ sx, Nsx ] = readvector( fid );

fprintf( '\n Number of source x coordinates, NSx   = %i \n', Nsx )
fprintf( '\n Source x coordinates (km) \n' )
fprintf( '%f ', sx )
fprintf( '\n' )


%%
% y coordinates

[ sy, Nsy ] = readvector( fid );

fprintf( '\n Number of source y coordinates, NSy   = %i \n', Nsy )
fprintf( '\n Source y coordinates (km) \n' )
fprintf( '%f ', sy )
fprintf( '\n' )

sx = 1000.0 * sx;   % convert km to m
sy = 1000.0 * sy;

