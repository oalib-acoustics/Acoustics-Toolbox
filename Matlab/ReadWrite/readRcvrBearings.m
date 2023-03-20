function [ theta ] = readRcvrBearings( fid )

% Read receiver bearings

Ntheta = fscanf( fid, '%i', 1 );
fprintf( '\nNumber of receiver bearings = %i \n', Ntheta )
fgetl( fid );

theta = fscanf( fid, '%f', Ntheta );

fprintf( '\nReceiver bearings (degrees) \n' )
fprintf( '%f ', theta )
fprintf( '\n' )

if Ntheta > 2
   theta = linspace( theta( 1 ), theta( 2 ), Ntheta )'; % generate vector of receiver ranges
end

% full 360-degree sweep? remove duplicate angle
if ( theta( Ntheta ) == theta( 1 ) + 360.0D0 )
   Ntheta = Ntheta - 1;
end

fgetl( fid );
