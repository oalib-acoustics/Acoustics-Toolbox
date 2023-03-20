function [ Title, Opt, MLimit, R, Rmin, Rmax, Pos, sx, sy, Nsx, Nsy, theta, Ntheta, NNodes, NElts, x, y, ModeFileName, Node ] = ...
   read_flp3d( fileroot )

% Read the field3d.flp file
% Usage:
%    read_flp3d
% mbp November 2012

MaxM = 4000;

fid = fopen( [ fileroot '.flp' ], 'r' );    % open the field parameters file

if ( fid == -1 )
   error( [ fileroot '.flp file does not exist' ] )
end

% read Title
Title = fgetl( fid );
% Extract letters between the quotes
nchars = strfind( Title, '''' );   % find quotes
if ( ~isempty( nchars ) )
   Title   = Title( nchars( 1 ) + 1 : nchars( 2 ) - 1 );
end
disp( Title )

% read Opt
Opt = fgetl( fid );
% Extract letters between the quotes
nchars = strfind( Opt, '''' );   % find quotes
Opt    = Opt( nchars( 1 ) + 1 : nchars( 2 ) - 1 );
disp( Opt )

if ( length( Opt ) <= 2 )
   Opt( 3 : 3 ) = 'O';   % default beampattern is omni
end

% read MLimit
MLimit   = fscanf( fid, '%i', 1 );
fprintf( 'Max. number of modes, MLimit = %i \n\n', MLimit )
fgetl( fid );

MLimit = min( MLimit, MaxM );

% Read source/receiver information
[ sx, sy, Nsx, Nsy ] = readsxsy( fid );          % Read source x-y coordinates
Pos = readszrz( fid );    % Read source/rcvr depths

R = readr( fid );   % read receiver ranges

R    = 1000.0 * R; % convert km to m
Rmin = R(  1  );
Rmax = R( end );

theta = readRcvrBearings( fid );  % Read angles for radials
Ntheta = length( theta );

% Read nodal coordinates
NNodes = fscanf( fid, '%i', 1 );
fprintf( '\nNNodes = %i \n\n', NNodes )
fgetl( fid );

x    = zeros( NNodes, 1 );
y    = zeros( NNodes, 1 );
Iset = zeros( NNodes, 1 );
ModeFileName = cell( NNodes, 1 );

for I = 1 : NNodes
   %[ x( I ), y( I ), foo ] = fscanf( fid, '%f %f %s', 1 )
   x( I ) = fscanf( fid, '%f', 1 );
   y( I ) = fscanf( fid, '%f', 1 );
   temp   = fscanf( fid, '%s', 1 );
   
   % Extract letters between the quotes
   nchars = strfind( temp, '''' );   % find quotes
   temp   = temp( nchars( 1 ) + 1 : nchars( 2 ) - 1 );
   ModeFileName{ I } = cellstr( temp );
   
   fgetl( fid );
end

x = 1000.0 * x;   % convert km to m
y = 1000.0 * y;

% Read in element definitions
NElts = fscanf( fid, '%i', 1 );
fprintf( 'NElts = %i \n\n', NElts )
fgetl( fid );

Node            = zeros( 3, NElts );

for IElt = 1 : NElts
   Node( :, IElt ) = fscanf( fid, '%i', 3 );
   fgetl( fid );
end

