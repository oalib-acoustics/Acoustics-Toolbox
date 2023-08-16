function [ TitleEnv, Opt, Comp, MLimit, NProf, rProf, Pos ] = read_flp( fileroot )

% Read the field.flp file
% Usage:
%    read_flp
% mbp 4/09

fid = fopen( [ fileroot '.flp' ], 'r' );    % open the field parameters file

if ( fid == -1 )
   fileroot
   error( 'flp file does not exist' )
end

% read Title
TitleEnv = fgetl( fid );
% Extract letters between the quotes
nchars = strfind( TitleEnv, '''' );   % find quotes
if ( ~isempty( nchars ) )
   TitleEnv   = TitleEnv( nchars( 1 ) + 1 : nchars( 2 ) - 1 );
end
disp( TitleEnv )

% read Opt
Opt = fgetl( fid );
% Extract letters between the quotes
nchars = strfind( Opt, '''' );   % find quotes
Opt    = Opt( nchars( 1 ) + 1 : nchars( 2 ) - 1 );
disp( Opt )

if ( length( Opt ) <= 2 )
   Opt( 3 : 3 ) = 'O';   % default beampattern is omni
end

if ( length( Opt ) <= 3 )
   Opt( 4 : 4 ) = 'C';   % default is Coherent not Incoherent TL
end

% select the component
if ( length( Opt ) >= 3 )
   Comp = Opt( 3 : 3 );
else
   Comp = 'P';
end

% read MLimit
MLimit   = fscanf( fid, '%i', 1 );
fprintf( 'MLimit = %i \n\n', MLimit )
fgetl( fid );

% read profile info

fprintf( '\n_______________________ \n' )

[ rProf, NProf ] = readvector( fid );

fprintf( '\n Number of profiles, NProf = %i \n', NProf )
fprintf( '\n Profile ranges, rProf (km) \n' )

if ( NProf < 10 )
   fprintf( '%8.2f  \n', rProf )   % print all the depths
else
   fprintf( '%8.2f ... %8.2f \n', rProf( 1 ), rProf( end ) ) % print first, last depth
end

% read receiver ranges
Pos.r.r = readr( fid );
Pos.r.r = 1000.0 * Pos.r.r; % convert km to m

PosTemp = readszrz( fid );
Pos.s.z = PosTemp.s.z;
Pos.r.z = PosTemp.r.z;

% read receiver range offsets (array tilt)

fprintf( '\n_______________________ \n' )

[ Pos.r.ro, Pos.Nro ] = readvector( fid );

fprintf( '\n Number of receiver range offsets (array tilt) = %i \n', Pos.Nro )
fprintf( '\n Receiver range offsets, Rro (m)\n' );

if ( Pos.Nro < 10 )
   fprintf( '%8.2f  \n', Pos.r.ro )   % print all the range offsets
else
   fprintf( '%8.2f ... %8.2f \n', Pos.r.ro( 1 ), Pos.r.ro( end ) ) % print first, last depth
end

disp( '  ' )

fprintf( '\n' )

if ( max( abs( Pos.r.ro ) ) > 0.0 )
   error( 'The routine field.m does not implement receiver range offsets' )
end
