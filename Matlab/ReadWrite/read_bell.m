function  [ Beam ] = read_bell( fid, Bdry, freq, depthB, depthT, Rmax )

% Read the rest of the environmental file
% This is the part containing control info specific to Bellhop
% sets:
% Beam.RunType
% Beam.Nbeams
% Beam.Ibeam (optionally)
% Beam.deltas
% Beam.Box
% Beam.Box.r

global alpha Nbeams

Beam.Isingl = 0;
c0          = 1500.0;    % reference sound speed used to compute default step size
Beam.Ibwin  = 5;
Beam.rLoop  = 0;
Beam.Nimage = 3;

% *** Read run type ***

Beam.RunType = fgetl( fid );

% Extract option letters between the quotes

nchars = strfind( Beam.RunType, '''' );   % find quotes
Beam.RunType   = [ Beam.RunType( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) blanks( 5 - ( nchars( 2 ) - nchars( 1 ) ) ) ];

disp( '    ' )
fprintf( '\n_______________________ \n' )


switch ( Beam.RunType(1:1) )
    case ( 'R' )
        disp( 'Ray trace run' )
    case ( 'E' )
        disp( 'Eigenray trace run' )
    case ( 'I' )
        disp( 'Incoherent TL calculation' )
    case ( 'S' )
        disp( 'Semi-coherent TL calculation' )
    case ( 'C' )
        disp( 'Coherent TL calculation' )
    case ( 'A' )
        disp( 'Arrivals calculation' )
    case ( 'a' )
        disp( 'Arrivals calculation' )
        Beam.RunType(1:1) = 'A';   % 'a' for binary mode is irrelevant in the Matlab version
    otherwise
        fclose all;
        error( 'Fatal Error: Unknown RunType selected' )
end

switch ( Beam.RunType(2:2) )
    case ( 'C' )
        disp( 'Cartesian beams' )
    case ( 'R' )
        disp( 'Ray centered beams' )
    case ( 'S' )
        disp( 'Simple gaussian beams' )
    case ( 'B' )
        disp( 'Geometric gaussian beams' )
    otherwise
        Beam.RunType(2:2) = 'G';
        disp( 'Geometric hat beams' )
end

Nbeams = fscanf( fid, '%i', 1 );
fprintf( '\nNumber of beams = %i \n', Nbeams )
Beam.Nbeams = Nbeams;

if ( Nbeams == 0 )
    Nbeams = max( ceil( 0.3 * 1000 * Rmax * freq / c0 ), 300 );   % automatically estimate NBeams to use
    fprintf( 'Nbeams calculated automatically, Nbeams = %i \n', Nbeams )
end

if ( Bdry.Top.Opt( 6:6 ) == 'I' )
    Beam.Ibeam = fscanf( fid, '%i', 1 );
end

fgetl( fid );

fprintf( 'Beam take-off angles (degrees) \n' );
alpha = fscanf( fid, '%f', Nbeams );

fprintf( '%f ', alpha )
fprintf( '\n' )

if Nbeams > 2
    alpha = linspace( alpha( 1 ), alpha( 2 ), Nbeams )';
end

% check for full 360 degree sweep and remove duplicate beams at the ends
if ( Nbeams > 1 )
   if ( abs( mod( alpha( Nbeams ) - alpha( 1 ), 360.0 ) ) < 10.0 * eps( 1.0D0 ) )
      Nbeams = Nbeams - 1;
   end
end

Beam.alpha = alpha;

fgetl( fid );

% *** Limits for tracing beams ***

Beam.deltas = fscanf( fid, '%f', 1 );
Beam.Box.z  = fscanf( fid, '%f', 1 );
Beam.Box.r  = fscanf( fid, '%f', 1 );

fprintf( '\nStep length,       deltas = %d m \n', Beam.deltas )
fprintf( 'Maximum ray depth, zBox   = %d m \n', Beam.Box.z )
fprintf( 'Maximum ray range, rBox   = %d km\n', Beam.Box.r )

% *** Automatic step size selection < 10 wavelengths

if ( Beam.deltas == 0.0 )
    Beam.deltas = ( depthB - depthT ) / 10.0;
    fprintf( 'Default step length,     deltas = %d m \n', Beam.deltas )
end

fgetl( fid );

% ****** Beam characteristics ******

% note: curvature change can cause overflow in grazing case
% suppress by setting BeamType(2:2) = 'Z'

switch Beam.RunType(2:2)
    case {'G', 'B', 'S' }
    Beam.Type(1:2) = 'MS';
    Beam.rLoop     = 1.0;
    Beam.epmult    = 1.0;
    Beam.Type(3:3) = Beam.RunType(3:3);
    switch ( Beam.Type(3:3) )
        case ( 'S' )
            disp( 'Beam shift in effect' )
        otherwise
            disp( 'No beam shift in effect' )
    end
    otherwise
    Temp = fscanf( fid, '%s', 1 );
    % Extract option letters between the quotes
    
    nchars = strfind( Temp, '''' );   % find quotes
    Beam.Type   = [ Temp( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) blanks( 5 - ( nchars( 2 ) - nchars( 1 ) ) ) ];
    
    Beam.epmult   = fscanf( fid, '%f',  1 );
    Beam.rLoop    = fscanf( fid, '%f',  1 );
    fprintf( '\n\nType of beam = %s', Beam.Type(1:1) )
    
    switch ( Beam.Type(2:2) )
        case ( 'D' )
            disp( 'Curvature doubling invoked' )
        case ( 'Z' )
            disp( 'Curvature zeroing invoked' )
        case ( 'S' )
            disp( 'Standard curvature condition' )
        otherwise
            fclose all;
            error( 'Fatal Error: Unknown curvature condition' )
    end
    
    fprintf( 'Epsilon multiplier %d \n', Beam.epmult )
    fprintf( 'Range for choosing beam width %d \n', Beam.rLoop )
    
    %  ****** Images, windows ******
    fgetl( fid );

    Beam.Nimage = fscanf( fid, '%i', 1 );
    Beam.Ibwin  = fscanf( fid, '%i', 1 );
    disp( '  ' )
    fprintf( 'Number of images, Nimage = %i \n', Beam.Nimage )
    fprintf( 'Beam windowing parameter,  Ibwin = %i \n', Beam.Ibwin )
    
end

if ( length( Beam.RunType ) < 4 )
    Beam.RunType(4:4) = ' ';
end

switch ( Beam.RunType(4:4) )
    case ( 'R' )
        disp( 'Point source (cylindrical coordinates)' )
    case ( 'X' )
        disp( 'Line source (Cartesian coordinates)' )
    otherwise
        Beam.RunType(4:4) = 'R';
        disp( 'Point source (cylindrical coordinates)' )
end
